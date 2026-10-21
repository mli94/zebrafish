#!/usr/bin/env python
"""
把 leaf 表（一行一个蛋白，go_terms 用 ; 连接）和 go_term_ic 表（一行一个 GO term）拼成长表：
一行 = 一个蛋白的一个叶子 GO term，带上这个 term 的 frequency / IC / IA。

用法:
  python merge_leaf_ic.py experimental_go_per_protein_leaf.tsv go_term_ic.tsv experimental_leaf_with_ic.tsv
"""
import sys
import pandas as pd

leaf_tsv, ic_tsv, out = sys.argv[1:4]
leaf = pd.read_csv(leaf_tsv, sep="\t", dtype={"go_terms": str})
ic = pd.read_csv(ic_tsv, sep="\t")

# 分母：每个分支根节点上的蛋白数
roots = {"F": "GO:0003674", "P": "GO:0008150", "C": "GO:0005575"}
n_root = {a: int(ic.loc[ic.go_id == g, "n_proteins"].iloc[0]) for a, g in roots.items()}
ic["n_aspect_proteins"] = ic["aspect"].map(n_root)
ic["frequency"] = (ic["n_proteins"] / ic["n_aspect_proteins"]).round(6)

# 拆开 go_terms：一个蛋白一行 -> 一个蛋白-term 一行；保留 leaf 表原来的蛋白顺序
leaf["protein_rank"] = range(1, len(leaf) + 1)
long = (leaf[["protein_rank", "db_object_id", "db_object_symbol", "n_go_terms", "go_terms"]]
        .assign(go_id=lambda d: d["go_terms"].str.split(";"))
        .explode("go_id")
        .drop(columns="go_terms"))

merged = long.merge(ic[["go_id", "name", "aspect", "n_proteins", "n_aspect_proteins",
                        "frequency", "IC", "IA"]], on="go_id", how="left", validate="many_to_one")
assert merged["IC"].notna().all(), "有叶子 term 在 IC 表里找不到"
assert len(merged) == leaf["n_go_terms"].sum(), "行数和叶子 term 总数对不上"

merged = (merged.sort_values(["protein_rank", "IC"], ascending=[True, False])
                .drop(columns="protein_rank")
                .rename(columns={"n_go_terms": "protein_n_leaf_terms",
                                 "name": "go_name",
                                 "n_proteins": "term_n_proteins"}))
merged.to_csv(out, sep="\t", index=False)
print(f"{out}: {len(merged)} 行（{merged.db_object_id.nunique()} 个蛋白）")
print(merged.head(6).to_string(index=False))
