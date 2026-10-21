#!/usr/bin/env python
"""
计算 GO term 的信息量（information content），以及每个蛋白注释的总信息量。

两种常见定义都算：
  IC(t) = -log2( n(t) / n(root) )
        n(t)    = 被注释到 t 的蛋白数（注释已沿 is_a/part_of 传播到所有祖先）
        n(root) = 该分支根节点的蛋白数（即在 F/P/C 里有任何注释的蛋白数）
        含义：在所有有注释的蛋白里，随机一个蛋白有 t 的概率有多小

  IA(t) = -log2( n(t) / n(所有父节点同时存在) )        —— CAFA 用的 information accretion
        含义：已经知道蛋白有 t 的所有父节点时，再知道它有 t 带来多少新信息

每个蛋白的注释信息量 = 它传播后所有 term 的 IA 之和（CAFA 加权评估用的就是这个量）

用法:
  python compute_ic.py experimental_go_per_protein_leaf.tsv go-basic.obo
  python compute_ic.py experimental_go_per_protein_leaf_noPB.tsv go-basic.obo --suffix _noPB
输出:
  go_term_ic{suffix}.tsv       每个 GO term 的 n_proteins / IC / IA
  protein_ic{suffix}.tsv       每个蛋白 F/P/C 各自的总 IA，以及叶子 term 里最大的 IC
"""
import sys, math, argparse
from collections import defaultdict
import pandas as pd
from remove_redundant_go import parse_obo, make_ancestor_fn   # 复用上一步写好的 obo 解析

ROOTS = {"F": "GO:0003674", "P": "GO:0008150", "C": "GO:0005575"}
ASPECT = {"molecular_function": "F", "biological_process": "P", "cellular_component": "C"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv"); ap.add_argument("obo")
    ap.add_argument("--suffix", default="", help="输出文件名后缀，例如 _noPB")
    args = ap.parse_args()
    tsv, obo, sfx = args.tsv, args.obo, args.suffix
    parents, namespace, name, alt2main, obsolete = parse_obo(obo)
    sys.setrecursionlimit(10000)
    ancestors = make_ancestor_fn(parents)

    df = pd.read_csv(tsv, sep="\t", dtype=str)

    # 1) 传播：每个蛋白的 term 集合 = 叶子 + 叶子的所有祖先
    prop = {}
    for _, r in df.iterrows():
        leaves = {alt2main.get(t, t) for t in r["go_terms"].split(";")}
        full = set(leaves)
        for t in leaves:
            full |= ancestors(t)
        prop[r["db_object_id"]] = (leaves, full)

    # 2) 每个 term 被哪些蛋白注释
    proteins_of = defaultdict(set)
    for pid, (_, full) in prop.items():
        for t in full:
            proteins_of[t].add(pid)

    # 3) 每个 term 的 IC 和 IA
    term_rows, ic, ia = [], {}, {}
    for t, prots in proteins_of.items():
        a = ASPECT[namespace[t]]
        n_t, n_root = len(prots), len(proteins_of[ROOTS[a]])
        ic[t] = -math.log2(n_t / n_root)
        ps = parents.get(t, set())
        if ps:
            both = set.intersection(*(proteins_of[p] for p in ps))   # 同时有所有父节点的蛋白
            ia[t] = -math.log2(n_t / len(both))
        else:
            ia[t] = 0.0                                               # 根节点
        term_rows.append({"go_id": t, "name": name.get(t, ""), "aspect": a,
                          "n_proteins": n_t, "IC": round(ic[t], 4) + 0.0, "IA": round(ia[t], 4) + 0.0})
    terms = pd.DataFrame(term_rows).sort_values(["aspect", "IC", "go_id"], ascending=[True, False, True])
    terms.to_csv(f"go_term_ic{sfx}.tsv", sep="\t", index=False)

    # 4) 每个蛋白的信息量
    prot_rows = []
    for _, r in df.iterrows():
        leaves, full = prop[r["db_object_id"]]
        row = {"db_object_id": r["db_object_id"], "db_object_symbol": r["db_object_symbol"]}
        for a in "FPC":
            row[f"IA_{a}"] = round(sum(ia[t] for t in full if ASPECT[namespace[t]] == a), 3)
        row["IA_total"] = round(row["IA_F"] + row["IA_P"] + row["IA_C"], 3)
        best = max(leaves, key=lambda t: ic[t])
        row["max_leaf_IC"] = round(ic[best], 3)
        row["most_specific_term"] = f"{best} {name.get(best, '')}"
        prot_rows.append(row)
    prots = pd.DataFrame(prot_rows).sort_values("IA_total", ascending=False)
    prots.to_csv(f"protein_ic{sfx}.tsv", sep="\t", index=False)

    for a in "FPC":
        print(f"{a}: 有注释的蛋白 {len(proteins_of[ROOTS[a]])}，出现过的 term {int((terms.aspect == a).sum())}")
    print("\n总信息量最高的 10 个蛋白:")
    print(prots.drop(columns="most_specific_term").head(10).to_string(index=False))
    print(f"\n已写出 go_term_ic{sfx}.tsv 和 protein_ic{sfx}.tsv")


if __name__ == "__main__":
    main()
