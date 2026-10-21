#!/usr/bin/env python
"""
按导师的定义计算 GO term 的频率和 IC，以及每个蛋白的 IC。

  frequency(t) = t 作为叶子注释出现的次数 / 所有叶子注释的总条数（F、P、C 合在一起）
  IC(t)        = -log2( frequency(t) )
  蛋白的 IC     = 该蛋白所有叶子 term 的 IC 之和

  - 输入是去冗余后的叶子表，每个 (蛋白, term) 只出现一次
  - 不沿 GO 树传播；没有作为叶子出现过的 term 不在表里

例：P1: GO1 GO2 GO3, P2: GO1 GO3 GO5  -> 总条数 6，freq(GO1)=2/6，freq(GO2)=1/6

用法:
  python compute_ic_annotfreq.py experimental_go_per_protein_leaf.tsv go-basic.obo
  python compute_ic_annotfreq.py experimental_go_per_protein_leaf_cafaPB.tsv go-basic.obo --suffix _cafaPB
输出:
  go_term_ic_annotfreq{suffix}.tsv   每个 term: 出现次数 / 总条数 / frequency / IC
  protein_ic_annotfreq{suffix}.tsv   每个蛋白: 叶子数、IC 之和（总计及 F/P/C 分开）、平均、最大
"""
import argparse, math
import pandas as pd
from remove_redundant_go import parse_obo

ASPECT = {"molecular_function": "F", "biological_process": "P", "cellular_component": "C"}


def term_freq_ic(pairs):
    """pairs: DataFrame，列 db_object_id, go_id（每行一条注释）。返回每个 term 的频率和 IC"""
    N = len(pairs)
    t = pairs.groupby("go_id").size().rename("n_annotations").reset_index()
    t["N_total_annotations"] = N
    t["frequency"] = t["n_annotations"] / N
    t["IC"] = -t["frequency"].map(math.log2) + 0.0
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("leaf_tsv"); ap.add_argument("obo")
    ap.add_argument("--suffix", default="")
    args = ap.parse_args()

    # 先用导师的例子自检
    demo = pd.DataFrame({"db_object_id": ["P1"] * 3 + ["P2"] * 3,
                         "go_id": ["GO1", "GO2", "GO3", "GO1", "GO3", "GO5"]})
    f = term_freq_ic(demo).set_index("go_id")["frequency"]
    assert abs(f["GO1"] - 2 / 6) < 1e-12 and abs(f["GO2"] - 1 / 6) < 1e-12

    _, namespace, name, _, _ = parse_obo(args.obo)
    leaf = pd.read_csv(args.leaf_tsv, sep="\t", dtype=str)
    pairs = (leaf[["db_object_id", "db_object_symbol", "go_terms"]]
             .assign(go_id=lambda d: d["go_terms"].str.split(";"))
             .explode("go_id").drop(columns="go_terms"))
    assert not pairs.duplicated(["db_object_id", "go_id"]).any()

    terms = term_freq_ic(pairs)
    terms.insert(1, "name", terms["go_id"].map(name))
    terms.insert(2, "aspect", terms["go_id"].map(lambda t: ASPECT[namespace[t]]))
    terms = terms.sort_values(["IC", "go_id"], ascending=[False, True])
    sfx = args.suffix
    terms.round({"frequency": 8, "IC": 4}).to_csv(f"go_term_ic_annotfreq{sfx}.tsv", sep="\t", index=False)

    # 每个蛋白：叶子 IC 之和
    p = pairs.merge(terms[["go_id", "aspect", "IC"]], on="go_id", how="left")
    prot = p.groupby(["db_object_id", "db_object_symbol"]).agg(
        n_leaf_terms=("IC", "size"), IC_sum=("IC", "sum"),
        IC_mean=("IC", "mean"), IC_max=("IC", "max")).reset_index()
    by_aspect = p.pivot_table(index="db_object_id", columns="aspect", values="IC",
                              aggfunc="sum", fill_value=0.0)
    for a in "FPC":
        prot[f"IC_sum_{a}"] = prot["db_object_id"].map(by_aspect[a] if a in by_aspect else {}).fillna(0.0)
    prot = prot[["db_object_id", "db_object_symbol", "n_leaf_terms", "IC_sum",
                 "IC_sum_F", "IC_sum_P", "IC_sum_C", "IC_mean", "IC_max"]]
    prot = prot.sort_values(["IC_sum", "db_object_id"], ascending=[False, True])
    prot.round(4).to_csv(f"protein_ic_annotfreq{sfx}.tsv", sep="\t", index=False)

    print(f"注释总条数 N = {len(pairs)}，出现过的 term {len(terms)} 个，蛋白 {len(prot)} 个")
    print("频率之和 =", round(terms["frequency"].sum(), 10))
    print("\n最常见的 5 个 term:")
    print(terms.tail(5).iloc[::-1].round(4).to_string(index=False))
    print("\nIC 之和最高的 10 个蛋白:")
    print(prot.head(10).round(2).to_string(index=False))
    print(f"\n已写出 go_term_ic_annotfreq{sfx}.tsv 和 protein_ic_annotfreq{sfx}.tsv")


if __name__ == "__main__":
    main()
