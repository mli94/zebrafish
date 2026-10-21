#!/usr/bin/env python
"""
按 Clark & Radivojac (2013, Bioinformatics 29:i53) 的方法计算：
  每个 GO term 的 information accretion  ia(v) = -log2 Pr(v | v 的所有父节点)
  每个蛋白注释的信息量                    i(T)  = sum_{v in T} ia(v)
其中 T 是该蛋白传播后的注释（含所有祖先），MF / BP / CC 三个本体分开计算。

步骤（与论文对应）：
  1. 从 GAF 取实验证据注释（论文第 4.1 节：EXP IDA IPI IMP IGI IEP TAS IC），去掉 NOT
  2. 旧编号(alt_id)换成主编号；沿 is_a + part_of 传播到所有祖先，得到每个蛋白的子图 T（第 3 节 consistency requirement）
  3. 对每个 term v，数 n(v) = T 里含 v 的蛋白数，n(Pa(v)) = T 里同时含 v 所有父节点的蛋白数
     Pr(v | Pa(v)) = n(v) / n(Pa(v))                         （第 3.1 节、Fig.1）
     加伪计数版： (n(v) + 1) / (n(Pa(v)) + 1)                 （第 4 节：每个 term 计数 +1，总数也 +1）
     根节点没有父节点，ia = 0
  4. ia(v) = -log2 Pr(v | Pa(v))；蛋白的 i(T) = T 里所有节点 ia 之和（Eq.1 下方）

用法:
  python paper_ic.py goa_zebrafish.gaf go-basic.obo
  python paper_ic.py goa_zebrafish.gaf go-basic.obo --codes EXP IDA IPI IMP IGI IEP --pseudocount 0 --suffix _6codes_nopc
输出:
  paper_term_ia{suffix}.tsv      每个 term 一行
  paper_protein_ic{suffix}.tsv   每个蛋白一行
"""
import argparse, math
from collections import defaultdict
import pandas as pd
from remove_redundant_go import parse_obo, make_ancestor_fn

PAPER_CODES = ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC"]
ASPECT = {"molecular_function": "F", "biological_process": "P", "cellular_component": "C"}


def read_gaf(path, codes):
    cols = list(range(17))
    df = pd.read_csv(path, sep="\t", comment="!", header=None, names=cols, dtype=str, quoting=3)
    df = df[df[6].isin(codes)]
    df = df[~df[3].str.contains("NOT", na=False)]
    return df[[1, 2, 4]].rename(columns={1: "pid", 2: "symbol", 4: "go"}).drop_duplicates()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gaf"); ap.add_argument("obo")
    ap.add_argument("--codes", nargs="+", default=PAPER_CODES)
    ap.add_argument("--pseudocount", type=float, default=1.0, help="论文用 1；设 0 即不加")
    ap.add_argument("--suffix", default="")
    args = ap.parse_args()
    pc = args.pseudocount

    parents, namespace, name, alt2main, obsolete = parse_obo(args.obo)
    ancestors = make_ancestor_fn(parents)

    ann = read_gaf(args.gaf, set(args.codes))
    ann["go"] = ann["go"].map(lambda t: alt2main.get(t, t))
    unknown = set(ann.go) - set(namespace)
    obs = set(ann.go) & obsolete
    ann = ann[~ann.go.isin(unknown | obs)]
    symbol = ann.drop_duplicates("pid").set_index("pid")["symbol"]

    # 2. 传播：T[pid] = 直接注释 + 所有祖先；同时记下每个蛋白的叶子
    direct = ann.groupby("pid")["go"].apply(set)
    T, leaves = {}, {}
    for pid, terms in direct.items():
        anc_union = set().union(*(ancestors(t) for t in terms))
        T[pid] = terms | anc_union
        leaves[pid] = terms - anc_union

    # 3. n(v)：每个 term 被哪些蛋白含有
    prots_of = defaultdict(set)
    for pid, ts in T.items():
        for t in ts:
            prots_of[t].add(pid)

    rows, ia_pc, ia_raw = [], {}, {}
    for v, pv in prots_of.items():
        pa = parents.get(v, set())
        n_v = len(pv)
        n_pa = len(set.intersection(*(prots_of[p] for p in pa))) if pa else n_v
        ia_raw[v] = -math.log2(n_v / n_pa) + 0.0 if pa else 0.0
        ia_pc[v] = -math.log2((n_v + pc) / (n_pa + pc)) + 0.0 if pa else 0.0
        rows.append({"go_id": v, "name": name[v], "aspect": ASPECT[namespace[v]],
                     "n_parents": len(pa), "n_v": n_v, "n_all_parents": n_pa,
                     "cond_prob": n_v / n_pa, "IA": ia_pc[v], "IA_no_pseudocount": ia_raw[v]})
    terms = (pd.DataFrame(rows).sort_values(["aspect", "IA", "go_id"], ascending=[True, False, True]))
    sfx = args.suffix
    terms.round({"cond_prob": 6, "IA": 4, "IA_no_pseudocount": 4}).to_csv(
        f"paper_term_ia{sfx}.tsv", sep="\t", index=False)

    # 4. 每个蛋白的 i(T)，三个本体分开
    prow = []
    for pid, ts in T.items():
        r = {"db_object_id": pid, "db_object_symbol": symbol[pid]}
        for a in "FPC":
            ta = [t for t in ts if ASPECT[namespace[t]] == a]
            r[f"i_{a}"] = sum(ia_pc[t] for t in ta)
            r[f"i_{a}_no_pseudocount"] = sum(ia_raw[t] for t in ta)
            r[f"n_leaf_{a}"] = sum(1 for t in leaves[pid] if ASPECT[namespace[t]] == a)
        r["i_total"] = r["i_F"] + r["i_P"] + r["i_C"]
        prow.append(r)
    prot = pd.DataFrame(prow)
    order = ["db_object_id", "db_object_symbol", "i_total", "i_F", "i_P", "i_C",
             "n_leaf_F", "n_leaf_P", "n_leaf_C",
             "i_F_no_pseudocount", "i_P_no_pseudocount", "i_C_no_pseudocount"]
    prot = prot[order].sort_values(["i_total", "db_object_id"], ascending=[False, True])
    prot.round(4).to_csv(f"paper_protein_ic{sfx}.tsv", sep="\t", index=False)

    # 汇总，与论文第 4.2 节、Fig.3 对比
    print(f"evidence codes: {' '.join(args.codes)} | pseudocount = {pc}")
    print(f"注释 {len(ann)} 条（去重后）, 蛋白 {len(T)} 个, term {len(terms)} 个; 跳过 obo 中没有的 {len(unknown)} 个、废弃的 {len(obs)} 个")
    print("\n每个本体（只统计在该本体有注释的蛋白）:")
    for a, full in zip("FPC", ["MFO", "BPO", "CCO"]):
        sub = prot[prot[f"n_leaf_{a}"] > 0]
        r = sub[f"i_{a}"].corr(sub[f"n_leaf_{a}"])
        print(f"  {full}: 蛋白 {len(sub)}, 平均信息量 {sub[f'i_{a}'].mean():.1f} bits (std {sub[f'i_{a}'].std():.1f}), "
              f"平均叶子数 {sub[f'n_leaf_{a}'].mean():.1f}, 信息量与叶子数 Pearson r = {r:.2f}")
    print("\n总信息量最高的 5 个蛋白:")
    print(prot.head(5)[["db_object_symbol", "i_total", "i_F", "i_P", "i_C"]].round(2).to_string(index=False))
    print(f"\n已写出 paper_term_ia{sfx}.tsv 和 paper_protein_ic{sfx}.tsv")


if __name__ == "__main__":
    main()
