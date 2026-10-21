#!/usr/bin/env python
"""
去掉每个蛋白 GO 注释里的冗余（祖先）term，只保留最具体的叶子 term。

规则：如果蛋白同时有 A 和 B，且 A 是 B 的祖先（沿 is_a / part_of 往上能走到 A），
那么 A 是冗余的——有 B 就已经隐含了 A——删掉 A。

用法:
  python remove_redundant_go.py experimental_go_per_protein.tsv go-basic.obo
  python remove_redundant_go.py experimental_go_per_protein.tsv go-basic.obo --drop GO:0005515 --out experimental_go_per_protein_leaf_noPB.tsv
输出:
  默认 experimental_go_per_protein_leaf.tsv
  python remove_redundant_go.py experimental_go_per_protein.tsv go-basic.obo --cafa-pb --out experimental_go_per_protein_leaf_cafaPB.tsv
  --drop 的 term 在去冗余之前就删掉；删完没有任何 term 的蛋白不出现在输出里
  --cafa-pb 按 CAFA3 规则：去冗余后 MF 只剩 protein binding (GO:0005515) 的蛋白，去掉它的全部 MF 注释，
            BP/CC 保留；PB 之外还有别的 MF 注释的蛋白不动
"""
import sys
import argparse
from collections import defaultdict
import pandas as pd


def parse_obo(path):
    """读 go-basic.obo，返回 parents / namespace / name / alt_id->主 id / 废弃 term 集合"""
    parents = defaultdict(set)   # term -> 直接父节点（is_a 和 part_of）
    namespace, name, alt2main, obsolete = {}, {}, {}, set()
    cur, in_term = None, False
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("["):              # 新的一节：[Term] 或 [Typedef]
                in_term = (line == "[Term]")
                cur = None
                continue
            if not in_term or ": " not in line:
                continue
            key, val = line.split(": ", 1)
            if key == "id":
                cur = val
            elif key == "name":
                name[cur] = val
            elif key == "namespace":
                namespace[cur] = val
            elif key == "alt_id":
                alt2main[val] = cur
            elif key == "is_obsolete" and val == "true":
                obsolete.add(cur)
            elif key == "is_a":
                parents[cur].add(val.split()[0])            # "GO:xxx ! name" -> "GO:xxx"
            elif key == "relationship" and val.startswith("part_of "):
                parents[cur].add(val.split()[1])
    return parents, namespace, name, alt2main, obsolete


def make_ancestor_fn(parents):
    """返回一个函数 ancestors(t)：t 的所有祖先（不含自己），带缓存"""
    cache = {}
    def ancestors(t):
        if t in cache:
            return cache[t]
        result = set()
        for p in parents.get(t, ()):
            result.add(p)
            result |= ancestors(p)
        cache[t] = result
        return result
    return ancestors


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv")
    ap.add_argument("obo")
    ap.add_argument("--drop", nargs="*", default=[], help="去冗余之前先删掉的 GO term，例如 GO:0005515")
    ap.add_argument("--out", default=None)
    ap.add_argument("--cafa-pb", action="store_true", help="按 CAFA3 规则处理 protein binding-only 的 MF 注释")
    args = ap.parse_args()
    tsv, obo = args.tsv, args.obo
    out = args.out or tsv.replace(".tsv", "_leaf.tsv")
    drop = set(args.drop)
    n_dropped_proteins = 0
    n_pb_only = 0

    parents, namespace, name, alt2main, obsolete = parse_obo(obo)
    sys.setrecursionlimit(10000)
    ancestors = make_ancestor_fn(parents)
    aspect_of = {"molecular_function": "F", "biological_process": "P", "cellular_component": "C"}

    df = pd.read_csv(tsv, sep="\t", dtype=str)
    unknown, obs_seen = set(), set()
    rows = []
    for _, r in df.iterrows():
        terms = {alt2main.get(t, t) for t in r["go_terms"].split(";")}   # 旧编号换成主编号
        terms -= drop
        if not terms:                     # 只有被 drop 的 term，这个蛋白整个去掉
            n_dropped_proteins += 1
            continue
        unknown |= {t for t in terms if t not in namespace}
        obs_seen |= terms & obsolete

        # 所有 term 的祖先合在一起；落在这个集合里的 term 就是冗余的
        redundant = set()
        for t in terms:
            redundant |= ancestors(t)
        leaves = sorted(terms - redundant)
        removed = sorted(terms & redundant)

        if args.cafa_pb:
            mf = lambda t: namespace.get(t) == "molecular_function"
            if {t for t in leaves if mf(t)} == {"GO:0005515"}:
                # 这个蛋白退出 MF 评估：MF 的叶子和被删的 MF 祖先都去掉
                n_pb_only += 1
                leaves = [t for t in leaves if not mf(t)]
                removed = [t for t in removed if not mf(t)]
                terms = {t for t in terms if not mf(t)}
                if not leaves:
                    n_dropped_proteins += 1
                    continue

        asp = [aspect_of.get(namespace.get(t), "?") for t in leaves]
        rows.append({
            "db_object_id": r["db_object_id"],
            "db_object_symbol": r["db_object_symbol"],
            "n_go_terms": len(leaves),
            "n_F": asp.count("F"), "n_P": asp.count("P"), "n_C": asp.count("C"),
            "n_original": len(terms),
            "n_removed": len(removed),
            "go_terms": ";".join(leaves),
            "removed_terms": ";".join(removed),
        })

    res = (pd.DataFrame(rows)
           .sort_values(["n_go_terms", "db_object_id"], ascending=[False, True]))
    res.to_csv(out, sep="\t", index=False)

    if drop:
        print(f"删掉的 term: {sorted(drop)}；因此整个去掉的蛋白: {n_dropped_proteins}")
    if args.cafa_pb:
        print(f"CAFA 规则：MF 只有 protein binding 的蛋白 {n_pb_only} 个，已去掉其 MF 注释；"
              f"其中没有 BP/CC、整个去掉的蛋白 {n_dropped_proteins} 个")
    print(f"蛋白数: {len(res)}")
    print(f"GO term 总数（蛋白-term 对）: 去冗余前 {res.n_original.sum()}  ->  去冗余后 {res.n_go_terms.sum()}"
          f"  (删掉 {res.n_removed.sum()})")
    print(f"obo 里找不到的 term: {len(unknown)} {sorted(unknown)[:10]}")
    print(f"已废弃(obsolete)的 term: {len(obs_seen)} {sorted(obs_seen)[:10]}")
    print("\n去冗余后 GO term 最多的前 10 个蛋白:")
    print(res.drop(columns=["go_terms", "removed_terms"]).head(10).to_string(index=False))
    print(f"\n已写出 {out}")


if __name__ == "__main__":
    main()
