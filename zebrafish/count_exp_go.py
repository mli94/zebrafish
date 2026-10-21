#!/usr/bin/env python
"""
统计 GAF 文件中指定 evidence code 的 GO 注释：
  1) 每个蛋白有多少个不同的 GO term（从大到小排序）
  2) 每个 GO term 被多少个蛋白注释（从大到小排序）

用法:
  python count_exp_go.py goa_zebrafish.gaf                 # 只要 evidence code = EXP
  python count_exp_go.py goa_zebrafish.gaf --experimental  # 所有实验证据: EXP IDA IPI IMP IGI IEP
  python count_exp_go.py goa_zebrafish.gaf --codes EXP IDA # 自定义
"""
import argparse
import pandas as pd

# GAF 2.x 的 17 列
GAF_COLS = [
    "db", "db_object_id", "db_object_symbol", "qualifier", "go_id",
    "db_reference", "evidence_code", "with_from", "aspect",
    "db_object_name", "db_object_synonym", "db_object_type", "taxon",
    "date", "assigned_by", "annotation_extension", "gene_product_form_id",
]
EXPERIMENTAL = ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gaf")
    ap.add_argument("--codes", nargs="+", default=["EXP"])
    ap.add_argument("--experimental", action="store_true",
                    help="使用全部实验证据码 " + " ".join(EXPERIMENTAL))
    ap.add_argument("--keep-not", action="store_true",
                    help="保留 qualifier 含 NOT 的注释（默认去掉）")
    ap.add_argument("--out", default="exp")
    args = ap.parse_args()
    codes = EXPERIMENTAL if args.experimental else [c.upper() for c in args.codes]

    df = pd.read_csv(args.gaf, sep="\t", comment="!", header=None,
                     names=GAF_COLS, dtype=str, quoting=3)

    sub = df[df["evidence_code"].isin(codes)]
    if not args.keep_not:
        sub = sub[~sub["qualifier"].str.contains("NOT", na=False)]

    # 同一个 (蛋白, GO term) 可能因不同文献出现多行，去重后再数
    pairs = sub[["db_object_id", "db_object_symbol", "go_id", "aspect"]].drop_duplicates()

    # 1) 每个蛋白的 GO term 数
    per_protein = (pairs.groupby(["db_object_id", "db_object_symbol"])
                   .agg(n_go_terms=("go_id", "nunique"),
                        n_F=("aspect", lambda a: (a == "F").sum()),
                        n_P=("aspect", lambda a: (a == "P").sum()),
                        n_C=("aspect", lambda a: (a == "C").sum()),
                        go_terms=("go_id", lambda g: ";".join(sorted(g))))
                   .reset_index()
                   .sort_values(["n_go_terms", "db_object_id"], ascending=[False, True]))

    # 2) 每个 GO term 的蛋白数
    per_term = (pairs.groupby(["go_id", "aspect"])
                .agg(n_proteins=("db_object_id", "nunique"))
                .reset_index()
                .sort_values(["n_proteins", "go_id"], ascending=[False, True]))

    per_protein.to_csv(f"{args.out}_go_per_protein.tsv", sep="\t", index=False)
    per_term.to_csv(f"{args.out}_proteins_per_go.tsv", sep="\t", index=False)

    print(f"evidence codes: {', '.join(codes)}")
    print(f"注释行数: {len(sub)}  去重后 (蛋白, GO) 对: {len(pairs)}")
    print(f"蛋白数: {per_protein.shape[0]}  不同 GO term 数: {per_term.shape[0]}\n")
    print("GO term 最多的前 10 个蛋白:")
    print(per_protein.drop(columns="go_terms").head(10).to_string(index=False))
    print("\n被最多蛋白注释的前 10 个 GO term:")
    print(per_term.head(10).to_string(index=False))
    print(f"\n已写出 {args.out}_go_per_protein.tsv 和 {args.out}_proteins_per_go.tsv")


if __name__ == "__main__":
    main()
