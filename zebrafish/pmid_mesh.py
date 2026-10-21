#!/usr/bin/env python3
"""
从 GAF 文件里提取所有 PMID，用 Biopython (Bio.Entrez) 去 PubMed 取每篇文献的 MeSH terms，
据此判断文献研究的物种（斑马鱼 / 人类 / 其他）。

注意：GAF 第 13 列 (taxon) 是“被注释蛋白”的物种，goa_zebrafish.gaf 里全部是 taxon:7955。
这里要回答的是“支持这条注释的文献实际研究的是什么物种”，所以要看文献的 MeSH。

用法：
    python pmid_mesh.py goa_zebrafish.gaf --email you@example.com
    python pmid_mesh.py goa_zebrafish.gaf --email you@example.com --experimental   # 只看实验证据码
    python pmid_mesh.py goa_zebrafish.gaf --email you@example.com --limit 20       # 先试 20 篇

输出：
    <out>_pmid_mesh.tsv        每个 PMID 一行：物种分类、物种类 MeSH、全部 MeSH、标题
    <out>_annotation_mesh.tsv  每条带 PMID 的注释一行，附上该文献的物种分类
    <out>_cache.json           已下载结果的缓存，中断后重跑会跳过已下载的 PMID
"""
import argparse
import json
import os
import sys
import time
from collections import Counter, defaultdict

import pandas as pd
from Bio import Entrez

EXP_CODES = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
             "HTP", "HDA", "HMP", "HGI", "HEP"}

# MeSH 里常见的物种类 descriptor（Organisms 分支 B01 的常用项）
ORGANISM_MESH = {
    "Humans", "Zebrafish", "Animals", "Mice", "Rats", "Xenopus laevis", "Xenopus",
    "Chick Embryo", "Chickens", "Drosophila melanogaster", "Drosophila",
    "Caenorhabditis elegans", "Saccharomyces cerevisiae", "Oryzias", "Cattle",
    "Swine", "Dogs", "Rabbits", "Cricetinae", "Cricetulus", "Macaca mulatta",
    "Escherichia coli", "Arabidopsis",
}


def read_gaf(path):
    cols = ["db", "db_object_id", "symbol", "qualifier", "go_id", "reference",
            "evidence", "with_from", "aspect", "name", "synonym", "type",
            "taxon", "date", "assigned_by", "ext", "gp_form"]
    df = pd.read_csv(path, sep="\t", comment="!", header=None, names=cols,
                     dtype=str, quoting=3)
    # reference 列可能是 "PMID:123|GO_REF:..." 这种多值
    df["pmid"] = df["reference"].str.split("|")
    df = df.explode("pmid")
    df = df[df["pmid"].str.startswith("PMID:", na=False)].copy()
    df["pmid"] = df["pmid"].str.replace("PMID:", "", regex=False)
    return df


def fetch_mesh(pmids, batch_size, sleep):
    """返回 {pmid: {"title":..., "mesh": [...], "major": [...]}}"""
    out = {}
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]
        for attempt in range(3):
            try:
                handle = Entrez.efetch(db="pubmed", id=",".join(batch),
                                       rettype="xml", retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                break
            except Exception as e:  # 网络抖动时重试
                print(f"  batch {i}: {e}, retry {attempt + 1}", file=sys.stderr)
                time.sleep(3 * (attempt + 1))
        else:
            print(f"  batch {i}: 放弃", file=sys.stderr)
            continue

        for art in records.get("PubmedArticle", []):
            cit = art["MedlineCitation"]
            pmid = str(cit["PMID"])
            mesh, major = [], []
            for mh in cit.get("MeshHeadingList", []):
                d = mh["DescriptorName"]
                mesh.append(str(d))
                if d.attributes.get("MajorTopicYN") == "Y":
                    major.append(str(d))
            out[pmid] = {"title": str(cit["Article"].get("ArticleTitle", "")),
                         "mesh": mesh, "major": major}
        print(f"  已获取 {min(i + batch_size, len(pmids))}/{len(pmids)}", file=sys.stderr)
        time.sleep(sleep)
    return out


def classify(mesh):
    s = set(mesh)
    if not s:
        return "no_mesh"   # 新文献尚未被 MEDLINE 标引，或非 MEDLINE 期刊
    z, h = "Zebrafish" in s, "Humans" in s
    if z and h:
        return "zebrafish+human"
    if z:
        return "zebrafish"
    if h:
        return "human"
    return "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gaf")
    ap.add_argument("--email", required=True, help="NCBI 要求提供联系邮箱")
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"),
                    help="可选，NCBI API key（速率从 3 次/秒提到 10 次/秒）")
    ap.add_argument("--experimental", action="store_true", help="只保留实验证据码")
    ap.add_argument("--limit", type=int, help="只取前 N 个 PMID（测试用）")
    ap.add_argument("--batch-size", type=int, default=200)
    ap.add_argument("--out", default="pmid")
    args = ap.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
    sleep = 0.12 if args.api_key else 0.4

    df = read_gaf(args.gaf)
    if args.experimental:
        df = df[df["evidence"].isin(EXP_CODES)]
    pmids = sorted(df["pmid"].unique(), key=int)
    if args.limit:
        pmids = pmids[:args.limit]
        df = df[df["pmid"].isin(pmids)]
    print(f"{len(df)} 条带 PMID 的注释，{len(pmids)} 个唯一 PMID", file=sys.stderr)

    cache_path = f"{args.out}_cache.json"
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    todo = [p for p in pmids if p not in cache]
    print(f"缓存中已有 {len(pmids) - len(todo)}，需下载 {len(todo)}", file=sys.stderr)
    if todo:
        cache.update(fetch_mesh(todo, args.batch_size, sleep))
        json.dump(cache, open(cache_path, "w"))

    # 每个 PMID 一行
    ev_by_pmid = df.groupby("pmid")["evidence"].apply(lambda s: ",".join(sorted(set(s))))
    n_by_pmid = df.groupby("pmid").size()
    rows = []
    for p in pmids:
        rec = cache.get(p)
        if rec is None:
            rows.append({"pmid": p, "species_class": "not_found"})
            continue
        mesh = rec["mesh"]
        rows.append({
            "pmid": p,
            "species_class": classify(mesh),
            "organism_mesh": "|".join(m for m in mesh if m in ORGANISM_MESH),
            "n_annotations": n_by_pmid.get(p, 0),
            "evidence_codes": ev_by_pmid.get(p, ""),
            "title": rec["title"],
            "major_mesh": "|".join(rec["major"]),
            "all_mesh": "|".join(mesh),
        })
    pm = pd.DataFrame(rows)
    pm.to_csv(f"{args.out}_pmid_mesh.tsv", sep="\t", index=False)

    ann = df.merge(pm[["pmid", "species_class", "organism_mesh"]], on="pmid", how="left")
    ann[["db_object_id", "symbol", "go_id", "evidence", "aspect", "pmid",
         "species_class", "organism_mesh"]].to_csv(
        f"{args.out}_annotation_mesh.tsv", sep="\t", index=False)

    print("\n== 按 PMID 统计 ==")
    print(pm["species_class"].value_counts().to_string())
    print("\n== 按注释条数统计 ==")
    print(ann["species_class"].value_counts().to_string())
    print(f"\n已写出 {args.out}_pmid_mesh.tsv, {args.out}_annotation_mesh.tsv")


if __name__ == "__main__":
    main()
