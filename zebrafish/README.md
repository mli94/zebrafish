# 统计 EXP 证据的 GO term

## 1. 建环境（只需一次，在 Mac 的终端里）
cd ~/Desktop/rotation2/statistics
conda env create -f environment.yml
conda activate goa-stats

## 2. 运行
python count_exp_go.py goa_zebrafish.gaf                  # 只要 evidence code = EXP
python count_exp_go.py goa_zebrafish.gaf --experimental --out experimental   # EXP IDA IPI IMP IGI IEP

## 输出
exp_go_per_protein.tsv   每个蛋白的 GO term 数（从大到小），含 F/P/C 分类和 GO 列表
exp_proteins_per_go.tsv  每个 GO term 被多少蛋白注释（从大到小）

说明：同一个 (蛋白, GO term) 出现多行只算一次；qualifier 带 NOT 的注释默认去掉（加 --keep-not 保留）。
