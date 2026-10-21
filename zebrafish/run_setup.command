#!/bin/bash
# 双击运行：创建 conda 环境 goa-stats 并跑统计脚本，输出同时写入 run_log.txt
cd "$(dirname "$0")"
exec > >(tee run_log.txt) 2>&1

# 找到 conda
for d in "$HOME/miniconda3" "$HOME/opt/miniconda3" "/opt/miniconda3" "$HOME/anaconda3" "/opt/homebrew/Caskroom/miniconda/base"; do
  if [ -f "$d/etc/profile.d/conda.sh" ]; then source "$d/etc/profile.d/conda.sh"; break; fi
done
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: 找不到 conda"; exit 1
fi
echo "conda: $(command -v conda)  $(conda --version)"

if conda env list | grep -qE '^goa-stats\s'; then
  echo "环境 goa-stats 已存在，更新中..."
  conda env update -n goa-stats -f environment.yml
else
  conda env create -f environment.yml
fi
conda activate goa-stats || { echo "ERROR: 激活失败"; exit 1; }
python -c "import sys, pandas, Bio; print('python', sys.version.split()[0], '| pandas', pandas.__version__, '| biopython', Bio.__version__)"

echo; echo "===== 只看 EXP ====="
python count_exp_go.py goa_zebrafish.gaf --out exp
echo; echo "===== 全部实验证据 ====="
python count_exp_go.py goa_zebrafish.gaf --experimental --out experimental
echo; echo "ALL DONE"
