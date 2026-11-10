#!/usr/bin/env bash
# 2026-10-28：决策 0026 之后的全量 external 测试，只靠用户配置文件（不改 HOME、不设 FFC_<工具>）。
set -eo pipefail
export FFC_MAMBA_ROOT=/home/ruigengji/miniforge3 FFC_MODULEFILES=/home/ruigengji/modulefiles
cd /home/ruigengji/CRUCIBLE
source scripts/env.sh >/dev/null 2>&1
set -u
export PYTHONPATH=src OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
env | grep '^FFC_' | sort
python -m ffcompiler.cli check-env
python -m pytest -m external -p no:cacheprovider -q -rf \
    --junitxml=experiments/20261028_tool_config_verify/junit_external.xml
