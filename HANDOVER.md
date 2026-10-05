# forcefield_compiler 交接状态 — 2026-10-05（Asia/Tokyo）

**版本 0.2.0。G0（MVP-0）与 G1（MVP-1）均已正式通过。** 下一阶段是 MVP-2（MLP 二面角），尚未开始。

文档：[README](README.md)（[中文](README.zh-CN.md)、[日本語](README.ja.md)）、[CHANGELOG](CHANGELOG.md)、
[决策索引](docs/decisions/README.md)、[理论文档](docs/theory/hessian_methods.md)。

## 正式运行

| Gate | 运行目录 | 结果 |
|---|---|---|
| G0 | `jobs/gate_mvp0_20261005_3/` | passed：11 固件 × 4 构型，8 条腿 2596 行 0 失败 |
| G1 | `jobs/gate_mvp1_20261005_4/` | passed：G1.1–G1.7 全部通过；1 条已知方法局限（见下） |

两次均为 AMBER 原生读写之后的运行（决策 0015）；pytest 非 external 355、external 24，
0 失败 0 跳过。历史运行全部保留不覆盖：`gate_mvp0_20261004*`、`gate_mvp0_20261005`（首个
通过）、`gate_mvp0_20261005_2` / `gate_mvp1_20261005_3`（一条测试的子进程漏传 PYTHONPATH
而失败，已修）、`gate_mvp1_20261005`、`gate_mvp1_20261005_2`。

## MVP-1 结论

- Sobtop oracle（12 case，§2.8 容差：k、键长、键角）：mSeminario、m2Seminario、DRIH 全部 12/12。
  - mSeminario：Allen 2018 本征投影 + 相邻键角缩放因子，bond 不取平均（决策 0011）。
  - m2Seminario = 计划 C2 候选，方法名 `m2seminario`，D2 结案（决策 0011）。
  - DRIH：内部使用原子单位（Bohr、E_h），用户已确认（决策 0012）。
- 全局拟合自洽恢复 6.0e-9；拟合目标是参考 Hessian `h.hessian`。
- 已知方法局限：环丙烷 / GFN-FF / m2Seminario 的环角 k ≈ −33.6，Sobtop 在同一 Hessian 上
  结果相同（决策 0014）。`derive_bonded` 对非正 k 抛 `NonPositiveForceConstantError`，不自动回退。
- 研究指标（频率 RMSE 中位数，cm⁻¹）：GFN-FF 下全局拟合最好（57），GFN2 下 mSeminario
  与 GAFF2 相当（≈110）。只报告、不设门槛。

## 验收约定与实现的变化（都有决策记录）

- 0013：mixed 精度 GROMACS 的 gmx 腿用 `1e-4 + 1e-5·|E|`（fp32 舍入；同一份 top/itp 在 fp64 的
  openmm-top 腿 0 失败）；double 版仍用严格容差。乙腈（近线性角）移出 G0/G1，只作诊断。
- 0015：AMBER prmtop/rst7 原生读写，ParmEd 不再是依赖（只在对照测试里可选使用）。
- 0014：方法给出非正 k 且 Sobtop 同 Hessian 复算确认 → 已知方法局限，不计入 G1.5。

## 环境

- 主机 node-A（40 物理核 / 80 线程，与他人共用）；`source scripts/env.sh`（openmm_dev，Python 3.12）。
- GROMACS 2026.3 mixed：`/home/user/gmx26.3/bin/gmx`（`FFC_GMX` 可替换）。
- xtb 6.7.1：`/home/user/software/xtb-dist/bin/xtb`（`FFC_XTB` 可替换，不查 PATH）。
- Sobtop 2026.1.16：复制在 `tools/sobtop/`（gitignore），运行需 `OMP_STACKSIZE=64M OMP_NUM_THREADS=4`。
- Gate G1：`python scripts/gate_mvp1.py --jobs 6`（每个 xtb 4 线程；按物理核留余量）。

## 待办

1. 理论文档审阅（G1.8，人工）；Allen 2018 原文公式编号待补（理论文档 §2）。
2. 在干净提交上重跑 G0/G1：目前通过的运行，manifest 记录的是 `c20d6a2` 加未提交改动。

## 本地 Git

工作区 `/home/user/forcefield_compiler`，元数据在 `/home/user/.local/share/git-sep/ffcompiler.git`。
挂载所有权不同，查看时用 `git -c safe.directory="$PWD" log --oneline`。全部提交仅在本地，
作者为仓库局部身份 `本机代理`（决策 0005）。MVP-1 验收后的提交从 `debcb30` 开始。
