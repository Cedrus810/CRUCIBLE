# 2026-10-10 AM1-BCC 复现性与 NAGL 量级冒烟（MVP-3 G3 判据依据）

node-A 本机，openmm_dev，单线程。

- `repro.py` → `repro.log`、`results.json`、`runs/<分子>/`：按 `scripts/make_fixtures.py` 的做法
  （`antechamber -c bcc -at gaff2 -nc Q -pf y`，再用同一个 `_equalize_mol2_charges` 均摊残差）在固件 SDF 上重跑，
  与固件 mol2 电荷逐原子比较。**12/12 固件均摊后最大差 0.0**（逐位相同）；均摊前（4 位小数）最大差 3.3e-4，
  即均摊量本身。每个分子 0.5–2.4 s。
- `nagl.py` → `nagl.log`、`nagl_results.json`：openff-gnn-am1bcc-1.0.0 与固件 AM1-BCC 的差：
  逐原子最大 4.3e-5（环丙烷）到 0.040 e（乙醇），RMS 2.5e-5 到 0.018 e；总电荷误差 ≤ 1.8e-7（float32 舍入量级），
  所以写回 IR 前要按同一规则均摊残差。每个分子 0.003–0.74 s。12/12 都在模型适用域内（`check_domains=True` 未报错）。
- 跨机器（node-B/183）的 sqm 复现性**未测**。
