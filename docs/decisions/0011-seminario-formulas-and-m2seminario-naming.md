# 决策 0011：mSeminario、m2Seminario 公式经 Sobtop 核对确定；方法名定为 m2seminario（取代 0009）

- 2026-10-05 验收审阅发现：决策 0009 尝试的 8 个 mSeminario 变体都不含经典
  Seminario 的本征投影式 k = Σλᵢ|u·vᵢ|；m2Seminario 的 C2 候选未实现。
- **mSeminario**（`bonded/mseminario.py`，公式见理论文档 §2）：
  - K_AB = −H_AB 不对称化，本征对取实部，k(u) = Σλᵢ|u·vᵢ|；
  - bond 取输入书写方向 (i, j) 的 k_ij，**不与 k_ji 平均**（平均或反向都不一致）；
  - angle 用 Allen 2018 的相邻键角缩放因子 f = 1 + Σ|u_PA·u_PA′|²/(n−1)，
    对两种端点顺序取平均。
- **m2Seminario**：按计划 §2.6 实现 C1、C2 两个候选（`bonded/m2seminario.py`，
  C1 为私有 `_c1`）。C2 = C1 × 同一缩放因子。
- 比对（`scripts/sobtop_oracle.py`，§2.8 容差：k、键长、键角，除 sf6_fchk）：

  | case | C1 bond 超差 | C1 angle 超差 | C2 超差 | mSeminario 超差 |
  |---|---|---|---|---|
  | acetate_gfn2 | 0/6 | 9/9 | 0/15 | 0/15 |
  | benzene_gfn2 | 0/12 | 18/18 | 0/30 | 0/30 |
  | chlorobenzene_gfn2 | 0/12 | 18/18 | 0/30 | 0/30 |
  | cyclobutane_gfn2 | 0/12 | 24/24 | 0/36 | 0/36 |
  | cyclopropane_gfn2 | 0/9 | 18/18 | 0/27 | 0/27 |
  | dmso_gfn2 | 0/9 | 15/15 | 0/24 | 0/24 |
  | ethanol_gfn2 | 0/8 | 12/13 | 0/21 | 0/21 |
  | methyl_benzoate_fchk | 0/18 | 27/28 | 0/46 | 0/46 |
  | methyl_benzoate_gfn2 | 0/18 | 27/28 | 0/46 | 0/46 |
  | methylammonium_gfn2 | 0/7 | 12/12 | 0/19 | 0/19 |
  | nma_gfn2 | 0/11 | 18/18 | 0/29 | 0/29 |
  | norbornane_gfn2 | 0/20 | 42/42 | 0/62 | 0/62 |

  最大相对 k 偏差：mSeminario 2.4e-6，C2 3.8e-6（param.txt 舍入量级）。
- 命名：按 §2.6 规则，C2 全部一致 ⇒ 公开方法名 `m2seminario`。D2 结案。
  模块、CLI、Gate、测试中的 `symseminario` 全部改名。
- oracle 测试（T17 用例 1、T18 用例 1）由 xfail 改为真实断言。
- ~~未完成：Allen et al. 2018 原文公式编号尚未对照填写~~ 2026-10-08 已补：用户批准联网查阅，
  对照 accepted version（`docs/theory/refs/allen2018_jctc_accepted.pdf`）写入理论文档 §2（bond Eq. 4，angle Eq. 5–7、12）。
  bond（Eq. 4）与相邻键角缩放（Eq. 12，求和含自身再减 1，与本项目写法等价）已与原文逐式对照一致；
  原文没写的细节（bond 不取 AB/BA 平均、angle 两种写法取平均、复本征值取实部）以 Sobtop 数值核对为准。
