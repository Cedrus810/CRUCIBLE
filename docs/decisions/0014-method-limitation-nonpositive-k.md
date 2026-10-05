# 决策 0014：方法给出非正力常数且经 Sobtop 同 Hessian 复算确认时，记为已知方法局限

- 2026-10-05 正式 G1（`jobs/gate_mvp1_20261005/`）唯一失败组合：环丙烷 / GFN-FF /
  m2seminario，三个环 C–C–C 角 k ≈ −33.6 kJ/mol/rad²，IR 校验拒绝。其余 109 个组合
  均有结果，G1.1–G1.4、G1.6、G1.7 全部通过。
- 同一 Hessian 交给 Sobtop 2026.1.16：mSeminario 179.4 / 179.4，m2Seminario
  **−33.6 / −33.6**，DRIH 56.5 / 56.5（我们 / Sobtop）。实现与参考一致；负值来自
  m2Seminario 的二次型 uᵀKu 在 GFN-FF 三元环上的投影，是方法在该 Hessian 上的局限，
  不是实现错误。同分子 GFN2 Hessian 上 m2Seminario 为正（oracle 12/12 之一）。
- 用户裁定："可以"（同意不让方法局限卡 Gate）。
- **决定**：
  1. `derive_bonded` 在推导结果中出现非正 k 时抛 `NonPositiveForceConstantError`
     （列出方法、项与 k，提示改用 mseminario 或 drih）。**不自动回退**到库参数。
  2. Gate 遇到该异常时，现场用 Sobtop 在同一 Hessian 上复算（`sobtop_confirm`，
     结果存 `sobtop_confirm/`）。每一个非正项 Sobtop 也为非正且与我们一致（§2.8
     容差）→ 记为 `method_limitation`，写进 `acceptance.json` 的 `method_limitations`
     和 report.md 的"已知方法局限"，不计入 G1.5。Sobtop 失败或任一项不符 → 仍是失败。
  3. 固件 `tests/fixtures/hessian/sobtop_known_limitations/cyclopropane_gfnff/` 保存
     该 Hessian 与 Sobtop 输出；回归测试断言我们与 Sobtop 一致且环角为负。
- 不改变：§2 的方法定义、oracle 矩阵与容差、其他失败类型的判定。
