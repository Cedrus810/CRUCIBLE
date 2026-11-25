# 决策 0027：xtb `--ohess` 停在鞍点时沿虚频自动重启

- 2026-11-25 用户同意。Claude 提议"让 `run_xtb_hessian` 在 xtb 报出显著虚频时从 `xtbhess.xyz`
  自动重启优化，重启次数设上限（比如最多 2 次），每次重启记进 provenance"，用户原话："可以"。
- **背景**：覆盖集诊断（`jobs/coverage_20261124/`）中 fast preset 的 bonded 拟合在 acetone、
  methyl_acetate、styrene、glycine_zwitterion 上被 `check_hessian` 以虚频 < −20 cm⁻¹ 拒绝。
  诊断（`experiments/20261125_saddle_diagnosis/`）：从 xtb 自己写出的 `xtbhess.xyz`（沿虚频扭曲）
  用同一 GFN-FF `--ohess vtight` 重启一次，四个都到达无虚频极小、能量更低（−0.18 ~ −9.4 kJ/mol）、
  通过同一质量门——问题是优化停在鞍点且管线不重启，不是势能面或 Hessian。
- **决定**（`hessian/xtb.py::run_xtb_hessian`）：
  1. 触发条件 = xtb 报告的最低频率低于 `check_hessian` 的同一报错线 `_IMAGINARY_ERROR_CM`
     （−20 cm⁻¹）**且** `xtbhess.xyz` 存在。不用 xtb 自己的 5 cm⁻¹ 阈值：那会改动现在只带警告、
     照常通过的分子。vibspectrum 解析失败 → 不重启（原行为）。
  2. 重启在 `workdir/imag_restart_<k>/` 跑同一 ohess，最多 `max_imaginary_restarts`（默认 2）次；
     返回最后一次的 Hessian。仍有虚频则原样返回，由 `check_hessian` 照常拒绝，不静默放行。
  3. 记录：`CartesianHessian.meta["imaginary_restarts"]`（逐次：起点文件、重启前最低频率、前后能量），
     `derive_bonded` 写入 IR `provenance.bonded.hessian_imaginary_restarts`。**只在发生重启时加键**。
- **不变**：不触发时返回值与改动前逐位相同（实测 ethanol：hessian 文件与覆盖集原运行逐字节相同，
  meta 不变）；第一次 ohess 的命令与文件位置不变；`check_hessian` 阈值与判据不变；
  G0–G3 的分子都能通过质量门，不会触发（未重跑 Gate）。
- **验证**：单测 5 例（不触发 / 小虚频不触发 / 重启成功 / 到上限返回最后一次 / 关闭）+ external
  acetone（从覆盖集 SDF 出发，重启后无虚频）；覆盖集 5 个失败分子重跑 fallback 阶段
  （`experiments/20261125_imag_restart_verify/coverage/`）：acetone、methyl_acetate、styrene、
  glycine_zwitterion 转为通过；oxirane 仍失败（负环角 k，另一原因，决策 0014 同类）。
