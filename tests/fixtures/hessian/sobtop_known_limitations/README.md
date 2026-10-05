# Sobtop 已知方法局限用例（决策 0014）

不进 oracle 一致性矩阵（`sobtop_refs/`），只用来证明"方法给出非正 k 时，Sobtop 在
同一 Hessian 上给出同样的非正 k"。

- `cyclopropane_gfnff/`：Gate G1 运行 `jobs/gate_mvp1_20261005/xtb/gfnff/cyclopropane/`
  的 xtb 6.7.1 GFN-FF 输出（`hessian`、`xtbopt.xyz`、`xtb.out`）；`input.mol2` 由
  `make_sobtop_refs.write_ir_mol2` 用 IR 键连 + xtbopt 几何写出；`param_m{2,3,4}.txt`
  为 Sobtop 2026.1.16 方法 2/3/4 的输出（§2.8 stdin 序列）。
  m2Seminario 的三个环 C–C–C 角 k ≈ −33.6 kJ/mol/rad²。
