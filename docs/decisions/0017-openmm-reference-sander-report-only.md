# 决策 0017：MLP 二面角 IR 以 OpenMM 为准，G2 的 sander 腿只报告

- 2026-10-10 用户决定。原话："如果 openmm 说力场没问题，那么就是没问题"；
  "这个问题交给后人测试……我没时间管 amber"。
- **背景**：G2 首次正式运行（`jobs/gate_mvp2_20261008_3/`）只有 G2.3 未过，唯一超差的是 nma 的
  sander-ours dihed 行（2.3e-4 > 容差 2e-4 kcal/mol）。拟合相位 3.141059（离 π 5.3e-4 rad）；
  假设 sander 读 prmtop 时把接近 π 的相位改成 π，numpy 重算与 sander 的打印值一致
  （`experiments/20261008_g2_nma_compare/snap_check.log`）。sander 的吸附窗口值**未在源码核实**。
- **决定**：
  1. MLP 二面角拟合出的任意相位原样写进 IR，不吸附成 0/π（计划 §7 第 5 条的修法不采用）。
  2. Gate G2 的 G2.3 只按 sander 以外的腿判定（numpy-ir、openmm-ir/xml/prmtop/top、GROMACS）。
     sander 腿（`sander-ours`、`sander-tleap`）照常运行，超差行写进报告的"已知局限"一节，不计失败。
     AMBER prmtop 的写出仍由 `openmm-prmtop` 腿（OpenMM 读我们的 prmtop）判定。
  3. 已知局限：在 sander / pmemd 里用导出的 prmtop，相位接近 0/π 的拟合项可能被吸附，
     能量偏差上限约 \(\sum_n k_n|\Delta\delta_n|\)。sander 的窗口值与 pmemd 是否同样吸附留给以后的 AMBER 测试。
  4. sander 自身崩溃（`compare` 抛异常）仍记为 G2.2 失败，这不是精度问题。
- **不变**：G0 的 sander 腿照常判定（GAFF2 相位都是精确的 0/π，不涉及吸附）。
