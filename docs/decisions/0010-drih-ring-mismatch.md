# 决策 0010：DRIH 与 Sobtop 在环系/高对称分子上的差异如实记录

> **已被决策 0012 取代**（2026-10-05 验收）：差异来自 Wilson B 的内部长度单位，以 Bohr 构造后 12/12 一致。下文保留为历史记录。

- DRIH 按 MVP-1 计划 §2.5 锁定定义实现（冗余集 = bond ∪ 非线性 angle ∪
  键图全部 proper 四元组；G=BBᵀ 特征分解伪逆，<1e-8·λmax 置零；A=G⁻B，
  Hq=AHAᵀ 取对角元）。
- oracle 比对（Sobtop 方法 4，param_m4.txt，§2.8 容差）：
  - **完全一致**：ethanol、acetate、dmso、methylammonium、nma（全部
    bond/angle 逐项在容差内；键值比值恒为 100.0000，即 Å²→nm² 单位换算）。
  - **不一致**：benzene、chlorobenzene、cyclobutane、cyclopropane、
    methyl_benzoate（苯环部分）、norbornane。偏差呈高度结构化的类别比值
    （苯：C–C 键 0.7103、C–C–C 角 7.985、C–C–H 角 2.724；C–H 键仍精确一致）。
- 已排除的假设：伪逆截断阈值（1e-8…1e-5 无变化）、环中心键二面角剔除、
  坐标集规模（Sobtop 报告的 12/18/24 与本实现一致）、单位换算。
- 判断：Sobtop 对高对称/环系的 DRIH 有未公开的额外处理（疑似对简并
  G 谱的不同正则化或额外坐标），无文档、无原文，无法在不改定义的前提下
  复现。按计划 §2.5"不一致就停下报告（不得改定义去凑数）"处理：
  实现保持锁定定义，差异记录于本决策与测试输出，等 Sobtop 文档/原文。
