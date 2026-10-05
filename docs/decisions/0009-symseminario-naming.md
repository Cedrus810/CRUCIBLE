# 决策 0009：m2Seminario 命名落为 symseminario（D2 待定），C2 因缺原文未实现

> **已被决策 0011 取代**（2026-10-05 验收）：C2 在全部 oracle case 上与 Sobtop 一致，方法名定为 `m2seminario`；mSeminario 的经典本征投影式 12/12 一致。下文保留为历史记录；文中的 `scripts/probe_oracle_variants.py` 已由 `scripts/sobtop_oracle.py` 取代。

- 按 MVP-1 计划 §2.6 的命名规则做 oracle 比对（C1 候选 vs Sobtop 方法 3 输出
  `param_m3.txt`，13 个 case，容差 §2.8：|Δk| ≤ 1e-4|k|+5e-3）。
- **结果：C1 的 bond 公式（k = uᵀKu）在全部 13 个 case 的全部 bond 上与
  Sobtop m2Seminario 完全一致**（超差 0 项）。angle 公式不一致（多数 case
  有一半以上 angle 超差；对称分子如苯的面内垂直投影曲率恰为 0，使 C1 的
  调和式坍缩）。
- C2 候选 = C1 × mSeminario 相邻键角缩放因子；该因子出自 Allen et al. 2018
  原文，原文未提供（决策 0008），无法核对实现，因此 C2 未实现。
- 命名：按计划规则，"哪个候选全部一致就命名 m2seminario；都不一致就把
  C1 实现为 symseminario、报告比对结果、等用户定（D2）"。现况是"C1 键一致、
  角不一致"→ 落到 **`symseminario`**，oracle 测试改为 `xfail(strict=True)`
  原因 "D2 待定"，比对数据保留在 `scripts/probe_oracle_variants.py` 及本记录。
- mSeminario（T17，方法 2）同样未与 Sobtop 一致：经典 Seminario 变体
  （eigenmax 键 / ±(λ1+λ2)/2 角等 8 种组合）全部 0/12 匹配，偏差 0.14–0.28
  （角度）与 >2（键，eigenmax vs 投影）；反向工程 Sobtop 私有实现的尝试
  （隐含缩放因子提取）因对称性除零和数据不足中止。待用户提供 Allen et al.
  2018 原文后重做。方法按决策 0008 照常交付。
