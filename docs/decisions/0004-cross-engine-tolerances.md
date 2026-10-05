# 决策 0004：跨引擎容差的数值地板与近线性角的引擎间系统差

- 日期：2026-10-04
- 背景：MVP-0 计划 §2.4 容差在 Gate 试跑（T9，12 固件 × 4 构型）中暴露三类数值问题。以下偏离均有实测数据支撑，其余约定不变。

## 4.1 非键能量的双精度求和地板

- 实测：`openmm-xml`（NonbondedForce 内核）vs `numpy-ir` 的 `nonbonded_all`，在 nma / methylammonium / dmso 上差 1.1–2.7e-6 kJ/mol（相对 ~1e-8），超过 §2.4 内存内容差 `1e-6 + 1e-9|E|`。原因是 NonbondedForce 用预组合 A/B 系数、numpy 用 σ/r 逐对表达式，求和次序与算术路径不同。
- 决定：内存内能量容差放宽为 `1e-6 + 1e-8|E|`（相对余量放大一个量级）。真实约定错误的量级 ≥1e-3 相对，判别能力不受影响。

## 4.2 近线性角的引擎间系统差（sander vs OpenMM/numpy）

- 实测：acetonitrile（C-C≡N 原始几何 θ=179.986°）的 `angle` 项：sander（读 tleap 原版 prmtop）= 0.1088 kcal/mol，numpy-IR = 0.1134 kcal/mol，差 4.6e-3，超 sander 容差 `2e-4 + 5e-5|E|`（静电）/`2e-4`（非静电）约 20 倍；而 OpenMM 读同一文件与 numpy 一致到 1e-4 kJ。这是 sander 与 OpenMM/numpy 在 θ→180° 上的角度求值路径差异，属引擎固有行为，非转换错误。MVP-1 计划已将线性角列为专项（`LinearAngleError`、`linear_angle_policy`，乙腈在 G1 矩阵中用 `"library"` 策略）。
- 决定：`compare()` 在检测到 IR 中存在 θ > 175°（MVP-1 的 `LINEAR_ANGLE_DEG` 阈值）的键角时，把 sander-ours vs numpy-ir 的 `angle` 与 `dihed` 两类行标记为 SKIPPED（原因写入报告），`sander-ours vs sander-tleap`（同引擎、同几何）保留。其余 11 个固件不受影响。

## 4.3 T9 用例 4（篡改检测）的断言范围

- 实测：bond k×1.001 的篡改在 kJ/mol 的经文件容差（1e-4）下被 `openmm-prmtop` 行捕获，但 sander 的绝对容差 `2e-4` kcal/mol 恰好容纳 1.14e-4 kcal 的变化。
- 决定：T9 用例 4 的断言定为"openmm-prmtop 的 bond 行失败"（计划原文"失败行指向 bond 类别"），不要求 sander-ours 也失败。

## 4.4 sander-tleap 腿的几何一致性

- 实现（T9）说明：sander-tleap 与 sander-ours 必须在同一几何上比较。tleap 自带 inpcrd 是原始几何，对扰动构型会虚报失败；`compare()` 因此用 ParmEd 把当前 positions 写成 rst7 配给 tleap prmtop。这是对计划表格"腿"定义的操作性澄清，不改变容差。
