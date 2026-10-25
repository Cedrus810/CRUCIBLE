# 决策 0022：吸收 Sobtop 逆向成果（感知层 + 数据文件解析）

## 背景

对 Sobtop 2026.1.16（Linux 版，未 strip）做了行为级逆向：Ghidra 反编译全部
86 个自有符号 + 12 组黑盒动态实验（乙醇/乙烯/苯/SiO₂，输出逐字段对照）。
完整报告入库 `docs/theory/sobtop_reverse_engineering.md`；反编译产物与沙箱
在 `/tmp/sobtop_re/`（会话内临时）。本决策记录吸收进代码的部分与两处对既有
理解的边界修正。

## 吸收内容

1. **`ffcompiler/perception/`（新包）**：几何化学感知，规则与 Sobtop 逐条对齐——
   - Cordero 2008 共价半径表（110 元素，从 .rodata 逐值抠出）；
   - 成键判据 `d < (r_i+r_j)×1.15`，bondcrit.dat 元素对**距离窗口**（注意：
     是 [下限,上限] Å 区间，不是缩放因子）覆盖；
   - 周期最近镜像（三斜，分数坐标 round 法，语义同 Sobtop getbondvec）；
   - 角 = 中心原子邻居两两组合（含氢、不去重、不剔线性）；
   - 二面 = 全部 i-j-k-l 路径，任一侧角 `180°−θ<1°` 剔除，规范 i<l；
   - improper 感知 = **数据驱动**：配位数 3 ∧ 类型 ∈ 库 $IMPROPER 第三列集合，
     邻居升序 (j₁,j₂,center,j₃)。
2. **`ffcompiler/fflib/`（新包）**：Sobtop 发行数据文件的解析与匹配——
   LJ_param.dat / bonded_param.dat / bondcrit.dat 解析；匹配优先级
   **先特异性分层、层内正序先于反序**（实证 oh-c3-c3-hc 命中反序精确条目
   hc-c3-c3-oh 而非正序通配 X-c3-c3-X）；IDIVF 除入势垒；负 pn 按 cos 偶
   函数等价改写（pn→|pn|、phase→−phase；phase=0 时不变，Sobtop 输出即此）；
   LJ 转换 σ=2R·2^(−1/6)·0.1 nm、ε×4.184；缺参"乱猜"常量（键 2.5e5
   kJ/mol/nm²、角 500 kJ/mol/rad²、improper X-X-ca-ha 等价 1.1 kcal/180°/2）。
3. **单测 24 个**：全部以沙箱 itp 输出为 ground truth（乙醇 8/13/12 列表、
   苯 6 条 improper、乙烯 improper 与 UFF 不触发、LJ 换算逐位、通配优先级、
   IDIVF/负 pn、PBC 正交/三斜镜像、bondcrit 窗口）。

## 对既有理解的两处修正（不改代码，改认知）

- **伪逆阈值是复现边界条件**：Sobtop 对 B 直接 SVD、奇异值绝对阈值
  |σ|>1e-10（等价 G 本征值 λ>1e-20）；本项目 DRIH 用 λ<1e-8·λmax 相对
  置零（决策 0012）。两者在 B 谱双峰（O(0.1–1) 与数值零）时切在空处——
  12/12 一致由此而来；**病态几何（近线性、近冗余 RIC）下两种策略可能分叉**，
  已在 theory §4.1 加注。
- **f 因子适用范围**：theory §2/§3 原文已正确地对 mSeminario 与 m2Seminario
  都乘相邻键角缩放因子（逆向确证 Sobtop 同此）；键参数下 Seminario ≡
  mSeminario（同一 Σλ|u·v| 投影），两方法差异只在角/二面的投影方式。

## 明确不吸收的

- 隐藏 k_method=5（局部 Hessian + 双伪逆投影）：Sobtop 菜单无入口的半成品；
- funct 3 键-键交叉项导出：仅 GROMACS 可表达，与三引擎等价性目标冲突
  （RIC Hessian 全矩阵仍在，需要时可作诊断输出复活）；
- "全方法"模式的谐振扭转（锁死构象）、缺失参数乱猜策略：作为**反面教材**
  记录在报告 §8，不进入管线。

## 顺带修正

早期会话把 bondcrit.dat 误记为"因子覆盖"，逆向确认为**距离窗口**
（'Fe Cl 2.0 2.4' = Fe–Cl 距离落在 [2.0,2.4] Å 内才成键）；实现与报告
均已按窗口语义书写。

## 关联

- 报告：`docs/theory/sobtop_reverse_engineering.md`
- 决策 0011（公式核对）、0012（DRIH 单位）——本决策为其补充边界条件
- 后续（未立项）：mol2 直读入口、周期体系 IR/导出、UFF LJ 进默认数据路径
