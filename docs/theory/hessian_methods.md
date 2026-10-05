# Hessian → bond/angle 参数推导方法理论文档（MVP-1 计划 T15）

> **状态：全文写完，待审阅（G1.8）。** §2 的 mSeminario 公式按 Allen et al.
> 2018 的公式结构写出，并以 Sobtop mSeminario 输出逐项核对（12 个 oracle case
> 全部在 §7 容差内，决策 0011）；**原文公式编号尚未对照原文填写**。§4 的
> Wilson B 长度单位（Bohr）是对计划 §2.5 的补全，待用户确认（决策 0012）。

本文档锁定四种方法用到的全部数学定义。计划 §2 的定义如与本文冲突，以发现
矛盾的方式停报，不得自行改。

---

## 1. 记号与总约定

- Cartesian Hessian \(\mathbf H_x\)（\(3N\times3N\)，kJ mol⁻¹ nm⁻²）在参考几何
  \(\mathbf x_0\)（nm）处；梯度 \(\mathbf g\)（kJ mol⁻¹ nm⁻¹）。
- 内坐标：bond \(q_t=|\mathbf r_j-\mathbf r_i|\)（nm）、angle（顶点 \(j\)，rad）、
  dihedral（MVP-0 §2.1 定义，顺式为 0）。
- Wilson B 矩阵 \(\mathbf B\)：行 \(t\) 为 \(\partial q_t/\partial\mathbf x\)
  （\(n_q\times3N\)），即 \(\mathbf b_t\)（\(3N\) 维行向量）。
- 所有方法输出的力常数 \(k\) 都是**能量对内坐标的二阶导数**（½k 约定：
  \(E=\tfrac12k(q-q^0)^2\)），与 ForceFieldIR 一致。

---

## 2. mSeminario（Allen et al. 2018；Sobtop 核对，决策 0011）

记号：\(\mathbf K_{AB}=-\mathbf H_{AB}\)（3×3 原子间块，**不对称化**），
\(\lambda_i^{AB},\mathbf v_i^{AB}\) 为其本征值与本征向量。\(\mathbf K_{AB}\) 一般非对称，
本征对可为复数，取实部。定义沿单位向量 \(\mathbf u\) 的投影力常数

\[
k_{AB}(\mathbf u)=\sum_{i=1}^{3}\lambda_i^{AB}\,\big|\mathbf u\cdot\mathbf v_i^{AB}\big| .
\]

**bond**（原文公式编号：待补）：\(k_{\text{bond}}=k_{AB}(\mathbf u_{AB})\)，
\(\mathbf u_{AB}\) 从 A 指向 B。Sobtop **不对 \(k_{AB}\)、\(k_{BA}\) 取平均**，取键在
mol2 / IR 里的书写方向 (A, B)；取平均或取反方向都与 Sobtop 不一致（决策 0011）。
因此 mSeminario 的键 k 依赖键的书写方向，对称等价键可以得到不同值（例：DMSO
的两条 S–C 键 1087.8 / 1289.8 kJ/mol/Å²），由 §2.7 的对称化消除。

**angle** A–B–C（原文公式编号：待补）：
\(\mathbf u_N=\widehat{\mathbf u_{CB}\times\mathbf u_{AB}}\)，
\(\mathbf u_{PA}=\mathbf u_N\times\mathbf u_{AB}\)，\(\mathbf u_{PC}=\mathbf u_{CB}\times\mathbf u_N\)
（\(\mathbf u_{AB}\) 从 B 指向 A，\(\mathbf u_{CB}\) 从 B 指向 C），

\[
\frac{1}{k_\theta}=\frac{f_A}{R_{AB}^2\,k_{AB}(\mathbf u_{PA})}+\frac{f_C}{R_{CB}^2\,k_{CB}(\mathbf u_{PC})},
\]

再对 (A,B,C) 与 (C,B,A) 两种写法取平均。**相邻键角缩放因子**：以 B 为顶点、
含键 A–B 的其余键角 A–B–D（D ≠ C）共 \(n-1\) 个，

\[
f_A=1+\frac{1}{n-1}\sum_{D}\big|\mathbf u_{PA}(A,B,C)\cdot\mathbf u_{PA}(A,B,D)\big|^2,
\]

没有其余键角时 \(f_A=1\)。\(f\) 依赖传入的 angle 列表（线性角在 library 策略下
先剔除，再算 \(f\)）。

核对：上式在 12 个 oracle case 上 k、键长、键角全部在 §7 容差内（最大相对 k
偏差 2.4e-6，为 param.txt 的舍入）。

---

## 3. m2Seminario：两个候选定义（计划 §2.6）

记号：\(\mathbf K_{AB}=-\mathbf H_{AB}\)（\(3\times3\) 原子间块，\(\mathbf H\)
按原子块切分），\(\mathbf u_{AB}\) 为从 A 指向 B 的单位向量，
\(R_{AB}=|\mathbf r_B-\mathbf r_A|\)。

### 3.1 C1（直接投影）

**bond**：

$$k_{AB}=\mathbf u_{AB}^{\mathsf T}\,\mathbf K_{AB}\,\mathbf u_{AB}.$$

**angle A-B-C**（顶点 B）：

$$\mathbf u_N=\widehat{\mathbf u_{CB}\times\mathbf u_{AB}},\qquad
\mathbf u_{PA}=\mathbf u_N\times\mathbf u_{AB},\qquad
\mathbf u_{PC}=\mathbf u_{CB}\times\mathbf u_N,$$

$$\frac{1}{k_\theta}=\frac{1}{R_{AB}^{2}\,\mathbf u_{PA}^{\mathsf T}\mathbf K_{AB}\mathbf u_{PA}}
+\frac{1}{R_{CB}^{2}\,\mathbf u_{PC}^{\mathsf T}\mathbf K_{CB}\mathbf u_{PC}}.$$

三个投影向量是相互正交的单位向量：\(\mathbf u_N\perp\) 角平面，
\(\mathbf u_{PA}\perp\mathbf u_{AB}\) 且在角平面内，\(\mathbf u_{PC}\perp\mathbf u_{CB}\)
且在角平面内。角平面上键方向的正交补即"垂直于键、带动角"的方向，
与 mSeminario 用本征向量取垂直分量的做法在数学目标上一致，但 C1 不做
本征分解，因此**没有本征向量取绝对值带来的方向歧义**。

**\(k_{AB}=k_{BA}\) 自动成立**：\(\mathbf K_{BA}=\mathbf K_{AB}^{\mathsf T}\)
（\(\mathbf H\) 对称），\(\mathbf u_{BA}=-\mathbf u_{AB}\)，故

$$\mathbf u_{BA}^{\mathsf T}\mathbf K_{BA}\mathbf u_{BA}
=\mathbf u_{AB}^{\mathsf T}\mathbf K_{AB}^{\mathsf T}\mathbf u_{AB}
=\mathbf u_{AB}^{\mathsf T}\mathbf K_{AB}\mathbf u_{AB}.$$

**等价键角相同**：C1 的全部量（\(\mathbf u_N,\mathbf u_{PA},\mathbf u_{PC}\)）
由参考几何唯一确定，\(\mathbf K\) 的二次型在这些方向上的值不依赖任何
人为选取的基。若分子的图自同构（T20 的等价类定义）同时保持几何与
\(\mathbf H\)，则等价类内两角的输入完全相同，输出必然相同。反向担忧
（等价角给出不同值）只出现在 Hessian 本身破坏对称时——那是输入数据
问题，由 T20 的对称化平均处理，不是公式歧义。

### 3.2 C2

bond 同 C1；angle 在 C1 的两项上各乘 §2 的相邻键角缩放因子：

\[
\frac{1}{k_\theta}=\frac{f_A}{R_{AB}^2\,\mathbf u_{PA}^\mathsf T\mathbf K_{AB}\mathbf u_{PA}}
+\frac{f_C}{R_{CB}^2\,\mathbf u_{PC}^\mathsf T\mathbf K_{CB}\mathbf u_{PC}} .
\]

除该乘子外与 C1 完全相同；\(k_{AB}=k_{BA}\) 与等价性论证不变（乘子对
等价角相同）。C2 = mSeminario 把本征投影 \(\sum\lambda|\mathbf u\cdot\mathbf v|\)
换成二次型 \(\mathbf u^\mathsf T\mathbf K\mathbf u\)。

### 3.3 命名规则（决策 D2）

对全部 oracle 用例（T16 `sobtop_refs`，除 `sf6_fchk`）按 §7 容差逐一比对
Sobtop 的 m2Seminario（方法 3）输出：

- 恰有一个候选全部一致 → 该候选实现为 `m2seminario`；
- 两个都不一致 → 实现 C1 为 `symseminario`，把比对表写进决策记录，等用户定。

**结果**（决策 0011）：C2 在 12 个 case 上全部一致，C1 的 bond 全部一致、angle
几乎全部超差（恰差因子 \(f\)，苯上 \(f=2\)）。方法名为 `m2seminario`（C2）。

---

## 4. DRIH（计划 §2.5，锁定）

### 4.1 冗余坐标变换

冗余内坐标集合 \(q\)（\(n_q\) 个）：

> IR 全部 bond ∪ IR 全部非线性 angle ∪ 键图生成的全部 proper dihedral
> 四元组（不含经过线性键角的四元组，不含 improper）。

在参考几何上算 \(\mathbf B\)（\(n_q\times3N\)），\(\mathbf G=\mathbf B\mathbf B^{\mathsf T}\)。
\(n_q>3N-6\) 时 \(\mathbf G\) 奇异，用特征分解求伪逆：特征值
\(\lambda<10^{-8}\lambda_{\max}\) 的方向置零，

$$\mathbf G^{-}=\sum_{\lambda_i\ge10^{-8}\lambda_{\max}}\frac{\mathbf v_i\mathbf v_i^{\mathsf T}}{\lambda_i}.$$

令 \(\mathbf A=\mathbf G^{-}\mathbf B\)（\(n_q\times3N\)），

$$\mathbf H_q=\mathbf A\,\mathbf H_x\,\mathbf A^{\mathsf T}.$$

每个 bond/angle 的力常数取 \((\mathbf H_q)_{tt}\)，平衡值取参考几何上的
\(q_t(\mathbf x_0)\)。只取对角元（Sobtop 的 DRIH 另有 bond–bond 耦合项，
AMBER 型函数表示不了，丢弃）。

### 4.1a 内部单位（决策 0012，用户已确认）

All Cartesian Hessians are internally normalized to atomic units (\(E_h/a_0^2\)) and
Cartesian coordinates to Bohr before the Wilson-\(B\) transformation. The resulting
redundant-internal-coordinate Hessian has coordinate-dependent units: \(E_h/a_0^2\)
for stretches, \(E_h/\mathrm{rad}^2\) for angular coordinates, and mixed units for
coupling terms.

| \(H_q\) 元素 | 坐标 | 自然单位 | 换到 IR |
|---|---|---|---|
| bond 对角元 | \(r\)（Bohr） | \(E_h/a_0^2\) | ×9.375829×10⁵ → kJ mol⁻¹ nm⁻² |
| angle 对角元 | \(\theta\)（rad） | \(E_h/\mathrm{rad}^2\) | ×2625.4996 → kJ mol⁻¹ rad⁻²（**不再除 \(a_0^2\)**） |
| dihedral 对角元 | \(\phi\)（rad） | \(E_h/\mathrm{rad}^2\) | 只进坐标集合，不输出 |
| bond–angle 耦合 | \(r,\theta\) | \(E_h/(a_0\,\mathrm{rad})\) | 丢弃（§2.5 第 5 条） |
| bond–bond 耦合 | \(r_1,r_2\) | \(E_h/a_0^2\) | 丢弃 |

能量单位对 \(\mathbf G^-\) 没有影响（\(\mathbf G=\mathbf B\mathbf B^\mathsf T\) 不含能量，\(H\) 线性进入）；
长度单位有影响（§4.4）。`CartesianHessian` 容器仍按计划 §2.1 用 kJ mol⁻¹ nm⁻²，
换算只发生在 DRIH 的入口和出口。

### 4.2 忽略梯度项的前提

把 Cartesian 位移用内坐标表达 \(\mathbf x=\mathbf x(q)\)，链式法则给出

$$\frac{\partial^2E}{\partial q_i\,\partial q_j}
=\underbrace{\big(\mathbf A\mathbf H_x\mathbf A^{\mathsf T}\big)_{ij}}_{\text{DRIH 保留}}
\;+\;\underbrace{\mathbf g^{\mathsf T}\,\frac{\partial^2\mathbf x}{\partial q_i\,\partial q_j}}_{\text{DRIH 丢弃}} .$$

第二项在 \(\mathbf g=0\)（驻点）时严格为零。因此 DRIH **必须**先通过
§2.3 的驻点检查（`check_hessian`，max|g| ≤ 25 kJ/mol/nm，
`allow_nonstationary` 只用于显式豁免并记录）。

### 4.3 无冗余时的精确恢复与有冗余时的固有混合

若坐标集恰好独立（\(n_q=3N-6\)）且 \(\mathbf H_x=\sum_t k_t\mathbf b_t\mathbf b_t^{\mathsf T}\)，
则 \(\mathbf G^{-}\mathbf G=\mathbf I\)，\(\mathbf A\mathbf b_t=\mathbf e_t\)，

$$\mathbf H_q=\sum_t k_t\,\mathbf e_t\mathbf e_t^{\mathsf T},\qquad (\mathbf H_q)_{tt}=k_t\ \text{精确}.$$

有冗余时 \(\mathbf P=\mathbf G^{-}\mathbf G\) 是非单位投影阵，
\((\mathbf H_q)_{tt}=\sum_u k_u\,\big(\sum_i A_{ti}B_{ui}\big)\) 混入其他坐标的
贡献。典型如 sp³ 中心：6 个键角只有 5 个独立自由度。**这是 DRIH 的固有
性质，不是 bug**（T19 用例 4 记录而不断言）。

### 4.4 冗余集合的选择会影响对角元

\(\mathbf P\) 依赖坐标集合：增删一个 dihedral 会改变 \(\mathbf G\) 的谱与
\(\mathbf A\)，从而改变 \((\mathbf H_q)_{tt}\)。所以 §4.1 的集合定义是锁定的，
实现不得增删坐标类型或改过滤规则。

**长度单位同样影响对角元。** bond 行的 B 元素无量纲，angle/dihedral 行的单位
是 rad/长度；换长度单位相当于对 B 的行做非均匀缩放。冗余关系只出现在同类
坐标之间时（无环分子：角–角、二面角–二面角），\(\mathbf G^{-}\) 的投影与该缩放
无关；冗余关系同时牵涉 bond 与 angle 时（环系），伪逆结果依赖长度单位。
因此单位是定义的一部分：B、G 以 **Bohr** 构造（决策 0012）。

### 4.5 与 Sobtop DRIH 的一致性

T15 审阅时要确认 §4.1 与 Sobtop DRIH（方法 4）结果一致（§7 容差）；
不一致就停下报告，不得改定义凑数。

**结果**（决策 0012）：B 以 Bohr 构造时，12 个 oracle case 全部一致（最大相对
k 偏差 1.3e-5）；以 nm 或 Å 构造时只有 5 个无环分子一致，环系偏差达 1–8 倍。
坐标集合、伪逆阈值、A、对角元均未改动。

---

## 5. 全局 Hessian 拟合（MVP-1b，计划 §2.9，锁定）

### 5.1 模型

待推导的 bond/angle 项为 \(t\)，等价类为 \(c\)（`symmetrize=False` 时每项
自成一类）。每项对 Hessian 的贡献为

$$\mathbf H^{(t)}=k_t\,\mathbf b_t^{\mathsf T}\mathbf b_t
\quad\Longleftrightarrow\quad \mathbf H^{(t)}_{ab}=k_t\,\frac{\partial q_t}{\partial x_a}\frac{\partial q_t}{\partial x_b}.$$

### 5.2 \(k\,\mathbf b^{\mathsf T}\mathbf b\) 形式只在 \(q^0=q(\mathbf x_0)\) 时成立

谐振子项 \(E_t=\tfrac12k_t\,(q_t-q_t^0)^2\) 对 Cartesian 求二阶导：

$$\frac{\partial^2E_t}{\partial x_a\,\partial x_b}
=k_t\,b_{t,a}b_{t,b}+(q_t-q_t^0)\,k_t\,\frac{\partial^2q_t}{\partial x_a\,\partial x_b}.$$

第一项是 \(k\,\mathbf b^{\mathsf T}\mathbf b\)；第二项在且仅在
\(q_t(\mathbf x_0)=q_t^0\) 时为零。这就是实现里**平衡值必须取参考几何上的
\(q_t(\mathbf x_0)\)**（对称化时类内平均）的原因，也是 T22 用例 2 的构造
能精确恢复的前提。若用别处（如优化后几何）的平衡值，模型对 \(k\) 不再
线性，拟合偏置且自洽恢复失效。

### 5.3 设计矩阵与求解

固定部分 \(\mathbf H^{\text{fixed}}\)：IR 去掉全部待推导项后（其余 bond/angle、
全部 torsion、improper、LJ、Coulomb、1-4 原样保留），在 \(\mathbf x_0\) 处的
FF Hessian（OpenMM Reference 解析力中心差分，\(h=10^{-5}\) nm，
\(\tfrac12(\mathbf H+\mathbf H^{\mathsf T})\)）。

以先验为尺度：\(k_c=s_c\,k_c^{\text{prior}}\)，每类一列

$$\mathbf D_{:,c}=k_c^{\text{prior}}\operatorname{vec}\!\Big(\sum_{t\in c}\mathbf b_t^{\mathsf T}\mathbf b_t\Big),
\qquad \mathbf y=\operatorname{vec}\big(\mathbf H^{\text{ref}}-\mathbf H^{\text{fixed}}\big),$$

vec 取全部 \(9N^2\) 个元素、权重 1（决策 D4）。求解

$$\min_{\mathbf s\ge0}\ \|\mathbf y-\mathbf D\mathbf s\|^2+\lambda\|\mathbf s-\mathbf 1\|^2,
\qquad \lambda=\alpha\,\operatorname{tr}(\mathbf D^{\mathsf T}\mathbf D)/n_c,$$

默认 \(\alpha=10^{-3}\)。正则项增广进最小二乘（行块
\(\begin{pmatrix}\mathbf D\\\sqrt\lambda\,\mathbf I\end{pmatrix}\mathbf s\) 对
\(\begin{pmatrix}\mathbf y\\\sqrt\lambda\,\mathbf 1\end{pmatrix}\)），下界 0，
用 `scipy.optimize.lsq_linear`。

**为什么是线性最小二乘、不需要 LM**：未知量只有各等价类的尺度
\(s_c\ge0\)（几何上是平衡值 \(q^0\)，已由几何直接给定，不是拟合变量），
残差 \(\mathbf y-\mathbf D\mathbf s\) 对 \(\mathbf s\) 线性。LM/Gauss–Newton
针对的是残差对参数非线性（例如同时拟合平衡几何、指数）的情形；这里
没有非线性参数，线性最小二乘即全局解。

### 5.4 先验与报告

- 先验 \(k^{\text{prior}}\)：默认 mSeminario（同一 Hessian、同样的对称化）；
  `prior_ir` 给出时直接取其 bond/angle 力常数（须同拓扑）。
- 报告字段：α、λ、\(\operatorname{cond}_2(\mathbf D)\)、拟合后相对残差
  \(\|\mathbf y-\mathbf D\mathbf s\|/\|\mathbf y\|\)、先验（\(\mathbf s=\mathbf 1\)）
  的相对残差、落在下界的类数、每类 \(s_c\)、拟合后 FF 在 \(\mathbf x_0\) 处
  的最大力。

---

## 6. 单位约定与各格式换算

所有方法输出的 \(k\) 都是二阶导（½k 约定）。IR 内部单位：键
kJ mol⁻¹ nm⁻²，角 kJ mol⁻¹ rad²，平衡值 nm / rad。

| 格式 | 键 \(k\) | 角 \(k\) | 换算到 IR |
|---|---|---|---|
| ForceFieldIR | kJ/mol/nm² | kJ/mol/rad² | — |
| Sobtop param.txt | kJ/mol/Å² | kJ/mol/rad²（角值 deg） | 键 ×100；角值 ×π/180 |
| GROMACS itp（funct 1/2、4） | kJ/mol/nm² | kJ/mol/rad² | 不变 |
| AMBER prmtop（\(E=K\Delta q^2\)，无 ½） | kcal/mol/Å² | kcal/mol/rad² | **×2**×4.184；键再 ×100 |

Hessian 输入单位：

| 来源 | 单位 | 换算到 kJ mol⁻¹ nm⁻² |
|---|---|---|
| Gaussian fchk `Cartesian Force Constants` | Ha/Bohr² | ×\(E_h/(a_0)^2\) = HARTREE_KJ_MOL / BOHR_NM² ≈ 937.583…×10³ |
| xtb `hessian` | Ha/Bohr² | 同上 |
| 坐标（fchk） | Bohr | ×BOHR_NM |
| 坐标（xyz/mol2，Å） | Å | ×0.1 |

谐振频率：特征值 λ（kJ mol⁻¹ nm⁻² amu⁻¹）→
\(\tilde\nu=\operatorname{sign}(\lambda)\sqrt{|\lambda|\times10^{24}}/(2\pi c)\)（cm⁻¹）。

交叉核对（SF6，Sobtop mSeminario）：param.txt bond k = 2928.942 kJ/mol/Å²
= 2.928942×10⁵ kJ/mol/nm² = GROMACS itp 值；AMBER 值应为其 ÷836.8
（= ÷(2×4.184×100) = 350.02 kcal/mol/Å²）。此核对在 T16 固件里体现。

**唯一规范约定**：IR 的全部谐振力常数按 \(E=\tfrac12 k\,\Delta q^2\) 储存，即
\(k=\partial^2E/\partial q^2\)（Hessian 对角元）。各后端只在导出时换算：

| 方向 | 换算 |
|---|---|
| DRIH / mSeminario / m2Seminario / 全局拟合 → IR | 不乘 ½（本身就是二阶导） |
| IR → GROMACS、OpenMM | \(k\) 原样 |
| IR → AMBER prmtop | \(K=k/2\)（再做 kcal、Å 换算） |
| AMBER prmtop → IR | \(k=2K\) |

---

## 7. Sobtop 比对容差（计划 §2.8，重抄）

- 力常数：\(|\Delta k|\le10^{-4}|k_{\text{Sobtop}}|+5\times10^{-3}\)（param.txt 单位）；
- 键长：≤ 2×10⁻⁶ Å；键角：≤ 2×10⁻⁴ deg。
- `sf6_fchk` 只作诊断不进 Gate（有 180° 键角，Sobtop 等价角本身不一致）。
