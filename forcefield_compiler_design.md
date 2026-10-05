# 通用小分子力场参数化编译器：设计草案

> 工作定位：不是简单复刻 Sobtop，而是把其已经成熟的拓扑/参数化流程重新拆成模块化、可测试、可扩展的 **force-field compiler**。  
> 核心新增部分：**GFN-FF / xTB / QM Hessian 后端 + 三种局部投影 + 一种全局拟合（QMDFF-inspired）的 Hessian 参数提取方法 + 已有 MLP 二面角模块 + 新型电荷模型 + GROMACS / AMBER / OpenMM 多后端输出**。
>
> 最直接的理论先例是 **QMDFF**（Grimme 2014，见 §6.0）；Sobtop 更像参数/拓扑工具箱。

---

## 1. 核心目标

我们要解决的不是“再写一个 `.itp` 生成器”，而是建立统一的参数化管线：

```text
Molecule
  │
  ├── chemical perception
  │     ├── bond order
  │     ├── aromaticity
  │     ├── ring perception
  │     └── formal charge
  │
  ├── atom typing / base FF
  │     ├── GAFF
  │     ├── GAFF2
  │     └── future/custom
  │
  ├── bonded parameterization
  │     ├── library parameters
  │     └── Hessian-derived parameters
  │
  ├── torsion parameterization
  │     └── our MLP
  │
  ├── charge assignment
  │     ├── RESP / RESP2
  │     ├── AM1-BCC
  │     ├── reference methods
  │     └── our new charge model
  │
  ├── nonbonded/LJ
  │     └── GAFF/GAFF2 initially
  │
  ▼
ForceFieldIR
  │
  ├── GROMACS writer
  ├── AMBER writer
  └── OpenMM writer
```

原则只有一句：

> **参数生成算法与模拟软件输出格式彻底解耦。**

GROMACS、AMBER、OpenMM 只是 backend；真正的核心产品是中间的 `ForceFieldIR`。

---

# 2. 第一阶段的物理模型边界

第一版不试图重新设计整个经典力场。

采用标准 AMBER/GAFF 风格 Hamiltonian：

\[
E =
E_\mathrm{bond}
+
E_\mathrm{angle}
+
E_\mathrm{torsion}
+
E_\mathrm{improper}
+
E_\mathrm{vdW}
+
E_\mathrm{Coulomb}
\]

其中：

| 项 | 第一版来源 |
|---|---|
| bond | GAFF/GAFF2 library 或 Hessian |
| angle | GAFF/GAFF2 library 或 Hessian |
| proper torsion | **我们的 MLP** |
| improper | library / rule-based；特殊刚性体系可 Hessian |
| partial charge | RESP / AM1-BCC baseline + **新 charge model** |
| LJ | GAFF/GAFF2 |
| exclusions / 1–4 | AMBER-family convention |

这样可以把研究工作集中到真正值得动刀的部分，而不是重新制造已经很成熟的 bond/LJ 基础设施。

---

# 3. 软件总体架构

建议拆成以下几个独立层：

```text
┌──────────────────────────────────────────────┐
│               Input / Molecule               │
│ MOL2 / SDF / PDB / XYZ + optional QM data   │
└──────────────────────┬───────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────┐
│             Molecular Graph Layer            │
│ connectivity / aromaticity / rings / charge │
└──────────────────────┬───────────────────────┘
                       │
          ┌────────────┼────────────┐
          ▼            ▼            ▼
     Atom typing    Bonded       Charge
                    engine       engine
                       │
                       ▼
                 Torsion engine
                       │
                       ▼
┌──────────────────────────────────────────────┐
│               ForceFieldIR                   │
│ engine-independent complete Hamiltonian     │
└──────────────────────┬───────────────────────┘
          ┌────────────┼────────────┐
          ▼            ▼            ▼
       GROMACS        AMBER        OpenMM
```

任何算法都不允许直接写 `.itp`、`prmtop` 或 OpenMM XML。

---

# 4. ForceFieldIR：整个项目最重要的一层

推荐内部数据结构大致如下：

```python
@dataclass
class Atom:
    index: int
    element: str
    atom_type: str
    formal_charge: int
    partial_charge: float
    sigma: float
    epsilon: float

@dataclass
class Bond:
    i: int
    j: int
    r0: float
    k: float
    source: str
    uncertainty: float | None = None

@dataclass
class Angle:
    i: int
    j: int
    k_atom: int
    theta0: float
    force_constant: float
    source: str
    uncertainty: float | None = None

@dataclass
class ProperTorsionTerm:
    i: int
    j: int
    k: int
    l: int
    periodicity: int
    phase: float
    force_constant: float
    source: str

@dataclass
class ImproperTorsion:
    ...

@dataclass
class ForceFieldIR:
    atoms: list[Atom]
    bonds: list[Bond]
    angles: list[Angle]
    proper_torsions: list[ProperTorsionTerm]
    impropers: list[ImproperTorsion]
    exclusions: list
    pairs14: list

    forcefield_family: str
    combination_rule: str
    coulomb_14_scale: float
    vdw_14_scale: float

    provenance: dict
```

### 4.1 每个参数必须带 provenance

例如：

```yaml
bond_12:
  source: m2seminario
  hessian_provider: gfnff
  reference_geometry: optimized
  value: ...

torsion_4:
  source: mlp
  model_version: 0.4.2
  uncertainty: 0.07

atom_7_charge:
  source: new_charge_model
  model_version: 0.1.0
```

这对：

- 调试
- benchmark
- 重现论文
- 参数回溯
- 自动 fallback

都会非常重要。

---

# 5. Hessian 系统

## 5.1 Hessian Provider 和 Hessian → FF 参数算法必须分开

这是整个设计中非常重要的一点。

### Hessian Provider

```text
HessianProvider
├── Gaussian
├── ORCA
├── xTB / GFN1-xTB
├── xTB / GFN2-xTB
├── GFN-FF
├── CP2K
└── user-provided Cartesian Hessian
```

它们的工作只有一个：

> 输出统一单位、统一原子顺序的 Cartesian Hessian。

即：

\[
\mathbf H_{ij}
=
\frac{\partial^2 E}
{\partial x_i \partial x_j}
\]

后续参数提取算法完全不知道 Hessian 来自 Gaussian、ORCA、GFN-FF 还是 CP2K。

---

# 6. Hessian → bonded 参数：三种局部投影 + 一种全局拟合

## 6.0 先例：QMDFF

QMDFF（Grimme, JCTC 2014）已经证明：

> **一个分子的 equilibrium geometry + Hessian + partial charges + bond orders，可以自动生成一个 molecule-specific classical force field。**

它的分工：

```text
bond/angle            ← 拟合 QM Hessian
rotatable torsion     ← 自动构造局部模型体系 → 完整 2π 转动能量曲线
                        （最小价基紧束缚 = extended Hückel 哈密顿量）→ 拟合 Fourier 项
nonrotatable torsion  ← Hessian 拟合
TB 失败的 torsion     ← 退回 Hessian 拟合
charge                ← CM5
dispersion            ← D3
repulsion             ← 经验物理模型
```

原文：

> In QMDFF, torsions about rotatable bonds are generally parameterized from full \(2\pi\) rotational energy profiles of automatically constructed local model systems, evaluated with a minimal-valence-basis tight-binding (extended Hückel) Hamiltonian and fitted to Fourier terms; nonrotatable torsions, as well as failed TB cases, are treated through Hessian fitting.

与我们自己走出来的结论一致：**Hessian 对 bond/angle 和刚性 torsion 很好，可旋转 torsion 要另外处理。** 分工也几乎一一对应：

| | QMDFF | 本项目 |
|---|---|---|
| 可旋转 torsion | EHT 模型体系 2π 扫描 → Fourier | MLP（§10） |
| 刚性 torsion | Hessian 拟合 | m2Seminario / DRIH / global fit（§10.2） |
| torsion 失败兜底 | 退回 Hessian 拟合 | 需要定义（§16 fallback graph） |

最后一行提醒：MLP 推理失败或不确定度过高时，我们也需要一条明确的 fallback，而不是静默留空。

但 QMDFF 并不是 Seminario：它不逐项投影，而是让**整个 FF 的 Hessian 去拟合 QM Hessian**，同时优化全部力常数。这给我们提供了第四种 Hessian 方法（§6.4）。

QMDFF 也提醒了它的代价：它是 **non-transferable** 的（molecule A → FF_A，molecule B → FF_B），而且它自己指出参数化的瓶颈基本是 QM Hessian，而不是后面的 FF 生成。我们的定位是站在两者之间：

```text
GAFF2
 │
 ├─ LJ                         keep
 ├─ topology conventions       keep
 ├─ exclusions / 1-4           keep
 │
 ├─ bond/angle                 Hessian refine（molecule-specific）
 ├─ torsion                    MLP replace
 └─ charge                     new model
```

再加上 GFN-FF / GFN2-xTB 把 QMDFF 的瓶颈（Hessian）便宜了一到几个数量级（§8）。因此可以把本项目表述为：

> **modern modular reinterpretation of the QMDFF philosophy**

而不是 QMDFF clone。论文里第四种方法一律称 **QMDFF-inspired global Hessian fit**：函数形式、非键模型都和 QMDFF 不同，不能写成"实现了 QMDFF"。

## 6.0.1 方法清单

Sobtop 文档目前列出了四种 Hessian 派生方案：原始 Seminario、modified Seminario、m2Seminario 和 DRIH。

但我们不需要把原始 Seminario 作为正式推荐路线。原始方法已经基本被 modified 方法取代。

因此我们的核心实现定为四种，分两类：

```text
局部投影法（每个力常数只看局部 Hessian 块）
1. mSeminario
2. m2Seminario
3. DRIH

全局反演法（所有力常数一起拟合整张 Hessian）
4. Global Hessian Fit   ← QMDFF-inspired
```

原始 Seminario 可以：

```text
legacy / compatibility only
```

甚至第一版完全不暴露给普通用户。

---

## 6.1 mSeminario

定位：

> 默认通用 Hessian bond/angle 参数生成方法。

适合：

- 普通有机小分子
- 可旋转键较多的 drug-like molecule
- 与独立 torsion 模块混合使用

推荐默认：

```yaml
bonded:
  method: mseminario
  derive:
    - bond
    - angle
```

不默认从 Hessian 生成 flexible torsion。

---

## 6.2 m2Seminario

定位：

> 更强调对称性和方向一致性的 Seminario 改进版本。

重点解决：

- 对称等价 bond/angle 得到不同力常数
- `i-j` 与 `j-i` 顺序产生微小差异

适合：

- 高对称体系
- 刚性配体
- 无机/配位片段
- 对参数对称性要求较高的体系

推荐：

```yaml
bonded:
  method: m2seminario
```

---

## 6.3 DRIH

DRIH：

> Diagonal elements of Redundant Internal-coordinate Hessian

基本路线：

```text
Cartesian Hessian
      ↓
Wilson B matrix
      ↓
redundant internal-coordinate Hessian
      ↓
take relevant diagonal elements
      ↓
bond / angle / rigid torsion force constants
```

优点：

- 理论映射更直接
- 对高对称刚性体系有优势
- 可以保留一定的内部坐标耦合信息

缺点：

- 对某些分子的二面角可能过软
- 与 GAFF-style flexible torsion 混用时未必最佳

因此 DRIH 不做 drug-like molecule 的默认方案。

主要用途：

```text
rigid molecule
metal-containing fragment
high-symmetry fragment
special bonded topology
```

---

## 6.4 Global Hessian Fit（QMDFF-inspired）

定位：

> 全局反演法。研究 / benchmark 用，第一版不做默认。

与前三种的根本区别：前三种逐项局部投影；本法把所有待定 bonded 力常数 \(\mathbf k\) 当未知数，让 FF 的 Cartesian Hessian 整体逼近参考 Hessian。

### 6.4.1 拟合目标：必须先扣掉固定部分

我们的 FF 里不只有 bond/angle。在参考几何 \(\mathbf x_0\) 上：

\[
\mathbf H^\mathrm{FF}(\mathbf k)
=
\mathbf H^\mathrm{fixed}
+
\sum_t k_t\,\mathbf G_t
\]

其中 \(\mathbf H^\mathrm{fixed}\) 是**不参与拟合**的项在 \(\mathbf x_0\) 处的 Hessian：

```text
H_fixed = H_LJ (GAFF2) + H_Coulomb (当前 charge provider)
        + H_1-4 (按 IR 的 scale) + H_torsion (MLP) + H_improper (若不拟合)
```

对谐振项 \(E_t = \tfrac12 k_t (q_t - q_t^0)^2\)，取 \(q_t^0 = q_t(\mathbf x_0)\)，则

\[
\mathbf G_t = \mathbf b_t \mathbf b_t^\mathsf T,
\qquad
\mathbf b_t = \partial q_t / \partial \mathbf x \big|_{\mathbf x_0}
\]

（\(\tfrac12\) 的有无跟随 IR 锁定的约定，见 §4。）

因此在 \(q^0\) 固定时，问题对 \(\mathbf k\) 是**线性最小二乘**，不需要 Levenberg–Marquardt（QMDFF 用 LM 是因为它的函数形式非线性）：

\[
\min_{\mathbf k}
\left\|
\mathbf W \odot
\left(
\mathbf H^\mathrm{ref} - \mathbf H^\mathrm{fixed} - \sum_t k_t \mathbf G_t
\right)
\right\|_F^2
+
\lambda \left\| \mathbf k - \mathbf k^\mathrm{prior} \right\|^2
\]

直接拿 \(\mathbf H^\mathrm{ref}\) 拟合而不扣 \(\mathbf H^\mathrm{fixed}\)，会把非键和二面角的曲率重复算进 bond/angle。这一点其实是全局法相对 Seminario 的优势：Seminario 完全不管非键项对局部曲率的贡献。

### 6.4.2 带来的管线依赖

- 全局拟合依赖 charge 和 torsion 的结果，因此**必须排在 charge provider 和 torsion provider 之后**。前三种局部方法没有这个依赖。
- 换电荷模型或 MLP 版本，bonded 参数就要重新拟合；provenance 必须记录 \(\mathbf H^\mathrm{fixed}\) 来自哪个 charge/torsion 版本。
- \(\mathbf H^\mathrm{fixed}\) 由 IR 求值得到（解析，或对力做中心差分；小分子 6N 次力计算可忽略）。

### 6.4.3 必须锁定的细节

1. **平衡位置不一致**：\(\mathbf H^\mathrm{fixed}\) 对应的非键 + 二面角在 \(\mathbf x_0\) 一般有非零梯度，取 \(q^0 = q(\mathbf x_0)\) 时 FF 的极小点不在 \(\mathbf x_0\)。（四种方法都有这个问题，但只有全局法有条件处理。）
   - 第一版：接受，并报告 FF 优化后几何相对 \(\mathbf x_0\) 的偏差；
   - 留接口：联合拟合 \(q^0\)，加梯度匹配项 \(\mathbf g^\mathrm{FF}(\mathbf x_0) = \mathbf 0\)（变成非线性问题）。
2. **病态**：sp3 中心的角冗余（6 个角只有 5 个独立）、平移/转动零空间。必须加正则：\(\mathbf k^\mathrm{prior}\) 取 mSeminario 结果（或 GAFF2 库值），并报告设计矩阵的条件数。
3. **对称性**：对称等价项共享同一个 \(k\)（复用 `bonded/symmetry.py` 的等价类），否则拟合会打破对称性。
4. **权重 \(\mathbf W\)**：Cartesian 原样，还是质量加权（与频率更直接对应）—— 待定；先选一种做默认，另一种做 benchmark。
5. **不可表示的耦合**：\(\mathbf H^\mathrm{ref}\) 里的 stretch–bend 等交叉项，AMBER 型函数表示不了，残差会被分摊到对角 \(k\) 上。这正是要和局部投影比较的地方；必须报告拟合残差（Frobenius 范数）和频率误差。
6. **拟合变量范围**：flexible proper torsion 不进拟合（仍归 MLP）；rigid torsion / improper 可选进入拟合变量。这与 QMDFF 的分工相同：可旋转 torsion 另行处理，不可旋转 torsion 走 Hessian 拟合。

推荐：

```yaml
bonded:
  method: global_fit
  prior: mseminario
  derive:
    - bond
    - angle
```

---

# 7. Hessian 方法选择策略

默认策略可以非常简单：

```text
普通 drug-like molecule
        │
        ├── bond  → mSeminario
        ├── angle → mSeminario
        └── torsion → our MLP

高对称 / 刚性体系
        │
        ├── m2Seminario
        └── optionally DRIH benchmark

金属 / 特殊成键
        │
        ├── m2Seminario
        ├── DRIH
        └── higher-level Hessian if needed

研究 / benchmark（opt-in）
        │
        └── Global Hessian Fit（prior = mSeminario），与上面三种对比
```

未来也可以自动选择：

```python
if high_symmetry or many_equivalent_terms:
    method = "m2seminario"
elif rigid_special_fragment:
    method = "drih"
else:
    method = "mseminario"
```

---

# 8. GFN-FF：廉价 Hessian Provider

GFN-FF 的定位不是：

> “我们的最终 force field。”

而是：

> **cheap local-curvature oracle / Hessian provider**

管线：

```text
molecule
   ↓
GFN-FF optimization
   ↓
GFN-FF Cartesian Hessian
   ↓
mSeminario / m2Seminario / DRIH
   ↓
bond + angle parameters
```

GFN-FF 在 xTB 中支持 geometry optimization 和 frequency/Hessian 计算，因此可以直接成为 Hessian backend。

---

## 8.1 为什么 GFN-FF 特别适合这一层

它的优势不是“比 QM 更正确”，而是：

```text
very cheap
fully automated
broad element coverage
large-system friendly
already provides gradients/Hessian workflow
```

所以可以把 Hessian provider 做成三级：

```text
FAST
GFN-FF
   ↓
MEDIUM
GFN2-xTB
   ↓
HIGH
DFT / ORCA / Gaussian / CP2K
```

例如：

```yaml
hessian:
  provider: gfnff
  optimization: true
```

或者：

```yaml
hessian:
  provider: gfn2
  optimization: true
```

或者：

```yaml
hessian:
  provider: orca
  method: r2scan-3c
```

---

## 8.2 一个必须写进文档的物理解释

GFN-FF Hessian 是：

\[
H^\mathrm{GFNFF}
=
H_\mathrm{bond}
+
H_\mathrm{angle}
+
H_\mathrm{torsion}
+
H_\mathrm{NCI}
+\cdots
\]

因此经 Seminario/DRIH 投影得到的是：

> **effective local harmonic curvature**

而不是直接抄出 GFN-FF 内部原始 bond spring constant。

这与 QM Hessian → Seminario 的思想一致：

> 用完整势能面的局部二阶响应，映射回简化经典 Hamiltonian。

所以这一用法在逻辑上完全成立，但论文里必须明确称为：

```text
GFN-FF-derived effective bonded parameters
```

而不能包装成 QM-derived parameters。

---

# 9. GFN-FF direct topology：保留为实验接口，不做主路径

GFN-FF 自己会生成 system-specific topology 和 derived force constants。

因此理论上可以：

```text
GFN-FF topology
      ↓
extract native parameters
      ↓
convert to ForceFieldIR
```

但这不是主设计。

主设计仍然是：

```text
GFN-FF
  ↓
Cartesian Hessian
  ↓
our common Hessian pipeline
```

理由：

1. 与 ORCA/Gaussian/xTB/CP2K 共用相同算法；
2. 避免把整个框架绑定到 GFN-FF 内部实现；
3. provenance 更清楚；
4. benchmark 更容易；
5. 后续可以自由更换 Hessian level。

---

# 10. 二面角：直接接入现有 MLP

这部分不重新造轮子。

我们的 MLP 成为独立：

```python
class TorsionProvider:
    def parameterize(self, molecule, graph, context):
        ...

class MLPTorsionProvider(TorsionProvider):
    ...
```

整体逻辑：

```text
proper torsion detected
       ↓
MLP inference
       ↓
torsion representation
       ↓
ForceFieldIR
```

---

## 10.1 默认规则

对于普通有机小分子：

```text
bond   ← Hessian / library
angle  ← Hessian / library
torsion ← MLP
```

即：

> 不再让 Hessian 方法负责 flexible proper torsion。

这是和 Sobtop 原始工作流最大的结构性区别之一。

---

## 10.2 刚性二面角例外

例如：

- 金属配位几何
- 共轭刚性区域
- 特殊无机结构
- 非标准环系

可以允许：

```text
rigid torsion
   ← m2Seminario / DRIH
```

因此：

```yaml
torsion:
  flexible:
    provider: mlp
  rigid:
    provider: hessian
```

---

## 10.3 MLP 输出形式

支持两条路径。

### A. Classical-compatible mode

```text
MLP torsion PES
      ↓
fit/distill
      ↓
Σ k_n [1 + cos(nφ - δ_n)]
```

优点：

- GROMACS
- AMBER
- OpenMM

全部原生兼容。

### B. Native OpenMM mode

如果 MLP 本身表达的函数明显优于 Fourier 截断：

```text
MLP representation
      ↓
OpenMM CustomTorsionForce / custom force layer
```

同时可以再生成一个经典 Fourier approximation 供 GROMACS/AMBER 使用。

---

# 11. Partial Charge：第二个真正值得创新的模块

Charge 必须做成 provider，而不是写死：

```text
ChargeProvider
├── ExistingCharges
├── RESP
├── RESP2
├── AM1-BCC
├── optional reference models
└── NewChargeModel
```

---

# 12. 为什么不能只做一个“更准的 ESP 拟合”

传统 fixed-charge 小分子参数中真正长期主流的路线非常集中：

```text
ESP fitting family
└── RESP / RESP2 / variants

semiempirical + correction family
└── AM1-BCC / ABCG2-like development
```

已有工作已经说明：

> 对 hydration free energy 优化得更好，并不自动意味着 protein–ligand binding free energy 也会更好。

因此新 charge model 不应该只优化单一目标。

---

# 13. 新 charge model：推荐研究方向

最有吸引力的方向不是“ML 直接预测 RESP”。

后者主要是：

> surrogate / acceleration

更有方法学价值的是：

> **environment-aware, charge-conserving learned correction model**

基本形式：

\[
q_i
=
q_i^{(0)}
+
\sum_j \Delta q_{ji}
\]

对于每条 bond：

\[
\Delta q_{ij} = -\Delta q_{ji}
\]

因此总电荷天然守恒：

\[
\sum_i q_i = Q_\mathrm{formal}
\]

模型：

```text
base charge
   │
   ├── EEQ
   ├── semiempirical
   ├── cheap QM partition
   └── another simple prior
          │
          ▼
local chemical environment
          │
          ▼
small graph/equivariant model
          │
          ▼
bond charge-transfer corrections Δqij
          │
          ▼
final fixed charges
```

可以理解成：

> **generalized / learned BCC**

但 correction 不再只由离散 bond type 决定，而是连续感知完整局部化学环境。

---

# 14. Charge model 的训练目标

不要只优化一个 observable。

建议多目标：

\[
L =
w_\mathrm{ESP}L_\mathrm{ESP}
+
w_\mu L_\mathrm{dipole}
+
w_\mathrm{HFE}L_\mathrm{hydration}
+
w_\mathrm{solv}L_\mathrm{transfer}
+
w_\mathrm{int}L_\mathrm{interaction}
+
w_\mathrm{reg}L_\mathrm{regularization}
\]

训练数据可以分成：

```text
Gas phase
├── ESP
├── dipole
└── multipole information

Water
├── hydration free energy
└── water–solute interaction energies

Protein-like fragment environments
├── amide
├── carboxylate
├── ammonium
├── guanidinium
├── aromatic fragments
├── heterocycles
└── sulfur-containing groups

Organic solvent
└── transfer free energies
```

重点不是所有目标都必须第一版实现，而是架构上不能把方法锁死成：

```text
HFE-only charge fitting
```

---

# 15. Charge 与 LJ 的边界

必须承认：

> charge 与 Lennard-Jones 参数并不是严格独立的。

GAFF/AMBER 中存在历史上的误差补偿。

因此开发顺序应该是：

### Phase 1

```text
new charge + frozen GAFF2 LJ
```

先测：

- ESP/dipole
- HFE
- solvent transfer
- interaction energies
- RBFE

### Phase 2

如果发现系统性残差与 atom type / LJ 强相关：

```text
joint charge/LJ refinement
```

但不要第一天就把 charge、LJ、torsion 全部一起训练，否则很难知道性能改善来自哪里。

---

# 16. Parameter Assignment 的优先级

建议把所有参数来源设计成明确的 fallback graph。

例如 bond/angle：

```text
1. exact trusted library parameter
        ↓ if missing / low confidence
2. GFN-FF Hessian
        ↓ if suspicious
3. GFN2-xTB Hessian
        ↓ if suspicious / special chemistry
4. DFT Hessian
```

也可以反过来提供：

```yaml
bonded:
  policy: derive_all
```

直接全部从指定 Hessian 得到。

---

# 17. 推荐默认 preset

## 17.1 Fast

用于：

- 大批量筛选
- 自动 workflow
- 参数预生成

```yaml
preset: fast

typing: gaff2

bonded:
  provider: hessian
  hessian: gfnff
  method: mseminario

torsion:
  provider: mlp

charge:
  provider: new_charge

vdw:
  provider: gaff2
```

---

## 17.2 Balanced

```yaml
preset: balanced

typing: gaff2

geometry:
  provider: gfn2

bonded:
  provider: hessian
  hessian: gfn2
  method: mseminario

torsion:
  provider: mlp

charge:
  provider: new_charge

vdw:
  provider: gaff2
```

---

## 17.3 High-accuracy

```yaml
preset: high

geometry:
  provider: orca

bonded:
  provider: hessian
  hessian: orca
  method: mseminario

torsion:
  provider: mlp

charge:
  provider: new_charge

vdw:
  provider: gaff2
```

---

## 17.4 Rigid / high-symmetry

```yaml
preset: rigid

bonded:
  hessian: gfn2
  method: m2seminario

rigid_torsion:
  method: drih

charge:
  provider: new_charge
```

---

# 18. 输出层

## 18.1 GROMACS

输出：

```text
molecule.itp
topol.top
coordinates.gro
```

必要时：

```text
ffbonded.itp
ffnonbonded.itp
```

---

## 18.2 AMBER

输出：

```text
mol2
frcmod
```

并可进一步生成：

```text
prmtop
inpcrd / rst7
```

---

## 18.3 OpenMM

两种模式。

### Portable classical mode

```text
ForceFieldIR
   ↓
standard OpenMM forces
```

例如：

```text
HarmonicBondForce
HarmonicAngleForce
PeriodicTorsionForce
NonbondedForce
```

### Native enhanced mode

```text
MLP/custom torsion
        ↓
CustomTorsionForce / specialized force
```

---

# 19. Backend equivalence test

这是整个项目非常值得做的一套自动测试。

同一个 `ForceFieldIR`：

```text
             ForceFieldIR
          /       |        \
         /        |         \
    GROMACS     AMBER      OpenMM
         \        |         /
          \       |        /
             same geometry
                  ↓
         energy / force compare
```

比较：

```text
E_bond
E_angle
E_torsion
E_LJ
E_Coulomb
E_1-4
total energy
atomic forces
```

要求：

\[
\Delta E \rightarrow \text{numerical tolerance}
\]

以及：

\[
\mathrm{RMSD}(\Delta F)
\rightarrow \text{numerical tolerance}
\]

这可以成为项目本身很强的工程卖点：

> **cross-engine Hamiltonian equivalence validation**

---

# 20. 参数本身的验证层级

不能只验证“程序没崩”。

## Level 1：局部参数

```text
bond equilibrium
angle equilibrium
force constants
torsion profile
charge moments
```

## Level 2：分子局部性质

```text
optimized geometry
vibrational frequencies
bond distributions
angle distributions
torsion populations
```

## Level 3：凝聚相

```text
density
heat of vaporization
hydration free energy
solvent transfer free energy
```

## Level 4：药物相关

```text
RBFE
ABFE
pose stability
ranking
```

这样每次替换一个模块，都能知道它真正改善了哪一层。

---

# 21. Repository 结构建议

```text
ffcompiler/
├── pyproject.toml
├── README.md
├── src/
│   └── ffcompiler/
│       ├── molecule/
│       │   ├── graph.py
│       │   ├── perception.py
│       │   └── io.py
│       │
│       ├── typing/
│       │   ├── base.py
│       │   ├── gaff.py
│       │   └── gaff2.py
│       │
│       ├── hessian/
│       │   ├── base.py
│       │   ├── gaussian.py
│       │   ├── orca.py
│       │   ├── xtb.py
│       │   ├── gfnff.py
│       │   └── cp2k.py
│       │
│       ├── bonded/
│       │   ├── mseminario.py
│       │   ├── m2seminario.py
│       │   ├── drih.py
│       │   ├── global_fit.py
│       │   └── symmetry.py
│       │
│       ├── torsion/
│       │   ├── base.py
│       │   ├── mlp.py
│       │   └── fourier.py
│       │
│       ├── charge/
│       │   ├── base.py
│       │   ├── resp.py
│       │   ├── am1bcc.py
│       │   └── learned_bcc.py
│       │
│       ├── ir/
│       │   ├── schema.py
│       │   └── validation.py
│       │
│       ├── exporters/
│       │   ├── gromacs.py
│       │   ├── amber.py
│       │   └── openmm.py
│       │
│       └── cli.py
│
├── tests/
│   ├── unit/
│   ├── regression/
│   ├── hessian/
│   ├── torsion/
│   ├── charge/
│   └── cross_engine/
│
└── benchmarks/
    ├── bonded/
    ├── torsion/
    ├── freesolv/
    └── rbfe/
```

---

# 22. CLI 草案

最简单：

```bash
ffcompiler ligand.sdf \
    --preset fast \
    --output gromacs
```

显式：

```bash
ffcompiler ligand.sdf \
    --typing gaff2 \
    --hessian gfnff \
    --bonded-method mseminario \
    --torsion mlp \
    --charge learned-bcc \
    --output amber
```

高精度：

```bash
ffcompiler ligand.sdf \
    --hessian orca \
    --bonded-method mseminario \
    --torsion mlp \
    --charge learned-bcc \
    --output openmm
```

多后端一次输出：

```bash
ffcompiler ligand.sdf \
    --preset balanced \
    --output gromacs amber openmm
```

---

# 23. MVP 顺序

不要同时把所有东西写完。

## MVP-0：IR + exporters

```text
existing fully parameterized molecule
          ↓
ForceFieldIR
          ↓
GMX / AMBER / OpenMM
```

先证明三后端 Hamiltonian 一致。

---

## MVP-1：Hessian subsystem

实现：

```text
ORCA / xTB / GFN-FF Hessian parser
                ↓
mSeminario
m2Seminario
DRIH
```

先只做 bond + angle。

### MVP-1b：Global Hessian Fit

在三种局部方法之后做。此时 MLP 和新电荷还没接入，\(\mathbf H^\mathrm{fixed}\) 先用 GAFF2 torsion + AM1-BCC（或载入的电荷）+ GAFF2 LJ；接口按 §6.4.2 设计，后续换 provider 时只需重拟合。

---

## MVP-2：MLP torsion 接入

```text
molecule
  ↓
our existing MLP
  ↓
ForceFieldIR torsion terms
```

然后做 torsion-profile regression。

---

## MVP-3：charge provider framework

先支持：

```text
load charges
RESP
AM1-BCC
```

把接口稳定下来。

---

## MVP-4：new charge model

实现：

```text
base charge
+
environment-aware conservative bond correction
```

然后独立 benchmark。

---

## MVP-5：自动 fallback / uncertainty

最终才加入：

```text
GFN-FF
  ↓ uncertain
GFN2-xTB
  ↓ uncertain
DFT
```

不要一开始就搞复杂调度。

---

# 24. 最终科学故事

这个项目可以有三层贡献，而不是把所有东西塞成一个无法解释的“大模型”。

## Contribution A：基础设施

**Engine-independent force-field compiler**

```text
one parameterization
     ↓
GROMACS / AMBER / OpenMM
```

带 cross-engine energy/force equivalence test。

---

## Contribution B：bonded 参数化

**Multi-level Hessian parameterization**

```text
GFN-FF → GFN-xTB → DFT
```

统一经过：

```text
mSeminario / m2Seminario / DRIH        （局部投影）
Global Hessian Fit                     （全局反演，QMDFF-inspired）
```

把 GFN-FF 变成一个真正实用的 cheap Hessian backend。

由此可以直接回答两个 benchmark 问题：

1. **local Hessian projection vs global Hessian reconstruction**：哪个对现代小分子 FF 更好？
2. **QMDFF 式 Hessian 参数化到底需要多高等级的 Hessian？**

```text
reference:  r2SCAN-3c / DFT Hessian
compare:    GFN2-xTB Hessian
            GFN-FF Hessian
```

在 §20 各层级上比较最终 FF：geometry、vibrational frequencies、bond/angle distributions、conformer energies、凝聚相性质。如果 GFN-FF 已经能恢复大部分局部 curvature，这本身就是一个结果。

---

## Contribution C：真正的新方法

```text
MLP torsion
+
new fixed-charge assignment model
```

尤其 charge 不只是：

```text
ML predicts RESP
```

而是：

> **learned, environment-aware, charge-conserving correction model optimized against several physically relevant objectives.**

---

# 25. 一句话总结整个方案

最终管线可以压缩成：

```text
                    molecular structure
                           │
                           ▼
                  graph / atom typing
                           │
          ┌────────────────┼──────────────────┐
          │                │                  │
          ▼                ▼                  ▼
      bond/angle        torsion             charge
          │                │                  │
  GFN-FF/xTB/QM       our MLP        new charge model
          │
 mSeminario /
 m2Seminario /
 DRIH /
 global fit
          │                │                  │
          └────────────────┼──────────────────┘
                           ▼
                      ForceFieldIR
                           │
             ┌─────────────┼─────────────┐
             ▼             ▼             ▼
          GROMACS        AMBER         OpenMM
```

也就是：

> **成熟算法负责基础 plumbing；GFN-FF 提供超便宜 Hessian；三种局部投影 + 一种 QMDFF-inspired 全局拟合负责 bond/angle；现有 MLP 接管 torsion；新模型接管 charge；最后通过统一 IR 同时输出 GROMACS、AMBER 和 OpenMM。**

这已经不再是 Sobtop clone，而是一个模块化的现代小分子经典力场参数化平台。

---

# 26. 当前明确不做的事情

第一版应主动排除：

- 从零重新训练完整 ML force field；
- 同时重新拟合所有 LJ；
- 自动混用 AMBER/CHARMM/OPLS 参数；
- 把 GFN-FF native topology 直接当最终 FF；
- 用 Hessian 拟合普通 flexible torsion，替代已有 MLP；
- 为了“统一”而强迫所有 backend 使用同一种文件表示；
- 第一版就解决 residue template / protein force-field generation。

范围控制很重要，否则项目会迅速膨胀成另一个完整 force-field ecosystem。

---

# 27. 参考与实现依据

1. Sobtop documentation：Sobtop 当前支持从 Gaussian、ORCA、xTB 和 CP2K Cartesian Hessian 派生 bonded force constants，并提供 Seminario、mSeminario、m2Seminario 和 DRIH 四种方法。本文档将原始 Seminario 视作 legacy，只把后三者作为核心实现。
2. xTB documentation：GFN-FF 支持 optimization、frequency/Hessian workflow；`--hess` 可对当前结构计算 Hessian，`--ohess` 可先优化再计算。
3. GFN-FF documentation：GFN-FF 是自动化 generic force field，包含 bond/angle/torsion 与非共价相互作用，并能生成 system-specific topology。
4. Behera S., Gapsys V., de Groot B. L. *Evaluation of the ABCG2 Charge Model in Protein–Ligand Binding Free-Energy Calculations*. J. Chem. Inf. Model. 2025, 65, 11505–11509. DOI: 10.1021/acs.jcim.5c02161. 该工作说明针对 hydration free energy 优化的 charge model 并不自动提高 protein–ligand RBFE，支持我们采用多目标 charge-model validation 的设计。
5. Grimme S. *A General Quantum Mechanically Derived Force Field (QMDFF) for Molecules and Condensed Phase Simulations*. J. Chem. Theory Comput. 2014, 10, 4497–4514. DOI: 10.1021/ct500573f. 从平衡结构 + QM Hessian + 部分电荷 + 共价键级自动生成 molecule-specific FF；bonded 力常数通过让 FF Hessian 拟合 QM Hessian 得到；可旋转键的 torsion 由自动构造的局部模型体系的完整 2π 转动曲线（extended Hückel 型最小价基 TB）拟合 Fourier 项，不可旋转 torsion 和 TB 失败的情况用 Hessian 拟合。§6.0 / §6.4 的直接先例。
