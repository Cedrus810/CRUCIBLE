# 通用小分子力场参数化编译器：设计草案

> 工作定位：不是简单复刻 Sobtop，而是把其已经成熟的拓扑/参数化流程重新拆成模块化、可测试、可扩展的 **force-field compiler**。  
> 核心新增部分：**GFN-FF / xTB / QM Hessian 后端 + 三种实用 Hessian 参数提取方法 + 已有 MLP 二面角模块 + 新型电荷模型 + GROMACS / AMBER / OpenMM 多后端输出**。

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

# 6. Hessian → bonded 参数：我们实际实现三种方法

Sobtop 文档目前列出了四种 Hessian 派生方案：原始 Seminario、modified Seminario、m2Seminario 和 DRIH。

但我们不需要把原始 Seminario 作为正式推荐路线。原始方法已经基本被 modified 方法取代。

因此我们的核心实现定为三种：

```text
1. mSeminario
2. m2Seminario
3. DRIH
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
mSeminario / m2Seminario / DRIH
```

把 GFN-FF 变成一个真正实用的 cheap Hessian backend。

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
 DRIH
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

> **成熟算法负责基础 plumbing；GFN-FF 提供超便宜 Hessian；三种实用 Hessian 映射方法负责 bond/angle；现有 MLP 接管 torsion；新模型接管 charge；最后通过统一 IR 同时输出 GROMACS、AMBER 和 OpenMM。**

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
