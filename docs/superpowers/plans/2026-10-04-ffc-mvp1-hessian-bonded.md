# ffcompiler MVP-1 施工计划：Hessian 子系统与 bond/angle 参数推导

> **给施工方：** 按任务顺序施工，每个任务先写测试再实现（TDD），任务完成即提交。本计划只给规格、接口、锁定的数学定义、必测用例和验收判据，不给实现代码。**T15 是硬停点**：理论文档交上来、代理审过之前，不得开始 T17–T19、T22 的实现。发现计划自相矛盾或数学定义有问题，停下来报告，不要自行改定义。

**目标：** 从 Cartesian Hessian（GFN-FF、GFN2-xTB、或现成的 Gaussian fchk）推出 bond 和 angle 参数，四种方法：mSeminario、m2Seminario、DRIH（局部投影）和 QMDFF-inspired 全局 Hessian 拟合（MVP-1b），结果写回 `ForceFieldIR`，并给出各方法对比报告（设计 §5–§8、§23 MVP-1）。

**架构：** Hessian provider 只负责输出统一单位、统一原子顺序的 `CartesianHessian`。参数推导方法只看 `CartesianHessian` 和 bond/angle 列表，不知道 Hessian 从哪来。驱动函数 `derive_bonded` 负责把结果写回 IR 并记录 provenance。Sobtop 二进制只在生成参考数据时调用一次，参考数据提交进仓库，测试不依赖 Sobtop。

**技术栈：** Python 3.12、numpy、scipy（`linalg`、`optimize.lsq_linear`）、RDKit（图等价类）、OpenMM Reference（FF Hessian、几何优化）、xtb 6.7.1、Sobtop 2026.1.16（只用来生成参考数据）、pytest。

**依据：** `forcefield_compiler_design.md` §5–§8、§6.4、§20 Level 1–2、§23；路线图 `2026-10-04-ffc-roadmap.md`（决策 D2–D6）。

**前置：** MVP-0 的 T0–T3、T5、T6 已完成（IR、校验、numpy 能量、prmtop 导入、OpenMM 导出）。T9 完成后，T23 才能做跨引擎复检。

---

## 1. 全局约束

MVP-0 计划 §1 的全部条目（环境、禁止装包、float64、Reference 平台、测试分层、`tests` 包遮蔽、日志、运行目录、种子、git）在这里同样适用。另外：

1. **外部程序**：
   - xtb：`/home/user/software/xtb-dist/bin/xtb`（6.7.1）。
   - Sobtop：把 `/home/user/software/sobtop_2026.1.16/` 整个复制到仓库内的 `tools/sobtop/`（加进 `.gitignore`），再 `chmod +x tools/sobtop/sobtop`。不得改动原目录。
   - 运行 Sobtop 必须设 `OMP_STACKSIZE=64M OMP_NUM_THREADS=4`。用默认线程数时会报 `OMP: Error #34 ... Resource temporarily unavailable`（2026-10-04 实测）。
2. **调用外部程序的测试**全部标 `@pytest.mark.external`。
3. **MVP-1 只推 bond 和 angle。** 二面角、improper 一律保留 IR 原值（GAFF2）。

---

## 2. 锁定的定义

### 2.1 单位与常数（`ffcompiler/units.py`）

| 常数 | 值 |
|---|---|
| `HARTREE_KJ_MOL` | 2625.499639479 |
| `BOHR_NM` | 0.0529177210544 |
| `ANGSTROM_NM` | 0.1 |
| `KCAL_KJ` | 4.184 |
| `C_CM_S` | 2.99792458e10 |

`CartesianHessian` 内部单位：坐标 nm、Hessian kJ/mol/nm²、梯度 kJ/mol/nm、能量 kJ/mol。Hartree/Bohr² → kJ/mol/nm² 的换算系数是 `HARTREE_KJ_MOL / BOHR_NM**2`。

### 2.2 谐振频率

- 质量加权：\(\tilde{\mathbf H}=\mathbf M^{-1/2}\mathbf H\mathbf M^{-1/2}\)，质量用 IR 的 `mass`（amu）；没有 IR 时用 RDKit 周期表的平均原子量。
- 投影掉平移和转动（以质心为原点，按质量加权构造 6 个基向量，线性分子 5 个，Gram–Schmidt 正交化），然后对角化。
- 特征值 λ（kJ mol⁻¹ nm⁻² amu⁻¹）→ 波数：\(\tilde\nu=\operatorname{sign}(\lambda)\sqrt{|\lambda|\times10^{24}}/(2\pi c)\)（cm⁻¹）。负值表示虚频。
- 返回 3N−6（线性分子 3N−5）个频率，升序排列。

### 2.3 Hessian 质量检查

| 量 | 定义 | 超过就报错 | 超过就警告 |
|---|---|---|---|
| 不对称度 | \(\max|H-H^\mathsf T|/\max|H|\) | 1e-3 | 1e-6（之后用 \(\tfrac12(H+H^\mathsf T)\)） |
| 平移残差 | \(\max_{i,a,b}|\sum_j H_{3i+a,3j+b}|/\max|H|\) | 1e-3 | 1e-5 |
| 驻点 | \(\max|g|\)（有梯度时） | 25 kJ/mol/nm（≈5×10⁻⁴ Ha/Bohr），除非 `allow_nonstationary=True` | — |
| 虚频 | 最小频率 | < −20 cm⁻¹ | −20…0 cm⁻¹ |

### 2.4 内坐标与 Wilson B

- 坐标类型：bond (i,j)、angle (i,j,k)（j 为顶点）、dihedral (i,j,k,l)（定义同 MVP-0 §2.1）。
- B 的行是坐标对 3N 个 Cartesian 分量的一阶导数，单位 1（bond，nm/nm）或 rad/nm。用解析公式：bond 和 angle 用 Wilson 公式；dihedral 用 Blondel–Karplus 公式。
- **线性键角**：θ > 175° 的键角（以及含这种键角的二面角）记为线性。所有方法遇到线性键角都抛 `LinearAngleError`，由驱动函数按 `linear_angle_policy` 处理（§2.7）。

### 2.5 DRIH（锁定）

1. 冗余内坐标集合 q = IR 全部 bond ∪ IR 全部非线性 angle ∪ 由键图生成的全部 proper dihedral 四元组（不含经过线性键角的四元组，不含 improper）。
2. 在参考几何上算 B（\(n_q\times3N\)），\(\mathbf G=\mathbf B\mathbf B^\mathsf T\)，用特征分解求伪逆 \(\mathbf G^-\)：特征值小于 \(10^{-8}\times\lambda_{\max}\) 的方向置零。
3. \(\mathbf A=\mathbf G^-\mathbf B\)，\(\mathbf H_q=\mathbf A\mathbf H_x\mathbf A^\mathsf T\)。**忽略梯度项**，所以必须先通过 §2.3 的驻点检查。
4. 每个 bond / angle 的力常数取 \((\mathbf H_q)_{tt}\)；平衡值取参考几何上的 \(q_t\)。单位自然就是 kJ/mol/nm² 和 kJ/mol/rad²，与 IR 的 ½k 约定一致。
5. 只取对角元。Sobtop 的 DRIH 会另外输出 bond–bond 耦合项，AMBER 型函数表示不了，我们丢弃。

T15 要确认以上与 Sobtop DRIH 的结果一致；不一致就停下报告（不得改定义去凑数）。

### 2.6 mSeminario、m2Seminario

- **mSeminario**：按 Allen et al., *J. Chem. Theory Comput.* 2018, 14, 274 的原文（决策 D6）。T15 负责把 bond 公式、angle 公式（含相邻键角的缩放因子）逐条誊写进 `docs/theory/hessian_methods.md`，带原文公式编号。验收以 Sobtop 的 mSeminario 输出为准。
- **m2Seminario**：卢天提出，没有发表的公式（决策 D2）。T15 写出下面两个候选定义，T18 拿 Sobtop 的 m2Seminario 结果逐一比对：
  - 记号：\(\mathbf K_{AB}=-\mathbf H_{AB}\)（3×3 原子间块），\(\mathbf u_{AB}\) 是从 A 指向 B 的单位向量。
  - **C1（直接投影）**：bond \(k_{AB}=\mathbf u_{AB}^\mathsf T\mathbf K_{AB}\mathbf u_{AB}\)；angle A-B-C 用 \(\mathbf u_N=\widehat{\mathbf u_{CB}\times\mathbf u_{AB}}\)、\(\mathbf u_{PA}=\mathbf u_N\times\mathbf u_{AB}\)、\(\mathbf u_{PC}=\mathbf u_{CB}\times\mathbf u_N\)，\(1/k_\theta=1/(R_{AB}^2\,\mathbf u_{PA}^\mathsf T\mathbf K_{AB}\mathbf u_{PA})+1/(R_{CB}^2\,\mathbf u_{PC}^\mathsf T\mathbf K_{CB}\mathbf u_{PC})\)。因为 \(\mathbf K_{BA}=\mathbf K_{AB}^\mathsf T\)，C1 天然满足 \(k_{AB}=k_{BA}\)，也没有本征向量取绝对值带来的方向歧义。
  - **C2**：C1 再乘 mSeminario 的相邻键角缩放因子（以 T15 誊写的原文公式为准）。
  - **命名规则**：哪个候选在全部 oracle 用例上都与 Sobtop m2Seminario 一致（§2.8 的容差），就用它，方法名叫 `m2seminario`。都不一致：把 C1 实现为 `symseminario`，报告比对结果，等用户定（D2）。
- 本节的 \(k\) 都是二阶导数，即 IR 的 ½k 约定。Sobtop 的 GROMACS 输出也是这个约定（SF6：param.txt 2928.942 kJ/mol/Å² = itp 2.928942E+05 kJ/mol/nm²）。

### 2.7 对称化与线性键角策略

- **原子等价类**：带颜色的键图的自同构轨道，颜色 = (`atomic_number`, `atom_type`, `formal_charge`)。只用图，不用 3D 几何。
- **项的等价类**：bond (i,j) 的键 = 无序对 {class(i), class(j)}，再要求存在一个图自同构把一个 bond 映到另一个（不能只比较类标签）。angle 同理，以顶点为中心的无序端点对。
- `symmetrize=True`（默认，决策 D5）：每个项等价类内，\(k\) 和平衡值都取算术平均；记录类内离散度 \((\max-\min)/|\text{mean}|\)。
- `linear_angle_policy`：
  - `"error"`（默认）：IR 里有线性键角就抛 `LinearAngleError`，信息里列出全部线性键角。
  - `"library"`：线性键角不推导，保留 IR 原参数，`meta["linear_angle"]=True`；从传给方法的 angle 列表里去掉。

### 2.8 Sobtop 参考数据与比对容差

- 已验证的 stdin 序列（2026-10-04，SF6 三种方法都跑通）：
  ```text
  <mol2 路径>
  -1
  <方法号：2=mSeminario，3=m2Seminario，4=DRIH>
  6
  <Hessian 路径：.fchk 或 xtb 的 hessian 文件>
  2
  0
  ```
  输出是当前目录下的 `param.txt`。格式：bond 行 `# 原子 原子 长度(Å) k(kJ/mol/Å²) k(kcal/mol/Å²)`，angle 行 `# 原子 原子 原子 角度(deg) k(kJ/mol/rad²) k(kcal/mol/rad²)`。
- 用 xtb Hessian 时，Sobtop 从 mol2 取几何，所以 mol2 必须写 `xtbopt.xyz` 的坐标。
- 比对容差：\(|\Delta k|\le10^{-4}|k_\text{Sobtop}|+5\times10^{-3}\)（param.txt 的单位）；键长 ≤ 2×10⁻⁶ Å；键角 ≤ 2×10⁻⁴ deg。
- **SF6 不进 Gate**，只作诊断：它有 180° 键角。实测 Sobtop 在 SF6 上即使用 m2Seminario/DRIH，对称等价的 90° 键角也给出不同的力常数（608–692、636/653 kJ/mol/rad²）。

### 2.9 全局 Hessian 拟合（MVP-1b，锁定）

按设计 §6.4。记待推导的 bond/angle 项为 t，等价类为 c（`symmetrize=False` 时每项自成一类）。

1. **参考几何** \(\mathbf x_0\) = Hessian 的几何；平衡值 \(q_t^0=q_t(\mathbf x_0)\)（对称化时类内平均）。
2. **固定部分** \(\mathbf H^\text{fixed}\)：IR 去掉全部待推导项后（其余 bond/angle、全部 torsion、improper、LJ、Coulomb、1-4 原样保留），在 \(\mathbf x_0\) 处的 FF Hessian。用 OpenMM Reference 的解析力做中心差分，步长 \(h=10^{-5}\) nm，再取 \(\tfrac12(H+H^\mathsf T)\)。
3. **设计矩阵**：\(\mathbf G_c=\sum_{t\in c}\mathbf b_t\mathbf b_t^\mathsf T\)（\(\mathbf b_t\) 是 B 的行）。以先验为尺度，令 \(k_c=s_c\,k_c^\text{prior}\)，列向量 \(\mathbf D_{:,c}=k_c^\text{prior}\operatorname{vec}(\mathbf G_c)\)，目标向量 \(\mathbf y=\operatorname{vec}(\mathbf H^\text{ref}-\mathbf H^\text{fixed})\)。vec 取全部 \(9N^2\) 个元素，权重为 1（决策 D4）。
4. **求解**：\(\min_{\mathbf s\ge0}\ \|\mathbf y-\mathbf D\mathbf s\|^2+\lambda\|\mathbf s-\mathbf 1\|^2\)，\(\lambda=\alpha\cdot\operatorname{tr}(\mathbf D^\mathsf T\mathbf D)/n_c\)，默认 \(\alpha=10^{-3}\)。用 `scipy.optimize.lsq_linear`（把正则项增广进最小二乘，下界 0）。
5. **先验** \(k^\text{prior}\)：默认 mSeminario（同一 Hessian、同样的对称化）。
6. **报告**：α、λ、\(\operatorname{cond}_2(\mathbf D)\)、拟合后相对残差 \(\|\mathbf y-\mathbf D\mathbf s\|/\|\mathbf y\|\)、先验（\(\mathbf s=\mathbf 1\)）的相对残差、落在下界的类数、每类的 \(s_c\)、拟合后 FF 在 \(\mathbf x_0\) 处的最大力（"平衡位置不一致"的量化）。

---

## 3. 文件结构

```text
src/ffcompiler/
├── units.py                      # T11
├── hessian/
│   ├── __init__.py
│   ├── base.py                   # T11：CartesianHessian、检查、频率
│   ├── gaussian.py               # T12
│   └── xtb.py                    # T13
├── bonded/
│   ├── __init__.py
│   ├── internal.py               # T14：内坐标、Wilson B、LinearAngleError
│   ├── mseminario.py             # T17
│   ├── m2seminario.py            # T18（或 symseminario，见 §2.6）
│   ├── drih.py                   # T19
│   ├── symmetry.py               # T20
│   ├── assign.py                 # T21：derive_bonded
│   └── global_fit.py             # T22
├── validation/
│   └── vibrations.py             # T22：ff_hessian、频率对比、FF 几何优化
docs/theory/hessian_methods.md    # T15
scripts/make_sobtop_refs.py       # T16
scripts/gate_mvp1.py              # T23
tests/fixtures/hessian/           # T12/T13/T16：fchk、mol2、xtb 输出、Sobtop param.txt
tests/unit/test_hessian_*.py、test_internal.py、test_mseminario.py、...
```

---

## 4. Review Focus

1. **原子顺序不一致，或 xtb 优化时断键/成键**：Hessian 的原子序列与 IR 不同，或优化后的几何与 IR 拓扑不符。必须报错，不能拿错位的 Hessian 继续推参数。→ T13、T21。
2. **线性键角**：乙腈（C-C≡N）、SF6 的 180° 角。默认报错，`"library"` 策略下保留原参数。→ T14、T17–T19、T21。
3. **非驻点几何**：用户给的 Hessian 不在极小点（fchk 有梯度）。DRIH 和全局拟合忽略梯度项，必须拦住。→ T11、T19、T22。
4. **参考 Hessian 是鞍点**：有虚频。→ T11、T13。
5. **带电分子**：xtb 必须收到正确的 `--chrg`；电荷从 IR 的 `total_charge` 取。→ T13（醋酸根、甲铵）。

---

## 5. 任务

### T11：单位、`CartesianHessian`、检查与频率

**Files:**
- Create: `src/ffcompiler/units.py`、`src/ffcompiler/hessian/__init__.py`、`src/ffcompiler/hessian/base.py`、`tests/unit/test_hessian_base.py`

**Interfaces:**
- Produces:
  ```python
  @dataclass(frozen=True)
  class CartesianHessian:
      atomic_numbers: np.ndarray          # (N,) int
      positions_nm: np.ndarray            # (N,3)
      hessian: np.ndarray                 # (3N,3N) kJ/mol/nm²
      gradient: np.ndarray | None         # (3N,) kJ/mol/nm
      energy_kj: float | None
      provider: str                       # "gaussian-fchk" | "xtb"
      method: str                         # 例如 "gfnff"、"gfn2"、fchk 里的方法名
      source_file: str | None
      source_sha256: str | None
      program_version: str | None

  @dataclass(frozen=True)
  class HessianDiagnostics:
      asymmetry_rel: float; translation_residual_rel: float
      max_gradient: float | None; frequencies_cm: np.ndarray
      n_imaginary: int; warnings: tuple[str, ...]

  class HessianQualityError(ValueError): ...

  def check_hessian(h: CartesianHessian, masses: np.ndarray | None = None, *, allow_nonstationary: bool = False) -> HessianDiagnostics
  def symmetrized(h: CartesianHessian) -> CartesianHessian
  def harmonic_frequencies(hessian: np.ndarray, positions_nm: np.ndarray, masses: np.ndarray) -> np.ndarray
  ```
- `check_hessian` 按 §2.3 报错或警告。

**必测用例**：
1. 单位：`HARTREE_KJ_MOL`、`BOHR_NM` 与 `scipy.constants` 推出的值相对差 ≤ 1e-9。
2. 双原子谐振子：两原子相距 0.1 nm，键 k=5×10⁵ kJ/mol/nm²，质量 1 和 19 amu。解析 Hessian 的频率 = \(\sqrt{k/\mu\times10^{24}}/(2\pi c)\)，误差 ≤ 1e-6 cm⁻¹；线性分子返回 1 个频率。
3. 平移残差：对一个合法 Hessian 加一个破坏平移不变性的扰动 → `HessianQualityError`，信息含 "translation"。
4. 不对称度 1e-4 → 报错；1e-8 → 通过，`symmetrized` 后严格对称。
5. 带梯度 max|g| = 100 kJ/mol/nm → 默认报错（信息含 "stationary"）；`allow_nonstationary=True` 时只记录。（Review Focus 3）
6. 构造一个有 −50 cm⁻¹ 虚频的 Hessian（在某个振动模式上取负曲率）→ 报错；−5 cm⁻¹ → 警告。（Review Focus 4）

---

### T12：Gaussian fchk 读取器

**Files:**
- Create: `src/ffcompiler/hessian/gaussian.py`、`tests/fixtures/hessian/sobtop_examples/{SF6.fchk, SF6.mol2, Methyl_benzoate.fchk, Methyl_benzoate.mol2}`（从 Sobtop examples 复制）、`tests/unit/test_hessian_gaussian.py`

**Interfaces:**
- Produces: `read_fchk(path: Path) -> CartesianHessian`

**规格**：
- 读取字段：`Atomic numbers`（I）、`Current cartesian coordinates`（R，Bohr）、`Cartesian Force Constants`（R，下三角按行展开：(0,0),(1,0),(1,1),(2,0),…，共 3N(3N+1)/2 个，Ha/Bohr²）、`Cartesian Gradient`（R，Ha/Bohr，没有就为 None）、`Total Energy`。
- 缺 `Cartesian Force Constants` → 报错，信息提示"不是 freq 任务的 fchk"。
- `method` 取 fchk 第二行的方法名。

**必测用例**：
1. SF6.fchk：7 个原子；S–F 距离 = 1.578599 Å（与 Sobtop param.txt 一致，误差 ≤ 1e-6 Å）；`check_hessian` 通过；频率中最小值 > 0。
2. Methyl_benzoate.fchk：原子数与 mol2 一致；`check_hessian` 通过。
3. 下三角展开顺序：手写一个 2 原子的小 fchk 文本（测试里内联），各元素位置正确，且结果对称。
4. 删掉 `Cartesian Force Constants` 段的 fchk → 报错。

---

### T13：xtb 运行器与读取器

**Files:**
- Create: `src/ffcompiler/hessian/xtb.py`、`tests/fixtures/hessian/xtb_ethanol_gfn2/`（一次真实运行的输出：`hessian`、`xtbopt.xyz`、`vibspectrum`、`xtb.out`）、`tests/unit/test_hessian_xtb.py`

**Interfaces:**
- Produces:
  ```python
  XTB_BINARY = Path("/home/user/software/xtb-dist/bin/xtb")

  def read_xtb_hessian(hessian_file: Path, xyz_file: Path, *, method: str, version: str | None = None) -> CartesianHessian
  def read_xtb_vibspectrum(path: Path) -> np.ndarray     # cm⁻¹，去掉 6 个（或 5 个）零模

  class XtbError(RuntimeError): ...

  def run_xtb_hessian(atomic_numbers: np.ndarray, positions_nm: np.ndarray, *, method: Literal["gfnff", "gfn2"],
                      charge: int, workdir: Path, threads: int = 4, binary: Path = XTB_BINARY) -> CartesianHessian
  ```

**规格**：
- 命令：`xtb input.xyz --gfnff --ohess vtight --chrg <Q> --uhf 0`；GFN2 用 `--gfn 2` 代替 `--gfnff`。环境变量 `OMP_NUM_THREADS=<threads>`、`OMP_STACKSIZE=500M`。stdout 写进 `workdir/xtb.out`。
- `hessian` 文件：`$hessian` 之后到 `$end`（或文件尾）的全部浮点数，按行主序共 (3N)² 个，单位 Ha/Bohr²。个数不对就报错。
- 几何取 `xtbopt.xyz`（Å）。
- 成功条件：退出码 0、`xtb.out` 含 `GEOMETRY OPTIMIZATION CONVERGED`、`hessian` 和 `xtbopt.xyz` 存在。任何一条不满足都抛 `XtbError`，信息附上 `xtb.out` 最后 30 行。
- `gradient=None`（`--ohess` 先优化到 vtight，视为驻点；`meta` 记 `"stationary_by": "xtb ohess vtight"`）。
- 开壳层：MVP-1 不支持，`uhf` 固定为 0。

**必测用例**：
1. （非 external）读 `xtb_ethanol_gfn2` 固件：由 `harmonic_frequencies` 算出的频率与 `vibspectrum` 逐个对比，max 差 ≤ 1 cm⁻¹。这一条同时验证格式、单位和质量。
2. （非 external）`hessian` 文件里少一个数 → 报错。
3. （external）ethanol 跑 GFN-FF 和 GFN2 都成功；`check_hessian` 通过；最小频率 > −20 cm⁻¹。
4. （external）醋酸根 `charge=-1`、甲铵 `charge=+1`：`xtb.out` 里的总电荷与之一致（解析 `xtb.out` 的电荷行断言）（Review Focus 5）。
5. （external）故意传错电荷使体系为奇电子（醋酸根给 `charge=0`）→ `XtbError`，信息提示电子数为奇数或 uhf 不匹配。
6. （非 external）`xtb.out` 里没有 `GEOMETRY OPTIMIZATION CONVERGED`（用改过的固件文本）→ `XtbError`。

---

### T14：内坐标与 Wilson B

**Files:**
- Create: `src/ffcompiler/bonded/__init__.py`、`src/ffcompiler/bonded/internal.py`、`tests/unit/test_internal.py`

**Interfaces:**
- Produces:
  ```python
  InternalCoord = tuple[Literal["bond", "angle", "dihedral"], tuple[int, ...]]
  LINEAR_ANGLE_DEG = 175.0

  class LinearAngleError(ValueError):
      angles: list[tuple[int, int, int]]

  def coord_value(positions_nm: np.ndarray, coord: InternalCoord) -> float      # nm 或 rad
  def wilson_b(positions_nm: np.ndarray, coords: Sequence[InternalCoord]) -> np.ndarray   # (n_q, 3N)
  def linear_angles(positions_nm: np.ndarray, angles: Sequence[tuple[int, int, int]]) -> list[tuple[int, int, int]]
  def proper_quartets(bond_graph: Mapping[int, frozenset[int]]) -> list[tuple[int, int, int, int]]   # 去重：(i,j,k,l) 与 (l,k,j,i) 只留 i<l 的那个
  ```
- `wilson_b` 遇到线性键角（或经过线性键角的二面角）抛 `LinearAngleError`，列出全部线性键角。

**必测用例**：
1. 对 ethanol 固件几何，每一行 B 与 `coord_value` 的中心差分（步长 1e-6 nm）一致：max 相对误差 ≤ 1e-7。
2. 平移不变性：B 的每一行在三个平移方向上的和 ≤ 1e-10。
3. 乙腈固件：`linear_angles` 返回 C-C≡N 那个角；把它放进 `wilson_b` → `LinearAngleError`（Review Focus 2）。
4. `proper_quartets`：正丁烷骨架只有 1 个；苯（含 H）是 24 个；没有重复。
5. 二面角在 ±π 附近（反式构型）时 B 行仍与中心差分一致（差分时要处理 2π 跳变）。

---

### T15：理论文档（硬停点）

**Files:**
- Create: `docs/theory/hessian_methods.md`

**内容要求**：
1. **mSeminario**：从 Allen et al. 2018 原文誊写 bond 和 angle 公式（含相邻键角缩放因子），标原文公式编号；说明原文是否对 \(k_{AB}\)、\(k_{BA}\) 取平均，以及 \(\mathbf K_{AB}\) 非对称时本征值如何取实部。原文由用户提供，或获准联网查阅（D6）。
2. **m2Seminario**：§2.6 的 C1、C2 两个候选，写出完整公式，并说明各自如何满足"等价键角相同、\(k_{ij}=k_{ji}\)"。
3. **DRIH**：§2.5 的推导，说明忽略梯度项的前提，以及冗余坐标集合的选择会影响对角元。
4. **全局拟合**：§2.9 的推导，包括 \(\mathbf b\mathbf b^\mathsf T\) 形式只在 \(q^0=q(\mathbf x_0)\) 时成立的原因，以及线性最小二乘为什么不需要 LM。
5. **单位约定**：每种方法输出的 k 都是二阶导数（½k 约定），并写出与 Sobtop param.txt、GROMACS itp、AMBER prmtop 之间的换算。

**完成标准**：文档提交后**停工**，通知 代理审阅。审阅通过前不得开始 T17–T19、T22；T16、T20 可以继续。

---

### T16：Sobtop 参考数据

**Files:**
- Create: `scripts/make_sobtop_refs.py`、`tests/fixtures/hessian/sobtop_refs/<case>/{input.mol2, hessian 或 .fchk, param_m2.txt, param_m3.txt, param_m4.txt}`、`tests/fixtures/hessian/sobtop_refs/manifest.json`、`src/ffcompiler/bonded/sobtop_ref.py`、`tests/unit/test_sobtop_ref.py`

**Interfaces:**
- Produces:
  - `read_sobtop_params(path: Path) -> SobtopParams`：dataclass，`bonds: dict[tuple[int,int], tuple[float, float]]`（键长 Å，k kJ/mol/Å²），`angles: dict[tuple[int,int,int], tuple[float, float]]`（角度 deg，k kJ/mol/rad²）。原子号从 param.txt 的 1 起算转成 0 起算。

**参考用例**（case 名即目录名）：
- `sf6_fchk`、`methyl_benzoate_fchk`：Sobtop 自带的 fchk + mol2。
- `<fixture>_gfn2`：MVP-0 的 12 个 GAFF2 固件中除乙腈以外的 11 个。先用 `run_xtb_hessian(method="gfn2")` 得到 Hessian，再把 `xtbopt.xyz` 的坐标连同 IR 的键连写成 mol2 交给 Sobtop。原子顺序与 IR 相同。
- 每个 case 对方法 2、3、4 各跑一次（stdin 序列见 §2.8）。

**规格**：
- 脚本把 Sobtop 复制到 `tools/sobtop/`（如果还没有），设置 §1 的环境变量，在每个 case 的临时子目录里运行，把 `param.txt` 改名存好。
- Sobtop 退出码非 0、没生成 `param.txt`、或 param.txt 里的 bond/angle 数与 mol2 键连推出的数不一致 → 脚本报错退出。
- `manifest.json` 记录 Sobtop 版本、每个 case 的输入 sha256、xtb 版本。
- 参考数据提交进仓库。

**必测用例**（非 external）：
1. `read_sobtop_params` 读 `sf6_fchk/param_m4.txt`：6 个 bond，k 全部 = 2928.942；15 个 angle。
2. 每个 case 的 bond 数、angle 数与 mol2 键连推出的一致。

---

### T17：mSeminario

**Files:**
- Create: `src/ffcompiler/bonded/mseminario.py`、`tests/unit/test_mseminario.py`

**Interfaces:**
- Consumes: T11、T14；T15 通过审阅的公式。
- Produces:
  ```python
  @dataclass(frozen=True)
  class DerivedTerm:
      atoms: tuple[int, ...]       # bond (i,j) 或 angle (i,j,k)，与输入顺序相同
      value: float                 # r0 nm 或 theta0 rad（参考几何上的值）
      k: float                     # kJ/mol/nm² 或 kJ/mol/rad²（½k 约定）

  def mseminario(h: CartesianHessian, bonds: Sequence[tuple[int, int]],
                 angles: Sequence[tuple[int, int, int]]) -> list[DerivedTerm]
  ```
  返回顺序：先全部 bond（按输入顺序），再全部 angle（按输入顺序）。输入里有线性键角 → `LinearAngleError`。

**必测用例**：
1. **Sobtop oracle**：`sobtop_refs` 里除 `sf6_fchk` 外的每个 case，全部 bond、angle 与 `param_m2.txt` 一致（§2.8 容差）。
2. 双原子：解析 Hessian 来自单根弹簧 k → 推出的 bond k 等于输入 k（相对误差 ≤ 1e-10）。
3. 结果与原子编号无关：对 ethanol 的 Hessian 做一次原子重排（同时重排 Hessian 的行列和项的下标），推出的参数按映射对应后相同。

---

### T18：m2Seminario（或 symseminario）

**Files:**
- Create: `src/ffcompiler/bonded/m2seminario.py`、`tests/unit/test_m2seminario.py`、`docs/decisions/<下一个编号>-m2seminario-naming.md`

**Interfaces:**
- Produces: `m2seminario(h, bonds, angles) -> list[DerivedTerm]`，签名、返回顺序与 `mseminario` 相同。如果按 §2.6 的命名规则落到 `symseminario`，函数和模块改名为 `symseminario`，其余不变。

**步骤要求**：
1. 先实现 C1、C2 两个候选（模块内私有函数），对全部 oracle case 与 `param_m3.txt` 比对，把比对表写进决策记录。
2. 按 §2.6 的命名规则确定公开的方法名，并在决策记录里写明。

**必测用例**：
1. 若命名为 `m2seminario`：全部 oracle case（除 `sf6_fchk`）与 `param_m3.txt` 一致。若为 `symseminario`：该测试改为 `xfail(strict=True)`，原因写 "D2 待定"，并保留比对数据。
2. **方向对称**：对每个 bond，交换输入顺序 (i,j)→(j,i)，k 完全相同（≤ 1e-12 相对误差）。
3. **等价项相同**：用解析构造的严格 D6h 苯几何（C-C 0.139 nm、C-H 0.108 nm，全部 120°），以及合成 Hessian \(\mathbf H=\sum_t k_t\mathbf b_t\mathbf b_t^\mathsf T\)（T14 的 `wilson_b`；全部 C-C bond 同一个 k，全部 C-H 同一个 k，C-C-C 和 C-C-H 各一个 k）。推导结果中 6 条 C-C bond、6 条 C-H bond、6 个 C-C-C、12 个 C-C-H 各自的 k 的相对离散度 ≤ 1e-8。等价集合在测试里手写，不依赖 T20。
4. 与 mSeminario 的 T17 用例 3 相同的原子重排不变性。

---

### T19：DRIH

**Files:**
- Create: `src/ffcompiler/bonded/drih.py`、`tests/unit/test_drih.py`

**Interfaces:**
- Produces: `drih(h, bonds, angles, *, dihedrals: Sequence[tuple[int,int,int,int]], allow_nonstationary: bool = False) -> list[DerivedTerm]`。`dihedrals` 只用来组成冗余坐标集合（§2.5 第 1 条），不返回。先调用 `check_hessian`（驻点检查按 `allow_nonstationary`）。

**必测用例**：
1. **Sobtop oracle**：全部 oracle case（除 `sf6_fchk`）与 `param_m4.txt` 一致。
2. **SF6 诊断**（不进 Gate，但测试必须存在）：`sf6_fchk` 的 bond 与 Sobtop 一致；用 `linear_angle_policy` 的等效做法（去掉 180° 键角）跑一遍，把 90° 键角的结果与 Sobtop 的差异写进测试输出，不做断言。
3. **无冗余时精确恢复**：合成 Hessian \(\mathbf H=\sum_t k_t\mathbf b_t\mathbf b_t^\mathsf T\)（T14 的 `wilson_b`），坐标集合恰好 3N−6 个且互相独立：
   - 水：2 个 bond + 1 个 angle；
   - H₂O₂（非平面构型，二面角约 110°）：3 个 bond + 2 个 angle + 1 个 dihedral。
   DRIH 推回的 bond/angle k 与输入的相对误差 ≤ 1e-8。
4. **有冗余时不能精确恢复（记录而不断言）**：对 ethanol 用同样的合成 Hessian，把 DRIH 推回值与输入值的最大相对偏差写进测试输出。sp3 中心有 6 个键角，却只有 5 个独立自由度，投影 \(\mathbf G^-\mathbf G\) 不是单位阵，对角元会混进别的坐标，这是 DRIH 的固有性质，不是 bug。
5. 带梯度的非驻点 Hessian（fchk 的梯度人为放大）→ 默认报错（Review Focus 3）。

---

### T20：对称等价类

**Files:**
- Create: `src/ffcompiler/bonded/symmetry.py`、`tests/unit/test_symmetry.py`

**Interfaces:**
- Produces:
  ```python
  def atom_classes(ir: ForceFieldIR) -> tuple[int, ...]                                   # 每个原子的类号
  def term_classes(ir: ForceFieldIR, terms: Sequence[tuple[int, ...]]) -> tuple[int, ...]  # bond 或 angle 的类号
  def symmetrize(terms: Sequence[DerivedTerm], classes: Sequence[int]) -> tuple[list[DerivedTerm], dict[int, float]]
      # 返回对称化后的项，以及每类 k 的相对离散度 (max-min)/|mean|
  ```
- 规格见 §2.7。实现可以用 RDKit（同位素编码颜色后求规范秩）或 networkx 的同构匹配，但必须满足下面的测试。

**必测用例**：
1. 苯：6 个 C 同类、6 个 H 同类；C-C bond 一类、C-H bond 一类；C-C-C、C-C-H 各一类。
2. 乙醇：CH3 的 3 个 H 同类，CH2 的 2 个 H 同类，OH 的 H 单独一类。
3. 降冰片烷：两个桥头 C 同类。
4. 带同样颜色但不等价的原子不得归为一类：用正戊烷骨架，C1/C5 同类、C2/C4 同类、C3 单独。
5. `symmetrize`：一类内 k 为 [1, 2, 3] → 都变为 2，离散度 = 1.0；单元素类不变，离散度 0。

---

### T21：驱动函数 `derive_bonded` 与 CLI

**Files:**
- Create: `src/ffcompiler/bonded/assign.py`、`tests/unit/test_assign.py`
- Modify: `src/ffcompiler/cli.py`（加子命令）

**Interfaces:**
- Consumes: T11–T14、T17–T20；MVP-0 的 `ForceFieldIR.replace_terms`、`validate`。
- Produces:
  ```python
  BondedMethod = Literal["mseminario", "m2seminario", "drih"]   # 若 D2 落到 symseminario，这里同步改名

  @dataclass(frozen=True)
  class BondedReport:
      method: str; hessian_provider: str; hessian_method: str; hessian_sha256: str | None
      diagnostics: HessianDiagnostics
      linear_angles_kept: tuple[tuple[int, int, int], ...]
      class_spread: dict[str, float]          # 项的标签 → 离散度
      max_bond_ratio: float                   # 优化后/输入 键长比的最大偏离

  class TopologyMismatchError(ValueError): ...

  def derive_bonded(ir: ForceFieldIR, h: CartesianHessian, *, method: BondedMethod,
                    symmetrize: bool = True, linear_angle_policy: Literal["error", "library"] = "error",
                    allow_nonstationary: bool = False,
                    input_positions_nm: np.ndarray | None = None) -> tuple[ForceFieldIR, np.ndarray, BondedReport]
  ```
  返回新 IR、Hessian 的几何（nm）和报告。

**规格**：
1. **原子顺序**：`h.atomic_numbers` 必须与 IR 的原子序列完全相同，否则 `TopologyMismatchError`。不做自动重排。
2. **拓扑一致性**（Review Focus 1）：在 Hessian 几何上，
   - 每个 IR bond 的长度与 `input_positions_nm` 上的长度之比在 [0.8, 1.25] 内（没给输入坐标时跳过这一条）；
   - 每个非 bond 原子对的距离 ≥ 1.1×(两原子 RDKit 共价半径之和)。
   违反任何一条就 `TopologyMismatchError`，列出原子对和距离。
3. 全部 bond 和非线性 angle 换成推导值，`source=<method>`，`meta` 含 `hessian_provider`、`hessian_method`、`hessian_sha256`、`symmetrized`、`class_spread`。线性键角按 §2.7 的策略处理。
4. torsion、improper、非键、exclusions、pairs14 不变。新 IR 的 `provenance["bonded"]` 记录方法、参数和 `BondedReport` 的摘要。
5. 结尾调用 `validate`。
6. CLI：`ffcompiler derive-bonded --prmtop P --inpcrd C --sdf S --hessian-provider {gfnff,gfn2,fchk} [--hessian-file F] --method M [--no-symmetrize] [--linear-angles {error,library}] --output {gromacs,amber,openmm}... --outdir D`。xtb provider 在 `outdir/xtb/` 下运行；输出用 Hessian 几何作坐标。

**必测用例**：
1. ethanol + `xtb_ethanol_gfn2` 固件 Hessian，三种方法都能跑；新 IR 通过 `validate`；torsion 项与原 IR 完全相同；全部 bond/angle 的 `source` 等于方法名。
2. Hessian 原子顺序打乱 → `TopologyMismatchError`（Review Focus 1）。
3. 把 Hessian 几何里某条 C-O 键拉长到 1.5 倍 → `TopologyMismatchError`，信息含该原子对（Review Focus 1）。
4. 乙腈：默认报 `LinearAngleError`；`"library"` 时线性角保留 GAFF2 参数且 `meta["linear_angle"]=True`，其余项正常推导（Review Focus 2）。
5. `symmetrize=True` 时乙醇 CH3 三个 C-H 的 k 完全相同；`False` 时允许不同，`class_spread` 记录离散度。
6. （external）CLI 对 ethanol 用 `gfnff` 跑通，三种格式文件都生成。

---

### T22：FF Hessian 工具与全局 Hessian 拟合（MVP-1b）

**Files:**
- Create: `src/ffcompiler/validation/vibrations.py`、`src/ffcompiler/bonded/global_fit.py`、`tests/unit/test_vibrations.py`、`tests/unit/test_global_fit.py`

**Interfaces:**
- Consumes: T14、T17、T20、T21；MVP-0 的 `to_openmm_system`。
- Produces:
  ```python
  # validation/vibrations.py
  def ff_hessian(ir: ForceFieldIR, positions_nm: np.ndarray, *, exclude_terms: Literal["none", "bonds_angles"] = "none",
                 step_nm: float = 1e-5) -> np.ndarray            # (3N,3N) kJ/mol/nm²，已对称化
  def minimize_ff(ir: ForceFieldIR, positions_nm: np.ndarray, *, tolerance_kj_nm: float = 1e-3) -> np.ndarray
  def kabsch_rmsd(a: np.ndarray, b: np.ndarray) -> float            # nm

  # bonded/global_fit.py
  @dataclass(frozen=True)
  class GlobalFitReport:
      alpha: float; lam: float; cond: float
      rel_residual_fit: float; rel_residual_prior: float
      n_at_bound: int; scale_by_class: dict[str, float]
      max_force_at_x0: float                                        # kJ/mol/nm

  def global_fit(ir: ForceFieldIR, h: CartesianHessian, *, prior_method: BondedMethod = "mseminario",
                 prior_ir: ForceFieldIR | None = None,
                 alpha: float = 1e-3, symmetrize: bool = True,
                 linear_angle_policy: Literal["error", "library"] = "error",
                 allow_nonstationary: bool = False, fd_step_nm: float = 1e-5,
                 input_positions_nm: np.ndarray | None = None,
                 subtract_fixed: bool = True) -> tuple[ForceFieldIR, np.ndarray, GlobalFitReport]
  ```
- 数学定义见 §2.9。
- `prior_ir` 不为 None 时，先验 k 直接取它的 bond/angle 力常数（必须与 `ir` 同拓扑），忽略 `prior_method`。
- `subtract_fixed=False` 只用于诊断（T22 用例 3）：令 \(\mathbf H^\text{fixed}=0\)。报告和 `meta` 里必须记下这个值。`exclude_terms="bonds_angles"` 只去掉待推导的 bond/angle（线性键角在 `"library"` 策略下保留在固定部分）。
- 新 IR 中拟合项的 `source="global_fit"`，`meta` 含 `prior_method`、`alpha`、`scale`。

**必测用例**：
1. `ff_hessian`：只有 bond+angle 项的 IR，在 \(q^0=q(\mathbf x_0)\) 的几何上，FD Hessian 与解析的 \(\sum_t k_t\mathbf b_t\mathbf b_t^\mathsf T\) 一致，相对误差 ≤ 1e-6。
2. **自洽恢复**：取 ethanol 的 GAFF2 IR，把全部 bond/angle 平衡值设成固件几何 \(\mathbf x_0\) 上的值，构造 \(\mathbf H^\text{ref}\)=该 IR 的完整 FF Hessian（含 torsion 和非键）。通过 `prior_ir` 传入先验：真值乘 U(0.5, 1.5)（种子 20261004）。`alpha=0`、`symmetrize=False` 时，拟合出的 k 与真值的相对误差 ≤ 1e-5。`symmetrize=True` 时同样成立（GAFF2 的等价项参数本来就相同）。这个构造是精确的：\(q^0=q(\mathbf x_0)\) 时，\(\mathbf H^\text{ref}-\mathbf H^\text{fixed}\) 严格等于 \(\sum_t k_t\mathbf b_t\mathbf b_t^\mathsf T\)。
3. **扣掉固定部分的必要性**：同一构造下用 `subtract_fixed=False`，恢复误差明显变大（> 1e-2）。这条说明设计 §6.4.1 的扣除确实起作用。
4. **下界**：构造一个先验让某一类的无约束解为负，拟合后该类 \(s_c=0\)，`n_at_bound ≥ 1`。
5. **报告字段**：`cond`、两个相对残差、`max_force_at_x0` 都是有限正数；`rel_residual_fit ≤ rel_residual_prior`。
6. 非驻点 Hessian → 默认报错（Review Focus 3）。

---

### T23：评估脚本与 Gate G1

**Files:**
- Create: `scripts/gate_mvp1.py`、`tests/unit/test_gate_mvp1_report.py`

**Interfaces:**
- `scripts/gate_mvp1.py [--providers gfnff gfn2] [--jobs N]`：在 `RunContext("gate_mvp1")` 下运行，写 `report.md`、`metrics.csv`、`acceptance.json`。

**评估矩阵**：
- 分子：MVP-0 的 12 个固件（乙腈用 `linear_angle_policy="library"`）。
- Hessian：`gfnff`、`gfn2`（各固件现场跑 xtb，按核数分线程并行；结果存进运行目录）。
- 方法：`gaff2`（原 IR，不推导，作基线）、`mseminario`、`m2seminario`（或 `symseminario`）、`drih`、`global_fit`。

**每个 (分子, Hessian, 方法) 的指标**（写进 `metrics.csv`）：
1. **频率**：FF Hessian 在参考几何上的频率 vs 参考 Hessian 的频率，两者都升序、逐个配对。指标：RMSE、MAE、最大绝对误差；再按参考频率 < 1000 和 ≥ 1000 cm⁻¹ 分开算 RMSE。
2. **几何**：从参考几何出发做 FF 能量最小化，指标：全原子 Kabsch RMSD（pm）、bond 长度 MAE（pm）、angle MAE（deg）。
3. 对称化离散度的最大值。
4. 全局拟合：`GlobalFitReport` 全部字段。
5. 推导后的 IR 跑 MVP-0 的 `compare`，记录 `passed`。

`report.md`：每个 Hessian 级别一张"方法 × 指标"的汇总表（12 个分子的中位数和最大值），外加 Sobtop 比对表和全部失败项清单。

**`acceptance.json`**：`gate: "G1"`、`passed`、`sobtop_match`（每个方法：一致的 case 数 / 总数）、`m2_naming`（`"m2seminario"` 或 `"symseminario"`）、`global_fit_self_consistency`（T22 用例 2 的误差）、`cross_engine_all_passed`、`failures`（列出每一个失败的组合及原因，不得静默丢弃）、`pytest_summary`。

**必测用例**：用两个假的指标行跑报告生成函数，`report.md` 表头固定，`acceptance.json` 字段齐全；有一个组合抛异常时，它出现在 `failures` 里且 `passed=false`。

---

## 6. Gate G1 验收判据

施工方提交 `jobs/gate_mvp1_<date>/` 全部材料和全套 pytest 日志；代理只审阅、不重跑。

| # | 判据 |
|---|---|
| G1.1 | mSeminario、DRIH 在全部 oracle case（除 `sf6_fchk`）上与 Sobtop 一致（§2.8 容差） |
| G1.2 | m2Seminario：要么全部一致并命名为 `m2seminario`，要么按 D2 命名为 `symseminario`，且决策记录里有完整比对表 |
| G1.3 | 全局拟合自洽恢复误差 ≤ 1e-5（T22 用例 2），扣除必要性测试（用例 3）通过 |
| G1.4 | 全部推导后的 IR 通过 MVP-0 的跨引擎比较 |
| G1.5 | 评估矩阵每个组合都有结果或列在 `failures` 里并说明原因；没有 NaN |
| G1.6 | Review Focus 1–5 对应的测试存在并通过（报告列出测试名） |
| G1.7 | pytest 全部通过；skip 只允许是 `gmx` 不存在 |
| G1.8 | T15 理论文档已经 代理审阅通过；§2 的定义没有被改动，偏离点都有决策记录 |

频率和几何指标**只报告、不设门槛**：哪种方法、哪个 Hessian 级别更好，是 MVP-1 要回答的研究问题（设计 §24 Contribution B），不是正确性判据。
