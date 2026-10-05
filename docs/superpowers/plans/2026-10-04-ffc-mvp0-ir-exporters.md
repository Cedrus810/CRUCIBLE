# ffcompiler MVP-0 施工计划：ForceFieldIR + 三后端导出 + 跨引擎等价

> **给施工方：** 按任务顺序施工，每个任务先写测试再实现（TDD），任务完成即提交。本计划只给规格、接口、锁定的约定、必测用例和验收判据，不给实现代码。发现计划自相矛盾或约定有问题，**停下来报告**，不要自行改约定。

**目标：** 把一个已经完整参数化的分子（GAFF2 + AM1-BCC 的 prmtop）读进引擎无关的 `ForceFieldIR`，再写出 OpenMM / AMBER / GROMACS 三种格式，并证明三者的能量和力在数值容差内一致（设计 §19、§23 MVP-0）。

**架构：** `ForceFieldIR` 是不可变的纯数据结构，带锁定单位和函数形式（§2）。导入器和导出器只做格式转换，不做任何参数推导。另写一个 numpy 参考能量函数，直接从 IR 计算分项能量，作为所有引擎的共同参照，用来抓住各导出器共有的约定错误。

**技术栈：** Python 3.12、numpy、ParmEd（prmtop 读写）、OpenMM 8.5.2（Reference 平台、`AmberPrmtopFile`、`GromacsTopFile`）、AmberTools（antechamber、parmchk2、tleap、sander）、RDKit（固件构象）、pytest。

**依据：** `forcefield_compiler_design.md` §3、§4、§18、§19、§21、§23；路线图 `2026-10-04-ffc-roadmap.md`（决策 D1）。

---

## 1. 全局约束

每个任务都隐含满足本节。

1. **环境**：Python 只用 mamba 环境 `openmm_dev`：
   ```bash
   export MAMBA_EXE=/home/user/miniforge3/bin/mamba
   export MAMBA_ROOT_PREFIX=/home/user/miniforge3
   source /home/user/miniforge3/etc/profile.d/mamba.sh
   mamba activate openmm_dev
   ```
   **不得创建、安装、升级、降级任何环境或包，包括 dry-run。** 缺包就停下报告。这个环境里 `import torch` 很慢，MVP-0 的代码不得 import torch。
2. **IR 单位**：长度 nm、能量 kJ/mol、角度 rad、电荷 e、质量 amu。文件中以度为单位的字段名必须以 `_deg` 结尾。
3. **数值**：全部 float64。OpenMM 只用 `Reference` 平台。
4. **非键**：MVP-0 只做气相 NoCutoff，没有周期性边界。
5. **测试分层**：单元测试不调用外部程序。凡是调用 antechamber、tleap、sander、gmx 的测试都标 `@pytest.mark.external`；`gmx` 测试在 `shutil.which("gmx")` 为空时 skip（决策 D1）。`pytest` 默认跑全部；`pytest -m "not external"` 必须在 1 分钟内跑完。
6. **测试包遮蔽**：`openmm_dev` 的 site-packages 里有一个无关的 `tests` 包。`tests/` 和每个子目录都要有 `__init__.py`，`pytest.ini` 里设 `pythonpath = src` 和 `rootdir`。
7. **日志**：每个脚本的日志每行带 `[YYYY-mm-dd HH:MM:SS +elapsed]`，开头写 START 行（主机、argv、git commit），结尾写 END 行（退出码、总 wall time，用 atexit 保证写出）。pytest 统一经 `scripts/run_logged.sh` 运行并加 `--durations=0`。
8. **运行目录**：正式运行（Gate）输出到 `jobs/<name>_<YYYYMMDD>/`，包含 `manifest.json`（字段见 T0）、完整日志、`report.md`、`acceptance.json`。不放 `/tmp`。
9. **随机数**：显式种子，默认 `20261004`，写进 manifest。
10. **版本控制**：工作树 `/home/user/forcefield_compiler`，git 目录 `/home/user/.local/share/git-sep/ffcompiler.git`（见 T0）。

---

## 2. 锁定的约定

### 2.1 IR 函数形式

| 项 | 能量 | 参数（IR 单位） |
|---|---|---|
| bond | \(E=\tfrac12 k(r-r_0)^2\) | `r0` nm，`k` kJ/mol/nm² |
| angle | \(E=\tfrac12 k(\theta-\theta_0)^2\) | `theta0` rad，`force_constant` kJ/mol/rad² |
| proper torsion | \(E=k[1+\cos(n\phi-\delta)]\) | `periodicity` n ∈ {1..6}，`phase` δ rad，`force_constant` k kJ/mol（可以为负，也可以为 0） |
| improper | 与 proper 完全相同的周期形式 | 同上 |
| LJ | \(4\varepsilon_{ij}[(\sigma_{ij}/r)^{12}-(\sigma_{ij}/r)^6]\)，\(\sigma_{ij}=(\sigma_i+\sigma_j)/2\)，\(\varepsilon_{ij}=\sqrt{\varepsilon_i\varepsilon_j}\) | `sigma` nm，`epsilon` kJ/mol |
| Coulomb | \(k_e q_iq_j/r\)，\(k_e=138.935456\) kJ·nm·mol⁻¹·e⁻² | `partial_charge` e |
| 1-4 | `vdw_14_scale`·LJ + `coulomb_14_scale`·Coulomb | 全分子统一的两个缩放因子 |

- **二面角定义**（proper 和 improper 都用）：\(\mathbf b_1=\mathbf r_j-\mathbf r_i\)，\(\mathbf b_2=\mathbf r_k-\mathbf r_j\)，\(\mathbf b_3=\mathbf r_l-\mathbf r_k\)，\(\phi=\operatorname{atan2}\big(|\mathbf b_2|\,\mathbf b_1\cdot(\mathbf b_2\times\mathbf b_3),\ (\mathbf b_1\times\mathbf b_2)\cdot(\mathbf b_2\times\mathbf b_3)\big)\)。顺式为 0，与 OpenMM 一致。
- **improper 的原子顺序就是求值顺序**：IR 里存的 (i, j, k, l) 直接按上式算 φ。AMBER 约定中心原子在第三位；IR 不重排，导入时保持 prmtop 的顺序，导出时原样写出。
- **非键对的分类**（i<j）：
  - `exclusions`：不计任何非键作用的对。必须包含拓扑上全部 1-2、1-3 对。
  - `pairs14`：按缩放计入的 1-4 对，即"最短路径恰为 3 条键"的对。
  - 其余对按全强度计入。
  - 一个对不能同时在 `exclusions` 和 `pairs14` 里。MVP-0 要求 `pairs14` 恰好等于拓扑 1-4 集合（AMBER 家族约定），`exclusions` 恰好等于 1-2 ∪ 1-3 集合。
- **一个 (i,j,k,l) 可以有多个 torsion 项**，但同一个四元组不得出现两个相同 periodicity 的项。四元组 (i,j,k,l) 与 (l,k,j,i) 视为同一个。
- **电荷守恒**：\(|\sum_i q_i - Q|\le 10^{-4}\) e，Q 为 `total_charge`（整数）。

### 2.2 单位换算（导入 / 导出 AMBER 时用）

AMBER 的 bond/angle 不带 ½：\(E=K(r-r_0)^2\)。

| 量 | AMBER → IR |
|---|---|
| bond K（kcal/mol/Å²） | \(k = 2K \times 4.184 \times 100\) |
| bond r0（Å） | ×0.1 |
| angle K（kcal/mol/rad²） | \(k = 2K\times4.184\) |
| torsion K（kcal/mol，prmtop 已除过 IDIVF） | ×4.184 |
| LJ（ParmEd 的 `atom.rmin` 是 Rmin/2，Å；`atom.epsilon` kcal/mol） | \(\sigma = 2\,r_{\min,1/2}\cdot 2^{-1/6}\times0.1\)；\(\varepsilon\times4.184\) |
| SCEE / SCNB | `coulomb_14_scale = 1/SCEE`，`vdw_14_scale = 1/SCNB` |

### 2.3 OpenMM 力组（`decompose=True` 时）

| group | 内容 |
|---|---|
| 0 | bond |
| 1 | angle |
| 2 | proper torsion |
| 3 | improper |
| 4 | LJ（非 1-4） |
| 5 | Coulomb（非 1-4） |
| 6 | LJ 1-4 |
| 7 | Coulomb 1-4 |

`decompose=False`（生产模式）时：一个 `NonbondedForce` 承担全部非键，1-2/1-3 作为零强度 exception，1-4 作为缩放 exception；propers 和 impropers 各放一个 `PeriodicTorsionForce`。

### 2.4 跨引擎容差

能量按"类别"比较。各引擎能分出的类别不同，见 T9。

| 比较 | 能量（每个类别） | 力 |
|---|---|---|
| 内存内（OpenMM-IR、OpenMM-XML vs numpy-IR） | \(|\Delta E|\le 10^{-6}+10^{-9}|E|\) kJ/mol | max \(|\Delta F|\le 10^{-6}\) kJ/mol/nm |
| 经文件（OpenMM 读 prmtop、OpenMM 读 top vs numpy-IR） | \(|\Delta E|\le 10^{-4}+10^{-7}|E|\) kJ/mol | max \(|\Delta F|\le 10^{-3}\) kJ/mol/nm |
| sander(我们的 prmtop) vs sander(tleap 的 prmtop) | \(|\Delta E|\le 2\times10^{-4}+10^{-7}|E|\) kcal/mol | — |
| sander(我们的 prmtop) vs numpy-IR（换算到 kcal/mol） | 非静电类：\(2\times10^{-4}+10^{-7}|E|\)；静电类（EEL、1-4 EEL）：\(2\times10^{-4}+5\times10^{-5}|E|\) kcal/mol | — |

静电类放宽的原因：sander 的库仑常数是 \(18.2223^2=332.0522\)，OpenMM 是 \(138.935456/4.184\times10=332.0637\)，相对差约 \(3.5\times10^{-5}\)。

---

## 3. 文件结构

```text
forcefield_compiler/
├── pyproject.toml
├── pytest.ini
├── .gitignore
├── scripts/
│   ├── env.sh                    # 激活 openmm_dev
│   ├── run_logged.sh             # 带时间戳运行任意命令
│   ├── make_fixtures.py          # T4
│   └── gate_mvp0.py              # T10
├── src/ffcompiler/
│   ├── __init__.py
│   ├── runtime.py                # T0：RunContext、日志、manifest
│   ├── ir/
│   │   ├── schema.py             # T1
│   │   ├── io.py                 # T1：JSON
│   │   ├── topology.py           # T2：键图、1-2/1-3/1-4
│   │   ├── validation.py         # T2
│   │   └── energy.py             # T3：numpy 参考能量
│   ├── importers/
│   │   └── amber_prmtop.py       # T5
│   ├── exporters/
│   │   ├── openmm.py             # T6
│   │   ├── amber.py              # T7
│   │   └── gromacs.py            # T8
│   ├── validation/
│   │   ├── engines.py            # T9：各引擎求能量/力的适配器
│   │   └── cross_engine.py       # T9：比较与报告
│   └── cli.py                    # T10
├── tests/
│   ├── __init__.py
│   ├── conftest.py
│   ├── fixtures/gaff2/           # T4 产物（提交进仓库）
│   ├── unit/ ...
│   └── cross_engine/ ...
├── docs/decisions/               # 决策记录 0001-...
└── jobs/
```

---

## 4. Review Focus

规格隐含、但单个任务的正常测试很容易漏掉的五类输入，按出问题的可能性排序。每一条都已经落到对应任务的必测用例里。

1. **环上的 1-4 对**：4 元环、5 元环和稠环（降冰片烷）里，同一个 1-4 对可以经过两条二面角路径到达；三元环里没有 1-4 对。每个引擎都必须恰好算一次。→ T2、T7、T9。
2. **没有 torsion 项的 1-4 对**：AMBER prmtop 的 1-4 作用挂在二面角条目上。某个四元组没有任何 torsion 项（或所有 k=0 被省略）时，1-4 不能丢。→ T7。
3. **带电分子与电荷精度**：经 JSON、itp、prmtop 往返后，总电荷和每个原子的电荷都要保住。→ T1、T7、T8、T9（醋酸根、甲铵固件）。
4. **同一四元组的多项 torsion，以及 improper 的原子顺序**：GROMACS 的 funct 9 是多行，AMBER 的多项共用一个四元组，OpenMM 是多次 `addTorsion`；improper 不得被重排。→ T6、T7、T8。
5. **IR 不变量被破坏时必须报错而不是静默写出**：同一 atom_type 却有不同的 LJ 参数（GROMACS `[atomtypes]` 按类型名存）；prmtop 里 SCEE/SCNB 不统一。→ T2、T5。

---

## 5. 任务

### T0：仓库骨架与运行时工具

**Files:**
- Create: `pyproject.toml`、`pytest.ini`、`.gitignore`、`scripts/env.sh`、`scripts/run_logged.sh`、`src/ffcompiler/__init__.py`、`src/ffcompiler/runtime.py`、`tests/__init__.py`、`tests/conftest.py`、`tests/unit/__init__.py`、`tests/unit/test_runtime.py`、`docs/decisions/0001-git-layout.md`

**Interfaces:**
- Produces:
  - `ffcompiler.runtime.RunContext(name: str, root: Path = Path("jobs"), seed: int = 20261004, inputs: Sequence[Path] = ())`：上下文管理器。进入时创建 `jobs/<name>_<YYYYMMDD>/`（已存在就加 `_2`、`_3` 后缀，不覆盖），配置带时间戳的日志写到 `run.log`；退出时写 `manifest.json`。属性：`dir: Path`、`log: logging.Logger`、`manifest: dict`。
  - `manifest.json` 字段：`git_commit`、`git_dirty`、`argv`、`host`、`env_name`、`python`、`packages`（numpy、scipy、openmm、parmed、rdkit 版本）、`seed`、`inputs`（路径 → sha256）、`started`、`ended`、`wall_s`、`exit_status`。
  - `scripts/run_logged.sh <cmd...>`：行为与 `/home/user/Amber-X/scripts/run_logged.sh` 相同（每行时间戳、START/END、透传退出码）。

**git 布局**：
```bash
git init --separate-git-dir=/home/user/.local/share/git-sep/ffcompiler.git /home/user/forcefield_compiler
```
工作树根目录会出现一个 `.git` 文件指向本地盘。`jobs/` 只提交 `manifest.json`、`acceptance.json`、`report.md`，其余忽略。设计文档 `forcefield_compiler_design.md` 和 `.v0.md` 一并提交。

**必测用例**（`tests/unit/test_runtime.py`）：
1. `RunContext` 在临时目录下创建运行目录；同名同日期第二次运行得到 `_2` 后缀，第一次的内容不被覆盖。
2. 正常退出和抛异常退出时都写出 `manifest.json`，后者的 `exit_status` 非 0，`wall_s` ≥ 0。
3. `run.log` 每一行都匹配 `^\[\d{4}-\d\d-\d\d \d\d:\d\d:\d\d \+\d+s\] `。
4. `inputs` 里给的文件，sha256 与 `hashlib` 直接算的一致。
5. `scripts/run_logged.sh false` 的退出码为 1，输出最后一行含 `END exit=1`。

**完成标准**：`pytest -m "not external"` 通过；第一次提交。

---

### T1：IR schema 与 JSON 读写

**Files:**
- Create: `src/ffcompiler/ir/__init__.py`、`src/ffcompiler/ir/schema.py`、`src/ffcompiler/ir/io.py`、`tests/unit/test_ir_schema.py`

**Interfaces:**
- Produces（全部是 `@dataclass(frozen=True)`）：
  ```python
  SCHEMA_VERSION = "ffir/0.1"

  class Atom:
      index: int; name: str; element: str; atomic_number: int; mass: float
      atom_type: str; formal_charge: int; partial_charge: float
      sigma: float; epsilon: float; residue_name: str = "MOL"

  class Bond:
      i: int; j: int; r0: float; k: float
      source: str; meta: Mapping[str, str | int | float | bool] = {}; uncertainty: float | None = None

  class Angle:
      i: int; j: int; k_atom: int; theta0: float; force_constant: float
      source: str; meta: Mapping[...] = {}; uncertainty: float | None = None

  class ProperTorsion:
      i: int; j: int; k: int; l: int; periodicity: int; phase: float; force_constant: float
      source: str; meta: Mapping[...] = {}; uncertainty: float | None = None

  class ImproperTorsion:  # 字段与 ProperTorsion 完全相同
      ...

  class ForceFieldIR:
      name: str
      atoms: tuple[Atom, ...]
      bonds: tuple[Bond, ...]
      angles: tuple[Angle, ...]
      proper_torsions: tuple[ProperTorsion, ...]
      impropers: tuple[ImproperTorsion, ...]
      exclusions: tuple[tuple[int, int], ...]   # i<j，已排序
      pairs14: tuple[tuple[int, int], ...]      # i<j，已排序
      total_charge: int
      forcefield_family: str                    # MVP-0 只有 "amber"
      combination_rule: str                     # MVP-0 只有 "lorentz-berthelot"
      coulomb_14_scale: float
      vdw_14_scale: float
      provenance: Mapping[str, Any]             # 全局来源：输入文件、sha256、生成工具版本
  ```
  - `ForceFieldIR.replace_terms(*, bonds=None, angles=None, proper_torsions=None, impropers=None, provenance_update=None) -> ForceFieldIR`：返回新 IR，原 IR 不变（MVP-1 用）。
  - `ffcompiler.ir.io.to_json(ir, path: Path) -> None`、`from_json(path: Path) -> ForceFieldIR`。JSON 顶层带 `"schema": "ffir/0.1"`；版本不符时报 `ValueError`。浮点数必须精确往返（用 `repr`/17 位有效数字，不得格式化截断）。
- 约定：IR 不含坐标。坐标单独以 `np.ndarray (N,3)`（nm）传递。

**必测用例**：
1. 构造一个 4 原子 IR（含 1 个 bond、1 个 angle、同一四元组 n=1,2,3 三项 proper、1 个 improper），JSON 往返后 `==` 原对象，且每个浮点数逐位相等。
2. 带电 IR：`partial_charge` 为 `[-0.8123456789012345, 0.3123456789012345, ...]`，往返后逐位相等（对应 Review Focus 3）。
3. 读到 `"schema": "ffir/0.2"` 时报 `ValueError`，信息里含两个版本号。
4. `replace_terms` 返回新对象；原对象的 bonds 不变；`provenance_update` 合并进新对象的 provenance，原对象不变。
5. dataclass 是 frozen：给 `atom.partial_charge` 赋值抛 `FrozenInstanceError`。

---

### T2：拓扑工具与 IR 校验

**Files:**
- Create: `src/ffcompiler/ir/topology.py`、`src/ffcompiler/ir/validation.py`、`tests/unit/test_topology.py`、`tests/unit/test_validation.py`

**Interfaces:**
- Consumes: T1 的 `ForceFieldIR`。
- Produces:
  - `bond_graph(ir) -> dict[int, frozenset[int]]`：邻接表，只由 `ir.bonds` 决定。
  - `topological_pairs(ir) -> tuple[frozenset[tuple[int,int]], frozenset[tuple[int,int]], frozenset[tuple[int,int]]]`：按**最短路径键数**分出 (1-2, 1-3, 1-4) 三个集合，元素 i<j，三者互不相交。
  - `expected_exclusions_and_pairs(ir) -> tuple[tuple[tuple[int,int],...], tuple[tuple[int,int],...]]`：按 §2.1 的约定生成 `exclusions`、`pairs14`（已排序）。导入器和测试都用它。
  - `class IRValidationError(ValueError)`，属性 `issues: list[str]`。
  - `validate(ir) -> None`：发现任何问题就一次性抛出全部 issues，不能只报第一个。

**校验项**（每一项都要有对应的反例测试）：
1. 原子 `index` 恰为 0..N-1 且有序；所有项的下标在范围内。
2. 每个 angle (i,j,k) 的 i-j 和 j-k 都是 bond；每个 proper (i,j,k,l) 的 i-j、j-k、k-l 都是 bond；每个 improper (i,j,k,l) 的第三个原子 k 与 i、j、l 都成键。
3. 没有重复的 bond / angle；同一四元组（含反向）没有重复的 periodicity。
4. bond `k>0`、`r0>0`；angle `force_constant>0`、`0<theta0≤π`；periodicity ∈ {1..6}。
5. `exclusions`、`pairs14` 与 `expected_exclusions_and_pairs(ir)` 完全相等。
6. \(|\sum q - Q|\le10^{-4}\)。
7. 同一个 `atom_type` 的 `sigma`、`epsilon`、`mass`、`atomic_number` 必须完全一致（Review Focus 5）。
8. `forcefield_family == "amber"`，`combination_rule == "lorentz-berthelot"`；两个 1-4 缩放都在 (0, 1]。

**必测用例**（`test_topology.py`，全部用手工构造的 IR，只需要原子和 bond）：
1. 正丁烷骨架 C0-C1-C2-C3：1-4 集合 = {(0,3)}。
2. 环丙烷 3 个 C：全部是 1-2，1-3 和 1-4 都为空。
3. 环丁烷 4 个 C：对角 (0,2)、(1,3) 是 1-3，没有 1-4。
4. 环戊烷 5 个 C：没有 1-4（环上任意两点最短路径 ≤2）。带上 H 后，H 与隔两根键的 C 是 1-4，且每对只出现一次。
5. 降冰片烷骨架（7 个 C，C1 和 C4 为桥头）：每个 1-4 对只出现一次，且没有对同时落在 1-3 和 1-4（Review Focus 1）。

**必测用例**（`test_validation.py`）：上面 8 个校验项各一个反例，断言 `issues` 里含对应信息；再加一个同时违反 3 项的 IR，断言 3 条 issue 全部报出。

---

### T3：numpy 参考能量

**Files:**
- Create: `src/ffcompiler/ir/energy.py`、`tests/unit/test_energy.py`

**Interfaces:**
- Produces:
  - `COULOMB_KJ_NM = 138.935456`
  - `ENERGY_KEYS = ("bond", "angle", "proper", "improper", "lj", "coulomb", "lj14", "coulomb14")`
  - `evaluate(ir, positions_nm: np.ndarray) -> dict[str, float]`：返回 `ENERGY_KEYS` 各项和 `"total"`，单位 kJ/mol，严格按 §2.1 实现。
  - `dihedral_angle(p0, p1, p2, p3) -> float`：§2.1 的定义，rad，范围 (-π, π]。

**必测用例**：
1. 两原子 bond：r=0.11 nm、r0=0.1、k=1000，E_bond = 0.05 kJ/mol（精确到 1e-12）。
2. 三原子 angle：θ=100°、θ0=90°、k=500，E = ½·500·(π/18)²。
3. `dihedral_angle`：顺式平面构型 = 0；反式 = π；绕 b2 右手旋转 +60° 的构型 = +π/3（与 OpenMM 符号一致：同一构型用 OpenMM `PeriodicTorsionForce(n=1, δ=0, k=1)` 求能量，E = 1+cos φ 一致）。
4. LJ 与 Coulomb：两原子相距 r，不在任何排除表里，结果与手算公式一致；同一对放进 `pairs14` 后只出现在 `lj14`/`coulomb14`，数值乘缩放因子；放进 `exclusions` 后为 0。
5. 一个 4 原子分子在 10 个随机构型（种子 20261004）上，`evaluate` 的总能量与 T6 之前用 OpenMM 手工搭的 System（同参数）一致到 1e-9 相对误差。这条测试在 T6 完成后改用 `to_openmm_system`，两者都保留。

---

### T4：GAFF2 固件

**Files:**
- Create: `scripts/make_fixtures.py`、`tests/fixtures/gaff2/<name>/{<name>.sdf, <name>.mol2, <name>.frcmod, <name>.prmtop, <name>.inpcrd}`、`tests/fixtures/gaff2/manifest.json`、`tests/cross_engine/test_fixtures_present.py`

**规格**：
- 分子清单（锁定，名字即目录名）：

  | name | SMILES | 覆盖什么 |
  |---|---|---|
  | ethanol | `CCO` | 基础 |
  | benzene | `c1ccccc1` | 芳香、improper |
  | nma | `CC(=O)NC` | 酰胺 improper |
  | cyclopropane | `C1CC1` | 三元环 |
  | cyclobutane | `C1CCC1` | 四元环 1-3/1-4 |
  | norbornane | `C1CC2CCC1C2` | 稠环、1-4 多路径 |
  | acetate | `CC(=O)[O-]` | 负电荷 |
  | methylammonium | `C[NH3+]` | 正电荷 |
  | chlorobenzene | `Clc1ccccc1` | 卤素 |
  | methyl_benzoate | `COC(=O)c1ccccc1` | 与 Sobtop 样例同一个分子 |
  | dmso | `CS(C)=O` | 硫 |
  | acetonitrile | `CC#N` | 线性键角（MVP-1 用） |

- 流程：RDKit 加氢 → ETKDGv3（种子 20261004）→ MMFF94 优化 → 写 SDF → `antechamber -fi sdf -fo mol2 -c bcc -at gaff2 -nc <Q>` → `parmchk2 -s gaff2` → tleap（`source leaprc.gaff2`，`loadamberparams`，`saveamberparm`）。
- `manifest.json` 记录：每个分子的 SMILES、总电荷、AmberTools 版本、各文件 sha256、生成时间。
- 固件**提交进仓库**；测试只读这些文件，不现场重跑 antechamber。重新生成要显式运行脚本。
- 任何一个分子 antechamber 或 tleap 失败：脚本以非 0 退出，并报告是哪个分子、哪一步。不得跳过后继续。

**必测用例**（`external`）：`make_fixtures.py --only ethanol --out <tmp>` 能从头生成；prmtop 原子数与 SDF（含 H）一致；每个分子 prmtop 的电荷总和与 `-nc` 一致到 1e-4。

**必测用例**（非 external）：12 个目录和 manifest 都存在，sha256 与 manifest 一致。

---

### T5：prmtop 导入器

**Files:**
- Create: `src/ffcompiler/importers/__init__.py`、`src/ffcompiler/importers/amber_prmtop.py`、`tests/unit/test_import_prmtop.py`

**Interfaces:**
- Consumes: T1、T2；T4 固件。
- Produces:
  - `from_amber(prmtop: Path, inpcrd: Path, sdf: Path, name: str | None = None) -> tuple[ForceFieldIR, np.ndarray]`：返回 IR 和坐标（nm，取自 inpcrd）。
    - 参数按 §2.2 换算；每一项的 `source="gaff2-library"`，`meta` 里记 AMBER 原始类型名。
    - `formal_charge` 来自 SDF（按原子顺序对应；元素序列不一致就报错）。
    - impropers 用 ParmEd 的 improper 标记识别，原子顺序保持 prmtop 原样。
    - `exclusions`、`pairs14` 由 `expected_exclusions_and_pairs` 生成，然后**与 prmtop 自身的排除表核对**：prmtop 排除表（1-2/1-3/1-4 合集）必须等于 IR 的 exclusions ∪ pairs14，不一致就报错。
    - 所有二面角条目的 SCEE、SCNB 必须统一，否则报 `ValueError`（Review Focus 5）。
    - 结尾调用 `validate(ir)`。

**必测用例**：
1. 12 个固件全部能导入并通过 `validate`。
2. ethanol：原子数、bond/angle/proper 数与 ParmEd 直接数的一致；某条 C-C bond 的 k 满足 \(k_{IR}=2\times K_{prmtop}\times418.4\)（精确到 1e-12 相对误差）。
3. 苯：improper 数 = 6（GAFF2 对每个芳香 C 一个），每个 improper 的第三个原子是芳香 C。
4. 醋酸根、甲铵：`total_charge` 分别为 -1、+1；`formal_charge` 非零的原子是 O⁻ 和 N⁺。
5. 篡改过的 prmtop（测试里用 ParmEd 把某个二面角的 SCEE 改为 1.0 再存到临时文件）→ 导入报错，信息含 "SCEE"。
6. SDF 元素顺序与 prmtop 不一致（测试里交换两个原子）→ 报错。
7. 对每个固件，`evaluate(ir, pos)` 与 OpenMM `AmberPrmtopFile(原 prmtop).createSystem(nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)` 的总能量差满足 §2.4 "经文件"容差。这条证明导入器和 numpy 参考都对。

---

### T6：OpenMM 导出器

**Files:**
- Create: `src/ffcompiler/exporters/__init__.py`、`src/ffcompiler/exporters/openmm.py`、`tests/unit/test_export_openmm.py`

**Interfaces:**
- Consumes: T1–T3。
- Produces:
  - `FORCE_GROUPS: dict[str, int]`：§2.3 的映射，键与 `ENERGY_KEYS` 相同。
  - `to_openmm_system(ir, *, decompose: bool = False) -> openmm.System`：先 `validate(ir)`。
  - `to_openmm_topology(ir) -> openmm.app.Topology`：一个残基，名字取 `atoms[0].residue_name`；元素来自 `atomic_number`；bond 来自 `ir.bonds`。
  - `write_openmm(ir, positions_nm, outdir: Path) -> dict[str, Path]`：写 `system.xml`（`XmlSerializer`，生产模式）、`molecule.pdb`；返回 `{"system": ..., "pdb": ...}`。

**必测用例**：
1. 对 T3 的手工 4 原子 IR 和全部 12 个固件，`decompose=True` 时各力组能量与 `evaluate` 的对应项满足"内存内"容差。
2. 同一 IR，`decompose=False` 的总能量与 `decompose=True` 的总能量一致到 1e-9 相对误差。
3. 同一四元组 n=1,2,3 三项 proper：OpenMM 里是 3 次 `addTorsion`，能量与 `evaluate` 一致（Review Focus 4）。
4. improper 的原子顺序：构造一个 improper 使 (i,j,k,l) 与 (l,j,k,i) 的能量不同；导出后 OpenMM 能量等于按 IR 原顺序算的值（Review Focus 4）。
5. `system.xml` 反序列化后，在 3 个构型上与内存中的 System 能量和力一致（"内存内"容差）。
6. 力：每个固件、每个构型，OpenMM 的力与 `evaluate` 中心差分力（步长 1e-6 nm）的 max 差 ≤ 1e-4 kJ/mol/nm。

---

### T7：AMBER 导出器

**Files:**
- Create: `src/ffcompiler/exporters/amber.py`、`tests/unit/test_export_amber.py`

**Interfaces:**
- Consumes: T1–T3。
- Produces:
  - `write_amber(ir, positions_nm, outdir: Path) -> dict[str, Path]`：经 ParmEd `Structure` 写 `molecule.prmtop` 和 `molecule.inpcrd`；返回 `{"prmtop": ..., "inpcrd": ...}`。先 `validate(ir)`。

**规格**：
- 参数按 §2.2 逆换算。参数按实例给出：ParmEd 会为每组不同取值建一个类型条目，这是允许的。
- improper 条目标为 improper，原子顺序保持 IR 原样。
- **1-4 规则**：`pairs14` 中每个对必须恰好挂在一个"计算 1-4"的二面角条目上；同一对出现在多个二面角上时，除一个之外全部标为不计 1-4（ParmEd 的 `ignore_end=True`）。如果某个 1-4 对没有任何 proper torsion 项，就补一条 k=0、n=1、δ=0 的二面角条目来承载它（Review Focus 2）。
- SCEE = 1/`coulomb_14_scale`，SCNB = 1/`vdw_14_scale`，写到每个二面角条目。
- 不写 frcmod（路线图 §3.4）。

**必测用例**：
1. 每个固件：`from_amber` → `write_amber` → OpenMM `AmberPrmtopFile` 读回，"经文件"容差内与 `evaluate` 一致（类别：bond、angle、torsion_all、nonbonded_all、total）。
2. **没有 torsion 的 1-4 对**：手工 IR，正丁烷骨架加 H，删掉 H-C-C-H 中某一个四元组的全部 torsion 项。写出的 prmtop 里该 1-4 对仍然存在：OpenMM 读回的 nonbonded 能量与 `evaluate` 一致（Review Focus 2）。
3. **多路径 1-4**：降冰片烷和环丁烷固件，写出的 prmtop 中每个 1-4 对恰好计算一次：OpenMM 读回的 nonbonded 能量与 `evaluate` 一致（Review Focus 1）。
4. 带电固件：读回的电荷总和与 `total_charge` 差 ≤ 1e-6（Review Focus 3）。
5. 同一四元组多项 proper 写出后 term 数与 IR 一致（ParmEd 读回计数）。

---

### T8：GROMACS 导出器

**Files:**
- Create: `src/ffcompiler/exporters/gromacs.py`、`tests/unit/test_export_gromacs.py`

**Interfaces:**
- Consumes: T1–T3。
- Produces:
  - `write_gromacs(ir, positions_nm, outdir: Path) -> dict[str, Path]`：写 `molecule.itp`、`topol.top`、`conf.gro`；返回 `{"itp": ..., "top": ..., "gro": ...}`。先 `validate(ir)`。

**格式规格**：
- `topol.top`：`[ defaults ]` = `1 2 yes <vdw_14_scale> <coulomb_14_scale>`；`#include "molecule.itp"`；`[ system ]`；`[ molecules ]`。不 include 任何 GROMACS 自带力场文件（自包含）。
- `molecule.itp`：
  - `[ atomtypes ]`：每个 atom_type 一行，`name at.num mass 0.0 A sigma epsilon`。
  - `[ moleculetype ]`：`nrexcl = 3`。
  - `[ atoms ]`：电荷用 `%.10f`；残基号 1。
  - `[ bonds ]` funct 1：`b0 kb`。
  - `[ pairs ]` funct 1：`pairs14` 全部列出。
  - `[ angles ]` funct 1：`theta0_deg k`。
  - `[ dihedrals ]` proper 用 funct 9，同一四元组每项一行：`phase_deg k n`。
  - `[ dihedrals ]` improper 用 funct 4，一个四元组只允许一项；多于一项就报错（GAFF2 的 improper 都是单项）。
  - `[ exclusions ]`：只写 `nrexcl=3` 覆盖不到的额外排除（MVP-0 下应为空，写出空段落也可以不写）。
  - 所有浮点数至少 10 位有效数字（`%.10g` 或更多）。
- `conf.gro`：标准固定列宽格式，坐标 3 位小数（只用于给人看；测试用 IR 坐标，不读 gro）。

**必测用例**：
1. 每个固件：`write_gromacs` → OpenMM `GromacsTopFile(top)` 读回（NoCutoff），用 IR 坐标求能量，"经文件"容差内与 `evaluate` 一致（类别同 T7）。
2. 同一四元组 n=1,2,3 三项 proper：itp 里有 3 行 funct 9，能量一致（Review Focus 4）。
3. 两项 improper 落在同一四元组 → 报错，信息含四元组下标。
4. 带电固件：itp 中 `[ atoms ]` 的电荷总和与 `total_charge` 差 ≤ 1e-8（Review Focus 3）。
5. 写出的 `[ pairs ]` 行数等于 `len(ir.pairs14)`；降冰片烷读回能量一致（Review Focus 1）。
6. `external`（`gmx` 不存在时 skip）：`gmx grompp` + `gmx mdrun -rerun` 单点，能量与 `evaluate` 一致（容差同"经文件"，静电类按 \(5\times10^{-5}\) 相对放宽）。

---

### T9：引擎适配器与跨引擎比较

**Files:**
- Create: `src/ffcompiler/validation/__init__.py`、`src/ffcompiler/validation/engines.py`、`src/ffcompiler/validation/cross_engine.py`、`tests/cross_engine/__init__.py`、`tests/cross_engine/test_engines.py`、`tests/cross_engine/test_cross_engine.py`

**Interfaces:**
- Consumes: T3、T5–T8。
- Produces:
  - `engines.openmm_energy_forces(system, positions_nm, *, groups: Mapping[str, int] | None) -> tuple[dict[str, float], np.ndarray]`：按力组或按力的类别（`HarmonicBondForce`→bond、`HarmonicAngleForce`→angle、`PeriodicTorsionForce`→torsion_all、`NonbondedForce`/`CustomBondForce`/`CustomNonbondedForce`→nonbonded_all）汇总能量，kJ/mol；力 (N,3) kJ/mol/nm。
  - `engines.sander_energies(prmtop: Path, inpcrd: Path, workdir: Path) -> dict[str, float]`：单点能量，键为 `bond, angle, dihed, vdw, eel, vdw14, eel14`，单位 kcal/mol。mdin 规格：`ntb=0, igb=0, cut=999.0`，气相，不做任何优化步；必须是**输入坐标处**的能量。解析 mdout 的第一个能量块；缺任何一个键就报错。
  - `cross_engine.compare(ir, positions_nm, workdir: Path, *, tleap_prmtop: Path | None, tleap_inpcrd: Path | None) -> CrossEngineReport`：跑下面的各条腿，套 §2.4 容差。
  - `CrossEngineReport`：dataclass，`rows: list[dict]`（每行：leg、reference、category、E_leg、E_ref、delta、tol、passed）、`force_rows`、`passed: bool`；`to_markdown() -> str`、`to_json() -> dict`。

**各条腿**：

| 腿 | 内容 | 能比的类别 | 参照 |
|---|---|---|---|
| `numpy-ir` | `evaluate` | 8 项 + total | 参照本身 |
| `openmm-ir` | `to_openmm_system(decompose=True)` | 8 项 + total | numpy-ir（内存内） |
| `openmm-xml` | 反序列化 `system.xml` | bond、angle、torsion_all、nonbonded_all、total | numpy-ir（内存内） |
| `openmm-prmtop` | `AmberPrmtopFile(我们的 prmtop)` | 同上 | numpy-ir（经文件） |
| `openmm-top` | `GromacsTopFile(我们的 top)` | 同上 | numpy-ir（经文件） |
| `sander-ours` | sander（我们的 prmtop） | bond、angle、dihed(=proper+improper)、vdw、eel、vdw14、eel14 | numpy-ir（sander 容差）；sander-tleap |
| `sander-tleap` | sander（tleap 原 prmtop） | 同上 | — |
| `gmx` | 真 GROMACS（没有就 skip，报告里标 SKIPPED） | total | numpy-ir |

力只比 OpenMM 系的四条腿，参照为 `openmm-ir`。

**必测用例**：
1. `sander_energies` 自检（`external`）：对 ethanol 的 tleap prmtop，sander 结果与 OpenMM `AmberPrmtopFile` 同一文件的能量满足 sander 容差。这证明 mdin 和解析是对的。
2. `sander_energies` 遇到截断或损坏的 mdout（测试里提供一个手工截断的样例文件）→ 报错，不返回部分结果。
3. 12 个固件 × 4 个构型（原始 + 3 个扰动：每个坐标加 N(0, 0.005 nm)，种子 20261004）全部 `passed`（`external`）。
4. 故意坏掉的 IR（把某个 bond 的 k 乘 1.001 再只写 prmtop）→ `openmm-prmtop` 与 `openmm-ir` 不一致，报告 `passed=False`，且失败行指向 bond 类别。用来证明比较器真的会报错。
5. `to_markdown()` 输出的表头固定、每行一个 (腿, 类别)；`to_json()` 能被 `json.dumps` 序列化。

---

### T10：CLI 与 Gate G0

**Files:**
- Create: `src/ffcompiler/cli.py`、`scripts/gate_mvp0.py`、`tests/unit/test_cli.py`；`pyproject.toml` 加 `[project.scripts] ffcompiler = "ffcompiler.cli:main"`

**Interfaces:**
- Produces:
  - `ffcompiler convert --prmtop P --inpcrd C --sdf S --output {gromacs,amber,openmm}... --outdir D [--ir-json F]`：导入后按 `--output` 写出；`--ir-json` 同时存 IR。
  - `ffcompiler check-equivalence --prmtop P --inpcrd C --sdf S --workdir D`：跑 `compare`，打印 markdown 表，退出码 0/1 表示通过/失败。
  - `scripts/gate_mvp0.py`：在 `RunContext("gate_mvp0")` 下对全部固件跑 `compare`，写 `report.md`（每个固件一节）和 `acceptance.json`。

**`acceptance.json` 字段**：`gate: "G0"`、`passed: bool`、`n_fixtures`、`n_geometries`、`legs`（每条腿：passed / failed / skipped 计数）、`max_rel_error_by_category`、`skipped_legs`（含原因，如 `gmx: not installed`）、`pytest_summary`（通过/失败/跳过数）。

**必测用例**：
1. `convert` 对 ethanol 输出三种格式的全部文件；`--output` 只给 `amber` 时不生成 GROMACS 文件。
2. `check-equivalence` 对 ethanol 退出码 0；对 T9 的坏 IR 场景退出码 1。
3. 参数缺失（如没给 `--sdf`）时退出码 2，打印用法。

---

## 6. Gate G0 验收判据

施工方提交以下材料，代理只审阅、不重跑：

1. `jobs/gate_mvp0_<date>/`：`manifest.json`、`run.log`（带时间戳、START/END）、`report.md`、`acceptance.json`。
2. 全套 pytest 日志（经 `run_logged.sh`，含 `--durations=0`）：非 external 和 external 都要跑。

通过条件（全部满足）：

| # | 判据 |
|---|---|
| G0.1 | `acceptance.json` 的 `passed=true`；12 个固件 × 4 个构型，除 `gmx` 外的全部腿通过 §2.4 容差 |
| G0.2 | `gmx` 腿标 SKIPPED 且注明原因（决策 D1）；如果有 GROMACS 可用，则必须通过 |
| G0.3 | pytest 全部通过，0 失败；skip 只允许是 `gmx` 不存在 |
| G0.4 | Review Focus 1–5 各自对应的测试存在并通过（报告中列出测试名） |
| G0.5 | T9 用例 4（故意坏掉的 IR）确实失败并被正确定位 |
| G0.6 | 没有修改 §2 的任何约定；如有偏离，`docs/decisions/` 里有记录并已获用户确认 |
