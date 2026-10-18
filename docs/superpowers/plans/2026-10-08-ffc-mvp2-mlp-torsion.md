# ffcompiler MVP-2 施工计划：MLP（MACE-POLAR-1-M）二面角接入

> **施工方（代理）（决策 0016，2026-10-08 用户决定由代理完全接手）。** 按任务顺序施工，每个任务先写测试再实现（TDD）。任务完成后写 commit 草稿（`.commit_drafts/`），不在 NFS 上跑 git。发现本计划自相矛盾或数学定义有问题，停下来报告，不自行改定义。

**目标：** 对每根可转键，用 MACE-POLAR-1-M 做约束松弛的 1D 二面角扫描（带总电荷，离子照常扫），扣掉 MM 的其余部分，拟合成 IR 的 AMBER 周期项 \(\sum_n k_n[1+\cos(n\phi-\delta_n)]\)，替换该键上原有的 GAFF2 二面角；失败一律退回 GAFF2，并逐键记录原因（设计 §10、§10.3 A、§23 MVP-2；决策 0016）。

**架构：** `torsion/` 子包分五层：可转键识别（`rotors.py`）→ 约束松弛扫描（`scan.py`，在 ASE 里做，对任意 ASE `Calculator` 通用：生产用 MACE-POLAR，测试用包了 OpenMM 的 MM 计算器）→ Fourier 拟合（`fourier.py`）→ MM 工具（`mm.py`：OpenMM→ASE 计算器、扣 MM）→ 驱动（`mlp.py`：拟合、写回 IR、兜底、provenance）。MACE 的加载单独放在 `mace.py`，只有它 import torch / mace。**不经 openmm-ml**（用户 2026-10-08 决定"直接 ase"）。

**技术栈：** Python 3.12、numpy、ASE 3.29.0（`FixInternals`、`LBFGS`）、mace-torch 0.3.16 + graph_longrange 0.4.0（`MACECalculator(model_type="PolarMACE")`）、OpenMM 8.5.2 Reference（MM 能量）、RDKit（SDF 键级、环信息）、pytest。全部已在 openmm_dev，不装任何包。

**依据：** `forcefield_compiler_design.md` §2、§4.1、§6.0、§10、§16、§20 Level 1、§23 MVP-2、§26；路线图 `2026-10-04-ffc-roadmap.md` §3、§3a（材料说明部分）；决策 0016。

**前置：** MVP-0、MVP-1 已完成（G0、G1 通过）。

---

## 0. 实测依据（2026-10-08，`experiments/20261008_mace_smoke/`）

| 事实 | 证据 |
|---|---|
| openmm-ml 1.6 不支持 PolarMACE：`MLPotential("mace", modelPath=MACE-POLAR-1-M)` 能加载，求能报 `'fermi_level'`（不传 `fermi_level` / `external_field`） | `polar_m_cpu.log` |
| mace 自带 ASE 计算器可用：`MACECalculator(model_paths=..., model_type="PolarMACE", default_dtype="float64", device="cpu")`，`atoms.info["charge"]`、`atoms.info["spin"]=1` | `polar_m_ase_cpu.log` |
| 总电荷生效：醋酸根 charge=−1 时原子电荷和 −1.0000（两个 O 约 −0.45）；同几何 charge=0 能量差 458 kJ/mol | 同上 |
| 模型：PolarMACE，权重 float32，覆盖 83 种元素，r_max 6 Å。`spin` 按官方示例（中性体系写 1）理解为自旋多重度——文档未明说，【未核实】 | 读模型；MACE 文档 polar_mace 页 |
| 乙醇一次能量+力约 340 ms（CPU 4 线程，float64） | `polar_m_ase_cpu.log` |
| `FixInternals` 二面角硬约束 + LBFGS（fmax 0.005 eV/Å），先用 `set_dihedral(mask=k 侧)` 刚性转 10° 再松弛：实际 φ 与目标一致到打印精度，每步 9–16 次迭代、2.4–4.2 s，键长最大变化 0.54 pm；醋酸根（−1）同样正常 | `polar_m_ase_scan.log` |
| 反例（MACE-OFF24 + OpenMM 软约束）：约束目标一次跳 60° 会把键拉坏（180° 点高 1446 kJ/mol），所以必须 10° 一步顺序扫 | `ethanol.log` |

## 1. 全局约束

MVP-0 计划 §1、MVP-1 计划 §1 的全部条目照常适用（openmm_dev、不装包、float64、Reference 平台、`-m external` 分层、日志、`RunContext` 运行目录、种子）。另外：

1. **模型**：`/home/user/MLP/mace/MACE-POLAR-1-M.model`，可用环境变量 `FFC_MACE_MODEL` 替换；不查 `~/.cache`。float64 求值，CPU（小分子 GPU 无收益，GPU 留给别人）。
2. **torch / mace 只在 `torsion/mace.py` 里 import**（函数内 import），其余模块与单元测试不得 import。ASE 可以在 `scan.py`、`mm.py` 里 import（轻量）。
3. 调用 MACE 的测试标新 marker `@pytest.mark.mlp`（在 `pytest.ini` 登记），同时也标 `external`。
4. **二面角形式只用 IR 已有的 AMBER 周期项**；不新增谐振二面角，Hessian 不推二面角（决策 0016 第 1 条，取代设计 §10.2 的"刚性二面角 ← Hessian"）。非可转键上的 GAFF2 二面角、全部 improper 原样保留。
5. **兜底**：任何失败都退回 GAFF2（保留该键原有项），在报告里写原因；不静默留空（设计 §6.0、决策 0016 第 2 条）。
6. **电荷与 1-4**：扣 MM 时用 IR 现有的（GAFF2 prmtop）电荷、LJ、1-4 缩放（决策 0016 第 3 条）。
7. **不做**：设计 §10.3 B（OpenMM `CustomTorsionForce` 原生模式）。A 模式三引擎原生兼容；B 模式等 A 的剖面复现不够时再议。
8. **线程**：每个扫描进程设 `OMP_NUM_THREADS=MKL_NUM_THREADS=4`，并 `torch.set_num_threads(4)`；Gate 并行时进程数 × 4 ≤ 机器物理核数减当前负载。
9. **长任务**：Gate 必须可续跑（§2.6）；启动前报计划表，等用户说"跑"。

---

## 2. 锁定的定义

### 2.1 可转键（`torsion/rotors.py`）

输入：IR 与同原子序的 RDKit 分子（来自固件 SDF，`removeHs=False`）。原子数或任一原子的原子序数与 IR 不一致 → `TopologyMismatchError`（复用 `bonded.assign` 里的类）。IR 的键集合与 RDKit 的键集合不一致也报这个错。

键 (j,k) 是**可转键** ⟺ 同时满足：
1. RDKit 键类型为 `SINGLE`（不是 aromatic、double、triple）；
2. 不在环内（`bond.IsInRing() == False`）；
3. j、k 的度数（含氢）都 ≥ 2；
4. j、k 都不是"线性中心"：在给定几何上，以 j（或 k）为顶点的任一键角 > 175°（与 MVP-1 §2.4 的线性阈值相同）。

可转键上的**四元组集合** = IR 中所有中心键为 {j,k} 的 proper torsion 的不同四元组（同一四元组的多个周期项算一个）。集合为空的可转键（GAFF2 没给项）也照样扫描，扫描四元组从图上生成。

**扫描四元组**（每根键一个）：在 (i, j, k, l) 中，i ∈ N(j)∖{k}、l ∈ N(k)∖{j}，按 (重原子数 i+l 降序, i 升序, l 升序) 取第一个。方向固定为 j < k（若 j > k，整体反序后再选）。

**等价类**：可转键按扫描四元组的图自同构轨道分类（`symmetry.term_classes` 扩展到四元组，(i,j,k,l) ≡ (l,k,j,i)），**只记录、不共享参数**：每根可转键各自扫描、各自拟合（T25 裁定）。原因：图自同构在几何上可能对应镜面（DMSO 两根 S–C 键），此时 \(E'(\phi)=E(-\phi)\)；拟合相位任意，直接复制 (k, δ) 会把剖面镜像错。同类键的参数差异记进报告。

### 2.2 约束松弛扫描（`torsion/scan.py`，ASE）

单位：ASE 内部 eV、Å；对外一律 kJ/mol、nm（1 eV = 96.485332123 kJ/mol，加进 `units.py` 为 `EV_KJ_MOL`）。

对任意 ASE `Calculator`（参考势）和扫描四元组 q = (i,j,k,l)：
1. **转动侧**：断开 j–k 后含 k 的连通分量（`mask`）。
2. **起点**：给定坐标无约束 LBFGS 松弛，得 \(\mathbf x_0\)、\(\phi_0\)。
3. **网格**：\(\phi_s=\phi_0+10°\cdot s\)，s = 0…35（36 点）。
4. **两遍顺序扫描**：正向从 \(\mathbf x_0\) 起 s = 1…35；反向从 \(\mathbf x_0\) 起 s = −1…−35（\(s\equiv s+36\)）。每步：清约束 → `set_dihedral(*q, φ_s, mask=mask)`（在上一步松弛几何上刚性转动）→ `FixInternals(dihedrals_deg=[[φ_s, q]])` → `LBFGS.run(fmax=0.005, steps=2000)`。
   - fmax 0.005 eV/Å = 0.48 kJ/mol/nm。依据：残余力 F 造成的能量误差量级 \(F^2/(2k)\)，最软模 k ≈ 100 kJ/mol/nm² 时 ≈ 1e-3 kJ/mol，远小于二面角能垒。
   - steps 上限 2000：冒烟实测每步 9–30 次迭代，上限是实测的 ~70 倍，只用来防死循环；到上限未收敛记为该点不收敛（§2.6）。
5. 每个网格点取两遍中**参考势能较低**的几何（硬约束下没有约束能），记录：目标 φ、实际 φ（IR 定义）、能量、坐标、两遍能量差（滞后量）、是否收敛。
6. 只报告、不设门槛的诊断量：两遍能量差的最大值、相对 \(\mathbf x_0\) 的最大键长变化（pm）。

### 2.3 扣 MM 与拟合目标（`torsion/mlp.py`）

对可转键 b，令 \(\text{IR}^{\setminus b}\) = IR 去掉所有中心键为 b 的 proper torsion（其余一切不变）。在每个网格几何 \(\mathbf x_s\) 上：
\[
y_s = E_\text{ref}(\mathbf x_s) - E_\text{MM}[\text{IR}^{\setminus b}](\mathbf x_s),
\]
\(E_\text{MM}\) 用 OpenMM Reference 求值（`to_openmm_system(..., validate_ir=False)`，在 `mm.py`）。拟合变量是扫描四元组的实际 φ（不是目标 φ）。

**已知近似（写进报告，不处理）**：同一键上其余四元组的能量全部并到扫描四元组；扫描时其他可转键上的二面角仍是 GAFF2，单遍、不迭代。

### 2.4 Fourier 拟合（`torsion/fourier.py`）

线性最小二乘（均匀权重）：
\[
y_s \approx c_0+\sum_{n=1}^{6}\big(a_n\cos n\phi_s+b_n\sin n\phi_s\big).
\]
换成 IR 形式：\(k_n=\sqrt{a_n^2+b_n^2}\)，\(\delta_n=\operatorname{atan2}(b_n,a_n)\)（\(a_n\cos n\phi+b_n\sin n\phi=k_n\cos(n\phi-\delta_n)\)，常数 \(\sum k_n\) 并入 \(c_0\)）。\(k_n \ge 0\) 自动成立；相位任意（不限 0/π）。n = 1…6 全部写入（拟合出的"零"系数是 1e-15 量级而不是 0，按大小删需要阈值，不设；T24 裁定）。

报告：拟合 RMSE、最大绝对残差（kJ/mol）、\(n_{\max}=6\)。

### 2.5 写回 IR 与 provenance

- 删除该可转键上的全部 proper torsion，在它的扫描四元组上写入 \(\{(n,\delta_n,k_n)\}\)，`source="mlp"`，`meta`：`model="MACE-POLAR-1-M"`、`scan_quartet`（字符串 "i-j-k-l"）、`fit_rmse`、`n_points=36`、`rotor_class`。
- IR `provenance["torsion"]`：模型路径、模型文件大小与修改时间、mace-torch / graph_longrange / ase 版本、总电荷与 spin、扫描参数（步长 10°、fmax、steps 上限）、每根可转键的结果或兜底原因。
- 写回后 `validate(ir)` 必须通过。

### 2.6 兜底（决策 0016 第 2 条）

| 原因代码 | 判定 | 范围 |
|---|---|---|
| `unsupported_element` | 分子含模型 `atomic_numbers` 以外的元素（运行时从模型读，83 种） | 整个分子 |
| `odd_electrons` | \(\sum Z-\text{total\_charge}\) 为奇数（本项目只做闭壳层，spin 固定 1） | 整个分子 |
| `mlp_error` | 扫描中参考势求值抛异常，或能量非有限 | 该可转键 |
| `scan_not_converged` | 任一网格点两遍都没在 2000 步内收敛 | 该可转键 |

兜底 = 保留该键原有的 GAFF2 项。带电分子**照常扫描**（总电荷传给模型，决策 0016 修订）。拟合残差、几何变形**不触发兜底**，只报告（全局规则：不自造门槛）。

### 2.7 续跑

Gate 把每根可转键的扫描结果存为 `<run>/scans/<molecule>/<j>-<k>.npz`（字段：目标 φ、实际 φ、参考能、坐标、滞后量、扫描参数）。重跑同一命令加 `--resume <run_dir>` 时，文件存在且其中记录的扫描参数与当前参数逐项相等就复用（不用 hash），否则重扫。

---

## 3. 文件结构

```text
src/ffcompiler/torsion/
├── __init__.py
├── rotors.py        # T25：可转键、扫描四元组、等价类
├── fourier.py       # T24：拟合与 IR 形式换算
├── scan.py          # T26：ASE 约束松弛扫描（对任意 ASE Calculator）
├── mm.py            # T26：OpenMMCalculator（IR→ASE 计算器）、去掉某键二面角的 MM 能量
├── mace.py          # T27：MACE-POLAR 计算器（唯一 import torch/mace 的模块）
└── mlp.py           # T28：MLPTorsionProvider、扣 MM、写回、兜底
src/ffcompiler/bonded/symmetry.py   # T25：term_classes 支持四元组
src/ffcompiler/cli.py               # T29：derive-torsions 子命令
scripts/gate_mvp2.py                # T30
tests/unit/test_torsion_*.py
```

---

## 4. Review Focus

1. **不该扫的键**：环内键（环丙烷、降冰片烷、苯环）、双键/芳香键、连着线性中心的键（乙腈 C–C≡N）。识别出的可转键必须为空。→ T25。
2. **带电分子**：醋酸根（−1）、甲铵（+1）。总电荷必须原样传给模型（`atoms.info["charge"]`），不能默认 0；奇电子数整分子兜底。→ T27、T28。
3. **元素不支持**：把一个原子的原子序数换成 Si 的合成 IR。整分子按 `unsupported_element` 兜底。→ T28。
4. **SDF 与 IR 原子序或键不一致**：必须报 `TopologyMismatchError`，不能拿错位的分子去扫。→ T25。
5. **对称等价键**：DMSO 两根 S–C 键同类（`rotor_class` 相同）、各自扫描；NMA 两个甲基各自成类。→ T25、T28。
6. **任意相位的跨引擎导出**：拟合出的 δ 不是 0/π，AMBER prmtop、GROMACS top 必须照样与 OpenMM 一致。→ T30（G2.3）。

---

## 5. 任务

### T24：Fourier 拟合

**Files:** Create `src/ffcompiler/torsion/__init__.py`、`src/ffcompiler/torsion/fourier.py`、`tests/unit/test_torsion_fourier.py`

**Interfaces（Produces）：**
- `@dataclass(frozen=True) FourierFit: k: tuple[float,...]（n=1..6）, phase: tuple[float,...], c0: float, rmse: float, max_abs_residual: float`
- `fit_fourier(phi_rad: np.ndarray, energy: np.ndarray, n_max: int = 6) -> FourierFit`
- `fourier_energy(fit: FourierFit, phi_rad) -> np.ndarray`（含 c0）
- `to_ir_terms(fit, quartet: tuple[int,int,int,int], *, source: str, meta: Mapping) -> list[ProperTorsion]`（n = 1..6 全写）

**必测用例：**
1. 合成数据：已知 \(k_n,\delta_n\)（含非 0/π 相位）+ 常数，在 36 个非均匀 φ 上取值 → 恢复 k、δ 的绝对误差 ≤ 1e-9（无噪声线性最小二乘，误差应在 1e-12 量级；1e-9 只用来抓公式错）。
2. `to_ir_terms` 生成的项按 `ir.energy` 的公式 \(\sum k_n[1+\cos(n\phi-\delta_n)]\) 求值，与 `fourier_energy` 在 100 个随机 φ 上之差为常数（差的标准差 ≤ 1e-9 kJ/mol）。
3. 少于 2·n_max+1 个点时报 `ValueError`。

### T25：可转键识别与四元组等价类

**Files:** Create `src/ffcompiler/torsion/rotors.py`、`tests/unit/test_torsion_rotors.py`；Modify `src/ffcompiler/bonded/symmetry.py`（`_normalize_term` 支持 4 元组，(i,j,k,l) 与 (l,k,j,i) 取字典序小者）、`tests/unit/test_symmetry.py`

**Interfaces（Produces）：**
- `@dataclass(frozen=True) Rotor: bond: tuple[int,int]（j<k）, scan_quartet: tuple[int,int,int,int], ir_quartets: tuple[tuple[int,int,int,int],...], rotor_class: int`
- `find_rotors(ir: ForceFieldIR, mol: Chem.Mol, positions_nm: np.ndarray) -> list[Rotor]`

**必测用例（12 个 GAFF2 固件，读 SDF）：**
1. 可转键数（按 §2.1 推出的预期；实现时先打印 SDF 键级核对，与预期不符就停下报告，不改规则凑数）：乙醇 2（C–C、C–O）、DMSO 2 且同类、NMA 3、苯甲酸甲酯 3、醋酸根 1、甲铵 1；苯、氯苯、环丙烷、环丁烷、降冰片烷、乙腈 0。
2. 原子序数换掉一个的 SDF / 少一根键的 IR → `TopologyMismatchError`。
3. `term_classes` 对四元组：DMSO 两根 S–C 键的扫描四元组同类；(i,j,k,l) 与 (l,k,j,i) 同类。

### T26：ASE 约束松弛扫描与 MM 工具

**Files:** Create `src/ffcompiler/torsion/scan.py`、`src/ffcompiler/torsion/mm.py`、`tests/unit/test_torsion_scan.py`；Modify `src/ffcompiler/units.py`（`EV_KJ_MOL = 96.485332123`）

**Interfaces（Produces）：**
- `mm.OpenMMCalculator(ase.calculators.calculator.Calculator)`：`__init__(ir, *, validate_ir=True)`，实现 `energy`、`forces`（OpenMM Reference，单位换成 eV、eV/Å）。
- `mm.mm_energy_without_bond(ir, bond: tuple[int,int], positions_nm) -> float`（kJ/mol）。
- `@dataclass ScanResult: quartet, phi_target_rad (36,), phi_rad (36,), energy (36,) kJ/mol, positions (36,N,3) nm, hysteresis (36,), converged (36,) bool, max_bond_change_pm: float, params: dict`
- `scan.relaxed_scan(calc, numbers, positions_nm, quartet, bonds, *, info=None, step_deg=10.0, fmax=0.005, max_steps=2000) -> ScanResult`；`info` 原样放进 `atoms.info`（电荷、spin）。
- 异常：参考势求值抛错或能量非有限 → `ScanError`（带四元组与网格点）。

**必测用例（参考势用 GAFF2 乙醇的 `OpenMMCalculator`，单元测试，不用 MACE）：**
1. 36 个实际 φ 与目标 φ 的差 ≤ 1e-6°（硬约束）。
2. 输出能量就是参考势本身：抽一个网格几何直接用 OpenMM 重算，差 ≤ 1e-8 kJ/mol。
3. 计算器返回 NaN → `ScanError`。
4. `OpenMMCalculator` 的力与 `ir.energy.forces_fd` 一致（换算单位后最大差 ≤ 1e-4 kJ/mol/nm；中心差分步长 1e-6 nm 的截断误差量级）。
5. `info={"charge": -1}` 原样出现在计算器收到的 `atoms.info` 里（用记录 info 的假计算器）。

### T27：MACE-POLAR 计算器

**Files:** Create `src/ffcompiler/torsion/mace.py`、`tests/unit/test_torsion_mace.py`（`@pytest.mark.mlp`）；Modify `pytest.ini`（登记 `mlp` marker）

**Interfaces（Produces）：**
- `default_model_path() -> Path`（`FFC_MACE_MODEL` 或 `/home/user/MLP/mace/MACE-POLAR-1-M.model`；不存在就 `FileNotFoundError`）
- `polar_calculator(*, model_path=None, device="cpu", threads=4) -> Calculator`（`MACECalculator(..., model_type="PolarMACE", default_dtype="float64")`，`torch.set_num_threads(threads)`）
- `supported_atomic_numbers(model_path=None) -> frozenset[int]`（从模型 `atomic_numbers` 读）
- `molecule_info(ir) -> dict`：`{"charge": ir.total_charge, "spin": 1}`

**必测用例（mlp）：**
1. 乙醇 fixture 坐标：能量 −407018.399 kJ/mol（冒烟值，差 ≤ 1e-2）。
2. 醋酸根 charge=−1：`calc.results["charges"]` 之和 −1（差 ≤ 1e-6）。
3. `supported_atomic_numbers` 含 1、6、7、8、9、15、16、17、35、53，共 83 个。

### T28：`MLPTorsionProvider`

**Files:** Create `src/ffcompiler/torsion/base.py`（`TorsionProvider` 抽象类，设计 §10：`parameterize(ir, positions_nm, mol) -> tuple[ForceFieldIR, TorsionReport]`）、`src/ffcompiler/torsion/mlp.py`、`tests/unit/test_torsion_mlp.py`

**Interfaces：**
- Consumes：T24 `fit_fourier`/`to_ir_terms`、T25 `find_rotors`、T26 `relaxed_scan`/`OpenMMCalculator`/`mm_energy_without_bond`、T27 `polar_calculator`/`supported_atomic_numbers`/`molecule_info`。
- Produces：
  - `@dataclass TorsionReport: rotors: list[dict]（bond、class、scan_quartet、status ∈ {"mlp","fallback"}、reason、fit_rmse、max_hysteresis、max_bond_change_pm）, molecule_fallback: str | None`
  - `class MLPTorsionProvider: def __init__(self, calculator_factory: Callable[[], Calculator] | None = None, *, supported: frozenset[int] | None = None, scan_kwargs=None, scan_cache: Path | None = None)`；`def parameterize(self, ir, positions_nm, mol) -> tuple[ForceFieldIR, TorsionReport]`
  - 默认 `calculator_factory=polar_calculator`、`supported=supported_atomic_numbers()`；测试里传 `lambda: OpenMMCalculator(ir)` 和显式元素集。

**必测用例：**
1. **自洽（核心）**：构造只有一根可转键、该键上只有一个四元组的玩具 IR（从乙醇 IR 删到 C–O 键只剩 H–O–C–C 一个四元组，原项 k、n、δ 已知），参考势 = 原 IR 自身的 `OpenMMCalculator` → 此时 \(y_s\) 恰好等于该四元组的二面角能（只依赖 φ，周期 ≤ 6，Fourier 可精确表示），拟合项与原项在 36 个扫描几何上之差为常数，标准差 ≤ 1e-8 kJ/mol（OpenMM Reference 双精度能量的舍入量级；与最小化收敛程度无关）。
2. 醋酸根（MM 参考势）：总电荷 −1 经 `molecule_info` 传入；奇电子数的合成 IR（删一个 H）→ `molecule_fallback="odd_electrons"`，输出 IR 的 proper torsions 与输入逐项相等。
3. 给 `supported` 传不含 Cl 的元素集，氯苯 → `unsupported_element`。
4. 参考势在第二根键上抛异常 → 该键 `fallback/mlp_error`，第一根键照常 `mlp`；`max_steps=1` → `scan_not_converged`。
5. DMSO（MM 参考势）：两根 S–C 键都是 `mlp`，各自的扫描四元组上写了 6 项；输出 IR 通过 `validate`。

### T29：CLI `derive-torsions`

**Files:** Modify `src/ffcompiler/cli.py`、`tests/unit/test_cli.py`

`ffcompiler derive-torsions --prmtop P --inpcrd C --sdf S --output {openmm,amber,gromacs}... --outdir D [--device cpu|cuda] [--ir-json F]`，写出各引擎文件、IR JSON 和 `torsion_report.json`。

**必测用例（mlp，external）：** 乙醇端到端跑通，`torsion_report.json` 里两根键都是 `mlp`（约 2 × 4 分钟，按 §0 实测每步 ~3.3 s × 70 步推算）。

### T30：Gate G2 脚本

**Files:** Create `scripts/gate_mvp2.py`、`tests/unit/test_gate_mvp2_report.py`

`python scripts/gate_mvp2.py [--jobs N] [--device cpu] [--resume RUN_DIR] [--skip-compare]`，`RunContext("gate_mvp2")`。分子：11 个 GAFF2 固件（乙腈按决策 0013 不进 Gate，与 G0/G1 相同；它没有可转键由 RF1 测试覆盖）。并行单位 = (分子, 可转键)，`ProcessPoolExecutor(max_workers=N)`，每进程 4 线程。

**每根可转键的指标（`metrics.csv`）：**
1. 拟合 RMSE（拟合项 vs \(y_s\)）。
2. 剖面复现：在扫描几何上，\(E_\text{MM}[\text{新 IR}]\) 与 \(E_\text{MLP}\) 各减去自身最小值后的 RMSE；同样算 \(E_\text{MM}[\text{GAFF2 IR}]\) vs \(E_\text{MLP}\)。
3. 滞后最大值、约束能最大值、最大键长变化（pm）。

另外：每个新 IR 跑 MVP-0 的 `compare`（OpenMM/AMBER/GROMACS），把扫描剖面存进运行目录供以后回归比对（设计 §23 MVP-2 "torsion-profile regression"）。

**必测用例：** 两行假指标生成报告，`report.md` 表头固定、`acceptance.json` 字段齐全；有一个组合抛异常时它进 `failures` 且 `passed=false`。

---

## 6. Gate G2 验收判据

| # | 判据 |
|---|---|
| G2.1 | T24 用例 1（Fourier 恢复）与 T28 用例 1（MM 自洽）通过 |
| G2.2 | 11 个固件每根可转键都有 `mlp` 结果或兜底原因；没有 NaN；可转键数与 T25 用例 1 一致 |
| G2.3 | 全部新 IR 通过 MVP-0 跨引擎比较（mixed GROMACS 按决策 0013 容差）。**2026-10-10 修订（决策 0017）**：以 OpenMM 为准，sander 腿只报告、不计失败 |
| G2.4 | Review Focus 1–6 对应的测试存在并通过（报告列出测试名） |
| G2.5 | pytest 全部通过；skip 只允许是 `gmx` 不存在或 CUDA 不可用 |

拟合 RMSE、剖面复现 RMSE、滞后、键长变化**只报告、不设门槛**：MLP 二面角比 GAFF2 好多少，是研究问题，不是正确性判据。

---

### G2 首次正式运行（2026-10-08，`jobs/gate_mvp2_20261008_3/`）

- 结果 **failed**：G2.1、G2.2、G2.4、G2.5 通过；G2.3 只有 nma 未过。12 根可转键全部 `mlp`，0 兜底，0 失败；pytest 非 external 410、external 28。
- G2.3 原因：唯一超差的腿是 sander-ours dihed，2.3e-4 > 容差 2e-4 kcal/mol。拟合项 0-1-3-4 n=2 相位 3.141059（离 π 5.3e-4 rad）。假设"sander 读 prmtop 时把 |δ−π| ≤ 4e-3 的相位改成 π"后，numpy 重算得 1.863826，打印为 1.8638，与 sander 相符；不做这个改动算出来是 1.863570，打印为 1.8636，与 sander 不符（`experiments/20261008_g2_nma_compare/snap_check.log`）。吸附阈值 4e-3 取自我对 sander 源码的记忆，**【未核实于源码】**。
- 只报告的观察：nma 1-3（酰胺 C–N）滞后 57.4 kJ/mol、键长变化 6.2 pm、剖面 RMSE 2.24（GAFF2 6.47）；methyl_benzoate 1-2 新 IR 剖面 RMSE 3.36，**比 GAFF2 的 2.70 还差**，拟合 RMSE 0.60。后者原因未查（候选：各键独立拟合时其余键仍是 GAFF2，最终 IR 里邻键已换成 MLP 项；也可能是"各减自身最小值"这个指标本身对单点敏感），**未核实**。

### G2 正式通过（2026-10-10，`jobs/gate_mvp2_20261010_2/`）

- 按决策 0017 修订后的 G2.3 重跑（`--resume jobs/gate_mvp2_20261008_3`，复用全部扫描；主机 node-B，762 s）：**passed**，G2.1–G2.5 全部通过。
- 12 根可转键全部 `mlp`，0 兜底，0 失败；11 个分子 sander 以外的跨引擎腿全部通过（GROMACS mixed）。pytest 非 external 413、external 28，0 失败 0 跳过。
- 已知局限（只报告）：nma sander-ours dihed 2.3e-4 > 2e-4 kcal/mol，与首跑相同。

## 7. 待用户确认

1. ~~带电分子整分子退回 GAFF2~~ 2026-10-08 改用 MACE-POLAR-1-M，离子照常扫（决策 0016 修订）。
2. 扫描几何变形（键长变化）是否要作为兜底条件：本计划只报告。若要设，需要先看 Gate 实测分布再定阈值。（2026-10-17 用户：先看分布，暂不设自动阈值。滞后分布见 `jobs/torsion_b3_sequential_20261017/report.md` 末节：nma 1-3 滞后 57.4 kJ/mol，单列为扫描质量问题。）
3. 同一键上全部能量并到一个四元组（§2.3 已知近似）：与 ABFE_IBS 的做法一致，未做多四元组分摊。
4. `spin=1` 按自旋多重度理解（官方示例），文档未明说；奇电子分子本期整分子兜底，所以只影响闭壳层是否该填 1。
5. ~~G2.3 的修法：拟合相位落在 sander 吸附窗口内时在 IR 里写成 0/π~~ 2026-10-10 用户决定不采用：IR 保留拟合相位，以 OpenMM 为准，G2 的 sander 腿只报告（决策 0017）。
6. methyl_benzoate 1-2 新 IR 比 GAFF2 差：是否要改成"其余键用已拟合的 MLP 项再拟合一轮"（迭代），或换剖面指标。**2026-10-10 用户：先记为已知缺陷，之后再定**（2026-10-17 B2 诊断：原因是单遍拟合的跨键耦合，见 `jobs/torsion_diagnosis_20261017/report.md`）（G2 正式运行中复现：剖面 RMSE 3.36 vs GAFF2 2.70，12 根键里唯一一根；已写进三语 README 已知限制）。
