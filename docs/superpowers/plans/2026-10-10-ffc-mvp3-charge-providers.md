# ffcompiler MVP-3 施工计划：电荷 provider 框架（载入 / AM1-BCC / NAGL）

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
>
> **施工方（代理）（决策 0016）。** 按任务顺序施工，每个任务先写测试再实现（TDD）。按 roadmap 约定，计划只给规格、接口、锁定的定义、必测用例和验收判据，不给实现代码。任务完成后写 commit 草稿（`.commit_drafts/`），不在 NFS 上跑 git。发现本计划自相矛盾或定义有问题，停下来报告，不自行改定义。

**目标：** 把电荷做成可替换的 provider（设计 §11），本期实现三种：载入、AM1-BCC（antechamber/sqm）、NAGL（openff-gnn-am1bcc-1.0.0）。电荷写回 IR 时总电荷严格守恒；给已有 MLP 二面角的 IR 换电荷时，用缓存的扫描自动重拟合二面角（决策 0018）。

**架构：** 新子包 `charge/`：provider 只算原始电荷（`compute`），`apply.assign_charges` 负责均摊残差、写回 IR、provenance、触发二面角重拟合。二面角侧在 `torsion/mlp.py` 记录被替换的原 GAFF2 项，并新增 `refit_torsions`。CLI 新增 `assign-charges`，`derive-torsions` 增加 `--scan-cache`、`--ir-in`。Gate G3 在 `scripts/gate_mvp3.py`。

**技术栈：** Python 3.12、numpy、RDKit、AmberTools（antechamber、sqm，openmm_dev/bin）、openff-nagl 0.6.1 + openff-nagl-models 2025.9.0（`openff-gnn-am1bcc-1.0.0.pt`）、openff-toolkit 0.19.0、OpenMM 8.5.2 Reference、pytest。全部已在 openmm_dev，不装任何包。

**依据：** `forcefield_compiler_design.md` §4.1、§6.4.2、§11、§15、§16、§23 MVP-3；roadmap §1、§3 第 3 条；决策 0016、0017、0018；冒烟 `experiments/20261010_am1bcc_repro/`。

**前置：** MVP-0、MVP-1、MVP-2 已完成（G0、G1、G2 通过；G2 = `jobs/gate_mvp2_20261010_2/`）。

---

## 0. 实测依据（2026-10-10，`experiments/20261010_am1bcc_repro/`，node-A 单线程）

| 事实 | 证据 |
|---|---|
| 按 `scripts/make_fixtures.py` 的做法（`antechamber -c bcc -at gaff2 -nc Q -pf y`，再用同一个残差均摊规则）在固件 SDF 上重跑：12/12 固件与固件 mol2 电荷逐位相同；每个分子 0.5–2.4 s | `repro.log`、`results.json` |
| 均摊前（antechamber 写的 4 位小数）与固件差 ≤ 3.3e-4 e，就是均摊量本身 | 同上 |
| NAGL `openff-gnn-am1bcc-1.0.0`：12/12 在适用域内；与 AM1-BCC 逐原子差 4.3e-5–0.040 e；总电荷误差 ≤ 1.8e-7（float32）；每个分子 0.003–0.74 s | `nagl.log`、`nagl_results.json` |
| 扫描缓存的键（`torsion/mlp.py: scan_cache_key`）只含模型、四元组、`info`（总电荷、spin）、扫描参数、原子序数、键；与部分电荷无关，换电荷可以复用 | 读代码 |
| CLI `derive-torsions` 目前不保存扫描（只有 Gate 保存） | 读代码 `cli.py` |

## 1. 全局约束

MVP-0/1/2 计划 §1 的全部条目照常适用（openmm_dev、不装包、float64、OpenMM Reference、`-m external` 分层、日志、`RunContext`、种子）。另外：

1. **RESP 本期不做**（决策 0018 第 1 条）。`ChargeProvider` 接口要能容纳以后的 RESP（provider 可以带外部输入文件），但不写 RESP 代码、不写占位类。
2. **不调 QM 程序**，不碰远程机器。
3. **外部程序与模型**：antechamber 调用的测试标 `external`；NAGL 会 import torch，测试也标 `external`。`charge/nagl.py` 只在函数内 import openff；`charge/` 其余模块与单元测试不得 import openff 或 torch。
4. **工作目录**：AM1-BCC 的中间文件写进调用方给的 `workdir`（CLI 默认 `<outdir>/am1bcc/`，Gate 写进运行目录），不用 `/tmp`，不用 `tempfile`。测试用 pytest 的 `tmp_path`。
5. **不加门**（全局规则第 8 章）：除了本计划写明的错误（拓扑不一致、外部程序失败、NAGL 适用域外、重拟合缺少缓存或缺少 `replaced` 记录），不设其他拒跑条件。电荷"好不好"只报告。
6. **不用 hash 判断能否复用**；扫描缓存沿用 MVP-2 的逐项比较。
7. **不做**：RESP；构象选择（AM1-BCC 用 antechamber 默认做法，与固件一致）；换电荷后对全局 Hessian 拟合（设计 §6.4.2）的 bonded 参数重拟合，只在 provenance 里记录电荷来源；LJ 不动（设计 §15 Phase 1）。

---

## 2. 锁定的定义

### 2.1 残差均摊（与固件规则相同）

输入：原始电荷 \(q\)（长度 N，float64），目标总电荷 \(Q\)（`ir.total_charge`）。
1. \(r = Q-\sum_i q_i\)。若 \(|r| \le 10^{-12}\)：原样返回。
2. 否则 \(q_i \leftarrow q_i + r/N\)；再令最后一个原子吸收剩余：\(q_{N-1} \leftarrow q_{N-1} + (Q-\sum_i q_i)\)。

这与 `scripts/make_fixtures.py::_equalize_mol2_charges` 的算术一致（该函数之后还把数值格式化成 10 位小数，所以和固件比较的容差取 1e-10，见 §6）。返回均摊后的电荷与 `r`。

### 2.2 写回 IR 与 provenance

- `ForceFieldIR.replace_terms` 增加关键字参数 `atoms: tuple[Atom, ...] | None`；电荷写回只改 `Atom.partial_charge`，其余字段不变。写回后 `validate(ir)` 必须通过。
- `provenance["charge"]`：`method`（`"loaded"`/`"am1bcc"`/`"nagl"`）、`raw_sum`、`residual`（§2.1 的 r）、`source`（loaded：`"ir"` 或文件路径；am1bcc：antechamber 路径与 AmberTools 版本行；nagl：模型文件名、大小、修改时间与 openff-nagl 版本）。
- MLP 二面角拟合时，`provenance["torsion"]` 新增 `charge_method`：拟合当时 `provenance["charge"]["method"]`，没有则为 `"prmtop"`。

### 2.3 三种 provider

接口：`ChargeProvider.compute(ir, positions_nm, mol) -> np.ndarray`（长度 N，单位 e，未均摊）。

| provider | 计算 | 拓扑检查（不一致 → `TopologyMismatchError`，复用 `bonded.assign`） |
|---|---|---|
| `LoadedCharges(source=None)` | `source=None`：取 IR 现有电荷。`.mol2`：读 `@<TRIPOS>ATOM` 第 9 列。其他后缀：文本文件，每行一个数（忽略空行和 `#` 开头的行） | 原子数相等；mol2 还要求逐原子 atom name 与 IR 相同 |
| `AM1BCCCharges(workdir)` | 把 `mol`（坐标换成 `positions_nm`×10 Å）写成 `workdir/input.sdf`，在 `workdir` 里运行 `antechamber -i input.sdf -fi sdf -o am1bcc.mol2 -fo mol2 -c bcc -at gaff2 -nc <ir.total_charge> -pf y`（命令行写进 `workdir/command.txt`，输出进 `workdir/antechamber.log`），读输出 mol2 第 9 列 | 原子数相等；输出 mol2 的元素序列（atom name 去掉数字，不区分大小写，如 `Cl1` → Cl）与 IR 的 `element` 相同 |
| `NAGLCharges(model="openff-gnn-am1bcc-1.0.0.pt")` | `Molecule.from_rdkit(mol, allow_undefined_stereo=True, hydrogens_are_explicit=True)`，`GNNModel.load(<模型路径>, eval_mode=True).compute_property(m, check_domains=True)`；模型路径 = `openff/nagl_models/models/am1bcc/<model>` | 原子数与逐原子原子序数与 IR 相同 |

错误：antechamber 返回非零或没有产出 mol2 → `ChargeError`（带 `antechamber.log` 末 20 行）；NAGL 适用域外 → `ChargeError`（带 openff 的原始报错）。

### 2.4 二面角重拟合（决策 0018 第 3 条）

- `MLPTorsionProvider.parameterize` 替换某根键的项时，把被删掉的原项按键记进报告行 `replaced_terms`（每项：i, j, k, l, periodicity, phase, force_constant, source, meta），并随 `report.to_dict()` 进入 `provenance["torsion"]["rotors"]`。兜底的键不记。
- `refit_torsions(ir, positions_nm, mol, *, scan_cache, model_path=None) -> tuple[ForceFieldIR, TorsionReport]`：
  1. IR 里没有 `source == "mlp"` 的项 → 原样返回（report 为空）。
  2. 对每个 `status == "mlp"` 的键：删掉该键上全部 proper torsion，放回 `replaced_terms`，得到基线 IR。任何一个 mlp 键缺少 `replaced_terms` → `TorsionRefitError`（提示从 GAFF2 IR 重跑 `derive-torsions`）。
  3. 用 provenance 里的 `model`、`scan`（扫描参数）构造 `MLPTorsionProvider(..., scan_cache=scan_cache, model_path=model_path, cache_only=True)`，对基线 IR 调用 `parameterize`。`model_path=None` 时取 `provenance["torsion"]["model_file"]["path"]`（没有记录则为 None，与原拟合一致），保证缓存键与原拟合相同。
  4. 原来是 mlp 的键若变成 `no_cached_scan` → `TorsionRefitError`（列出键）；不静默退回 GAFF2。
- `assign_charges` 写完电荷后，若 IR 有 mlp 项：给了 `scan_cache` → 调 `refit_torsions`；没给 → `TorsionRefitError`。

**自洽（核心）**：电荷不变时重拟合，输出的 `proper_torsions` 与重拟合前逐项相等（`==`，同一缓存、同一单遍流程、同一 IR 浮点数）。

---

## 3. 文件结构

```text
src/ffcompiler/ir/schema.py           # T31：replace_terms(atoms=...)
src/ffcompiler/charge/
├── __init__.py
├── base.py        # T31：ChargeProvider、ChargeReport、ChargeError、equalize_residual
├── apply.py       # T31（写回）/ T35（触发重拟合）：assign_charges
├── loaded.py      # T32
├── am1bcc.py      # T33
└── nagl.py        # T34（唯一 import openff 的模块）
src/ffcompiler/torsion/mlp.py         # T35：replaced_terms、charge_method、refit_torsions、TorsionRefitError
src/ffcompiler/cli.py                 # T36：assign-charges；derive-torsions --scan-cache / --ir-in
scripts/gate_mvp3.py                  # T37
tests/unit/test_charge_*.py、tests/unit/test_torsion_refit.py、tests/unit/test_gate_mvp3_report.py
```

---

## 4. Review Focus

1. **原子顺序不一致**：电荷文件 / SDF / mol2 与 IR 原子数、名字或元素对不上，必须报 `TopologyMismatchError`，不能把电荷错位写进 IR。→ T32、T33、T34。
2. **带电分子**：醋酸根（−1）、甲铵（+1）。antechamber 的 `-nc` 取 `ir.total_charge`；三种 provider 写回后总电荷都等于 Q。→ T31、T33、T34。
3. **换电荷但不能重拟合**：有 mlp 项却没给缓存、缓存缺键、缺少 `replaced_terms`（例如 G2 产出的旧 IR），必须明确报错，不能静默保留旧二面角或退回 GAFF2。→ T35。
4. **外部失败**：antechamber 失败、NAGL 适用域外（含 Si 的分子），报 `ChargeError` 并带原始信息。→ T33、T34。
5. **重复操作**：同一 IR 连续两次 `assign_charges(LoadedCharges())`，电荷与二面角都不变（幂等）。→ T35。

---

## 5. 任务

### T31：电荷写回框架

**Files:** Modify `src/ffcompiler/ir/schema.py`；Create `src/ffcompiler/charge/__init__.py`、`base.py`、`apply.py`、`tests/unit/test_charge_apply.py`

**Interfaces（Produces）：**
- `ForceFieldIR.replace_terms(..., atoms: tuple[Atom, ...] | None = None, ...)`
- `class ChargeError(RuntimeError)`
- `class ChargeProvider(ABC)`：`method: str`（类属性）；`compute(self, ir, positions_nm, mol) -> np.ndarray`；`source_info(self) -> dict`
- `equalize_residual(q: np.ndarray, total_charge: int) -> tuple[np.ndarray, float]`（§2.1）
- `@dataclass ChargeReport: method: str, raw_sum: float, residual: float, source: dict, torsion_refit: dict | None = None`；`to_dict()`
- `assign_charges(ir, provider, positions_nm, mol, *, scan_cache: Path | None = None) -> tuple[ForceFieldIR, ChargeReport]`（T31 只写电荷与 provenance；T35 加重拟合）

**必测用例：**
1. `equalize_residual`，合成数据：N=9、Q=0、Σq=1e-3 → 均摊后 \(|\sum q - Q| \le 10^{-15}\)，且前 N−1 个原子各加 r/N；\(|r| \le 10^{-12}\) 时原样返回。与固件规则的一致性由 T33 用例 1 检验。
2. 用一个返回固定数组的假 provider：写回后 `partial_charge` 逐原子等于均摊结果，其他 Atom 字段、全部项、`total_charge` 不变；`validate` 通过；`provenance["charge"]` 字段齐全。
3. 醋酸根（Q=−1）、甲铵（Q=+1）：假 provider 给出和为 Q+2e-7 的电荷 → 写回后 \(|\sum q - Q| \le 10^{-10}\)。

### T32：`LoadedCharges`

**Files:** Create `src/ffcompiler/charge/loaded.py`、`tests/unit/test_charge_loaded.py`

**Interfaces（Produces）：** `LoadedCharges(source: Path | None = None)`，`method = "loaded"`

**必测用例：**
1. `source=None` → 返回 IR 现有电荷（乙醇固件）。
2. 固件 mol2（`tests/fixtures/gaff2/ethanol/ethanol.mol2`）→ 与 IR 现有电荷逐原子差 ≤ 1e-8 e（IR 电荷来自 prmtop 的 q×18.2223，E16.8 格式有 9 位有效数字，|q×18.2223| ≤ 20 时舍入 ≤ 5e-8/18.2223 ≈ 3e-9 e；实测乙醇 1.9e-9、醋酸根 8.8e-10、苯甲酸甲酯 5.5e-10）。
3. 文本文件（含空行、`#` 注释）读入正确。
4. 少一行的文本文件、atom name 改掉一个的 mol2 → `TopologyMismatchError`。

### T33：`AM1BCCCharges`

**Files:** Create `src/ffcompiler/charge/am1bcc.py`、`tests/unit/test_charge_am1bcc.py`（`external`）

**Interfaces（Produces）：** `AM1BCCCharges(workdir: Path)`，`method = "am1bcc"`；`source_info()` 含 antechamber 路径与 `antechamber.log` 里的版本行。

**必测用例（external）：**
1. **12 个固件**：`assign_charges(ir, AM1BCCCharges(tmp_path/name), ...)` 后与固件 mol2 电荷逐原子差 ≤ 1e-10 e（§6 G3.1 的依据；冒烟实测 0）。参数化 12 个分子。
2. 醋酸根、甲铵：antechamber 命令行里的 `-nc` 等于 `ir.total_charge`（读 `workdir` 里记录的命令行）。
3. 把 `FFC_ANTECHAMBER` 指向一个总是失败的脚本（`tmp_path` 里写一个 `exit 3` 的 shell 脚本）→ `ChargeError`，信息里有日志末尾。为此 `am1bcc.py` 用环境变量 `FFC_ANTECHAMBER`（默认 PATH 里的 `antechamber`）定位程序，与 `FFC_XTB`、`FFC_GMX` 一致。

### T34：`NAGLCharges`

**Files:** Create `src/ffcompiler/charge/nagl.py`、`tests/unit/test_charge_nagl.py`（`external`）

**Interfaces（Produces）：** `NAGLCharges(model: str = "openff-gnn-am1bcc-1.0.0.pt")`，`method = "nagl"`；`nagl_raw_charges(mol: Chem.Mol, model: str) -> np.ndarray`（不经 IR，供适用域测试）

**必测用例（external）：**
1. 12 个固件都能算；写回后 \(|\sum q - Q| \le 10^{-10}\)；与固件 AM1-BCC 的逐原子最大差 ≤ 0.1 e（这一条只用来抓原子错位：冒烟实测最大 0.040 e，错位会把 H 与重原子的电荷对调，差在 0.3 e 以上）。
2. 四甲基硅烷（RDKit 由 SMILES `C[Si](C)(C)C` 加氢、嵌入）→ `nagl_raw_charges` 抛 `ChargeError`。
3. 把 mol 的原子顺序打乱（`Chem.RenumberAtoms`）后用原 IR → `TopologyMismatchError`。

### T35：二面角 `replaced_terms` 与重拟合

**Files:** Modify `src/ffcompiler/torsion/mlp.py`、`src/ffcompiler/charge/apply.py`；Create `tests/unit/test_torsion_refit.py`

**Interfaces：**
- Consumes：T31 `assign_charges`、T32 `LoadedCharges`；MVP-2 `MLPTorsionProvider`、`OpenMMCalculator`。
- Produces：`class TorsionRefitError(RuntimeError)`；`refit_torsions(ir, positions_nm, mol, *, scan_cache: Path, model_path: Path | None = None, calculator_factory=None) -> tuple[ForceFieldIR, TorsionReport]`（`calculator_factory` 只给测试用；cache_only 下不会被调用）。

**必测用例（参考势 = GAFF2 乙醇的 `OpenMMCalculator`，`model_name="mm-self"`，单元测试，不用 MACE）：**
1. `parameterize` 后两根 mlp 键的报告行都有 `replaced_terms`，内容与输入 IR 该键上的项逐项相等。
2. **自洽**：`assign_charges(ir_mlp, LoadedCharges(), ..., scan_cache=缓存)` → `proper_torsions` 与 `ir_mlp` 逐项相等（`==`）；再做一次仍相等（Review Focus 5）。
3. 用文本文件把两个原子的电荷各 ±0.05 e（总和不变）→ 重拟合后两根键仍是 mlp，二面角项与原来不同；`provenance["torsion"]["charge_method"] == "loaded"`。
4. 有 mlp 项但 `scan_cache=None` → `TorsionRefitError`；缓存目录里删掉一根键的 `.npz` → `TorsionRefitError` 并列出该键；删掉 provenance 里的 `replaced_terms` → `TorsionRefitError`。

### T36：CLI

**Files:** Modify `src/ffcompiler/cli.py`、`tests/unit/test_cli.py`

- `ffcompiler assign-charges --prmtop P --inpcrd C --sdf S [--ir-in F] --method {loaded,am1bcc,nagl} [--charges-file F] [--scan-cache D] --output {openmm,amber,gromacs}... --outdir D [--ir-json F]`：写各引擎文件、IR JSON、`charge_report.json`。`--ir-in` 给了就从 IR JSON 读 IR（坐标仍取 `--inpcrd`）。am1bcc 的工作目录为 `<outdir>/am1bcc/`。
- `derive-torsions` 增加 `--scan-cache D`（默认 `<outdir>/scans`）与 `--ir-in F`。

**必测用例：**
1. `assign-charges --method loaded --charges-file <文本>`（乙醇，非 external）：`charge_report.json` 与 IR JSON 里的电荷一致，退出码 0。
2. `derive-torsions --reference mm --scan-cache D` 后 D 里有每根可转键的 `.npz`；接着 `assign-charges --ir-in <该 IR> --method loaded --scan-cache D` 退出码 0，二面角不变。
3. `--ir-in` 指向有 mlp 项的 IR 但不给 `--scan-cache` → 退出码 1，stderr 有 `TorsionRefitError` 的提示。

### T37：Gate G3、文档

**Files:** Create `scripts/gate_mvp3.py`、`tests/unit/test_gate_mvp3_report.py`；Modify `CHANGELOG.md`、三语 README、`HANDOVER.md`、roadmap、决策索引

`python scripts/gate_mvp3.py [--scan-cache DIR] [--skip-compare] [--skip-pytest] [--molecules ...]`（后三个只用于调试），`RunContext("gate_mvp3")`，默认 `--scan-cache jobs/gate_mvp2_20261008_3/scans`（只读，`cache_only`）。
1. G3.1：12 个固件跑 `AM1BCCCharges`，与固件 mol2 比较。
2. G3.2、G3.3：11 个 Gate 分子（乙腈按决策 0013 不进 compare）× 3 种 provider（`LoadedCharges()`、AM1-BCC、NAGL）写回，记录总电荷偏差，跑 MVP-0 `compare`，按 `gate_mvp2.judge_compare` 判定（决策 0017），每个比较存 `compare/<provider>/<分子>/report.json`。
3. G3.4：有可转键的 6 个分子（共 12 根键）：GAFF2 IR 先换 AM1-BCC 电荷 → 用缓存派生 MLP 二面角（`cache_only`）→ 再跑一次 AM1-BCC（新的 workdir）并重拟合 → 二面角与派生结果逐项相等（同时检验 antechamber 的确定性）；从同一个 MLP IR 换 NAGL 电荷 → 12 根键全部 mlp。注意 IR 的初始电荷来自 prmtop，与固件 mol2 差约 1e-9，所以不能拿 GAFF2 IR 直接当"电荷不变"的基线。
4. G3.5：pytest（复用 `gate_mvp1._pytest_summary`）。
5. 只报告：每个分子 NAGL 与 AM1-BCC 的逐原子最大差、RMS；NAGL 电荷下重拟合后的剖面 RMSE（指标与 G2 相同）。

**必测用例：** 两行假结果生成报告，`report.md` 表头固定、`acceptance.json` 字段齐全；任一判据缺数据时为 False（fail-closed）。

文档：Gate 通过后更新 CHANGELOG、三语 README（功能表加电荷一行、验收表加 G3）、HANDOVER、roadmap 状态列、决策索引加 0018。

---

## 6. Gate G3 验收判据

| # | 判据 | 容差与依据 |
|---|---|---|
| G3.1 | 12 个固件：AM1-BCC provider 与固件 mol2 电荷逐原子一致 | ≤ 1e-10 e：固件 mol2 电荷写到 10 位小数，舍入误差 5e-11；冒烟实测差 0。只在 node-A 上测过，跨机器未测 |
| G3.2 | 三种 provider 写回后总电荷守恒 | \(|\sum q - Q| \le 10^{-10}\) e：均摊后只剩 float64 求和误差（约 N·1e-16），留足余量 |
| G3.3 | 3 种 provider × 11 个分子的 IR 通过跨引擎比较 | 判定按决策 0017（以 OpenMM 为准，sander 腿只报告） |
| G3.4 | 二面角重拟合：AM1-BCC 重跑一次后逐项相等；NAGL 电荷下 12 根键全部 mlp | 逐项相等（`==`） |
| G3.5 | pytest 全部通过；skip 只允许是 `gmx` 不存在或 CUDA 不可用 | 同 G2.5 |

NAGL 与 AM1-BCC 的差、换电荷后的二面角剖面指标**只报告、不设门槛**。

### G3 正式通过（2026-10-10，`jobs/gate_mvp3_20261010_2/`，node-A，650 s）

- G3.1–G3.5 全部通过。G3.1 最大差 4.3e-11 e（12 个固件）；G3.2 最大 8.3e-17 e（33 组）；G3.3 33/33；G3.4 6 个分子 12 根键：AM1-BCC 重跑后逐项相等，NAGL 电荷下全部 mlp；pytest 非 external 436、external 59，0 失败 0 跳过。sander 没有超差行。
- 只报告：NAGL 与 AM1-BCC 逐原子最大差 4.3e-5–0.040 e；换成 NAGL 电荷后剖面 RMSE 与 AM1-BCC 下相近（methyl_benzoate 1-2：3.10 vs 3.36 kJ/mol，仍是已知缺陷）。
- 执行中的偏离（审阅 Important 1）：§2.3 的拓扑检查从"只比元素"加强为"原子序数 + 键集合"（`torsion.rotors.check_mol_topology`），AM1-BCC 在调用 antechamber 前检查；`--ir-in` 与 `--prmtop/--sdf` 核对（审阅 Important 2）。执行记录与裁定见 `.superpowers/sdd/2026-10-10-ffc-mvp3-charge-providers/progress.md`。

## 7. 待用户确认

（无。范围与耦合已由决策 0018 定下；G3 判据 2026-10-10 用户确认。）
