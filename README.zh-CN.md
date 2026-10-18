# Crucible

[English](README.md) | **中文** | [日本語](README.ja.md)

引擎无关的小分子力场编译器。已参数化的分子先读进统一的中间表示 `ForceFieldIR`，
可选地用量子化学 Hessian 重新推导参数，再原生写出 **GROMACS**、**AMBER**、**OpenMM** 三种格式。

**版本 0.2.0**：MVP-0（IR + 三后端导出）、MVP-1（由 Hessian 推导 bond/angle 参数）、MVP-2（MLP 二面角，
以 OpenMM 为准，决策 0017）与 MVP-3（电荷 provider：载入、AM1-BCC、NAGL）均已完成；2026-10-18 在当前
源码上重跑 G0–G3，全部通过。2026-10-13 以来代码还新增：ORCA `.hess` 读取器、跨命令共享的 IR + 坐标
管线（`--ir-in` / `--coords-in`）、`global-fit` 与 `check-staleness` 命令、Hessian fallback 框架
（MVP-5，仅 Python API）、learned BCC 骨架（MVP-4，无训练模型），以及可选的顺序二面角重拟合
（`derive-torsions --refinement-rounds {0,1,2,3}`，默认 0 = 单遍）。2026-10-18 对照设计稿的结案验收
（34 条，现行偏离清单）见
[ACCEPTANCE_20261018_CLOSURE](docs/history/ACCEPTANCE_20261018_CLOSURE.md)；变更见
[CHANGELOG](CHANGELOG.md)。

## 功能

| 模块 | 内容 |
|---|---|
| IR | `ForceFieldIR`：原子、bond、angle、proper/improper 二面角、LJ 与电荷、排除表与 1-4 对；校验与 NumPy 参考能量 |
| 导入 | AMBER prmtop/inpcrd + SDF（取形式电荷），原生读取，不依赖 ParmEd |
| 导出 | GROMACS top/itp/gro、AMBER prmtop/inpcrd、OpenMM System/XML，全部原生写出 |
| Hessian | Gaussian fchk、ORCA `.hess`；xtb 6.7.1（GFN-FF、GFN2），带收敛、电荷与质量检查；GFN-FF → GFN2 → DFT 文件的 fallback 框架（仅 Python API） |
| 参数推导 | mSeminario、m2Seminario、DRIH、全局 Hessian 拟合；按分子图对称性平均 |
| 二面角（MVP-2） | 可转键识别、MACE-POLAR-1-M（经 ASE）约束松弛扫描、扣除 MM、拟合成 AMBER 周期项；失败退回 GAFF2 并记录原因；可选顺序重拟合轮数（`ffcompiler derive-torsions`） |
| 电荷（MVP-3） | 载入电荷、AM1-BCC（antechamber/sqm）、NAGL（openff-gnn-am1bcc-1.0.0）；总电荷严格守恒；换电荷时用缓存的扫描重拟合 MLP 二面角（`ffcompiler assign-charges`）。RESP 尚未做 |
| 管线 | IR 与坐标 sidecar 在命令之间传递（`--ir-in` / `--coords-in`）；记录拟合依赖并检查是否过期（`ffcompiler check-staleness`） |
| 验证 | 跨引擎能量/力比较（NumPy、OpenMM ×4、sander ×2、GROMACS）与四道验收 Gate（G0–G3） |

## 环境要求

使用现有科学计算环境，本项目不安装任何东西。

- Python ≥ 3.12，需要 `numpy`、`scipy`、`openmm`、`rdkit`；MLP 二面角与 NAGL 电荷功能另需
  ASE + MACE 和 OpenFF NAGL（模型 `openff-gnn-am1bcc-1.0.0`）
- 外部程序（只在用到对应功能或测试时需要）：
  - `xtb` 6.7.1：默认 `/home/user/software/xtb-dist/bin/xtb`，可用 `FFC_XTB` 替换
  - GROMACS：按 `FFC_GMX`、`PATH`、`/home/user/gmx26.3/bin/gmx` 的顺序查找；mixed 精度即可
  - AmberTools `sander` / `tleap`：跨引擎检查与固件生成；AM1-BCC 电荷另用 `antechamber`
    （`FFC_ANTECHAMBER` 或 `PATH`）
  - 二面角扫描用的 MACE-POLAR-1-M 模型：默认 `/home/user/MLP/mace/MACE-POLAR-1-M.model`，可用 `FFC_MACE_MODEL` 替换
  - Sobtop 2026.1.16：只在重新生成参考数据时使用（`scripts/make_sobtop_refs.py`）

```bash
source scripts/env.sh
```

## 快速上手

把 AMBER 参数化的分子转换到三种引擎：

```bash
PYTHONPATH=src python -m ffcompiler.cli convert \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --output gromacs amber openmm --outdir jobs/ethanol_conversion
```

用 GFN2 Hessian 重新推导 bond 和 angle（xtb 在 `<outdir>/xtb/` 下运行）：

```bash
PYTHONPATH=src python -m ffcompiler.cli derive-bonded \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --hessian-provider gfn2 --method mseminario \
  --output gromacs amber openmm --outdir jobs/ethanol_mseminario
```

- `--hessian-provider {gfnff,gfn2,fchk,orca-hess}`（`fchk`、`orca-hess` 需要 `--hessian-file`）
- `--method {mseminario,m2seminario,drih}`
- `--no-symmetrize`：不对对称等价项取平均
- `--linear-angles {error,library}`：默认 `error`；`library` 对 > 175° 的键角保留原参数

用 MACE-POLAR 松弛扫描拟合可转键二面角（扫描缓存在 `<outdir>/scans/` 下）：

```bash
PYTHONPATH=src python -m ffcompiler.cli derive-torsions \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --output openmm --outdir jobs/ethanol_torsions
```

- `--step-deg`（默认 10°）、`--device`（默认 cpu）、`--scan-cache`（默认 `<outdir>/scans`）
- `--refinement-rounds {0,1,2,3}`：单遍之后再做若干轮顺序更新，每根键按其他键的当前拟合项重新拟合；
  默认 0 保持单遍。换电荷（`assign-charges`）时按记录的协议用缓存扫描重放。

检查导出文件在每个引擎里都能复现 IR 能量：

```bash
PYTHONPATH=src python -m ffcompiler.cli check-equivalence \
  --prmtop ... --inpcrd ... --sdf ... --workdir jobs/ethanol_check
```

在同一个 IR 上串联：每一步写出 `ir.json` 和 `coords.json`（本步实际使用或返回的几何，例如
Hessian 几何）；下一步用 `--ir-in` / `--coords-in` 读入（原子序数、元素、原子名必须一致，不做重排）：

```bash
# 1) 电荷
PYTHONPATH=src python -m ffcompiler.cli assign-charges --prmtop P --inpcrd C --sdf S \
  --method am1bcc --output openmm --outdir jobs/mol/q
# 2) 由 Hessian 推导 bond/angle，接着第 1 步
PYTHONPATH=src python -m ffcompiler.cli derive-bonded \
  --ir-in jobs/mol/q/ir.json --coords-in jobs/mol/q/coords.json \
  --hessian-provider gfn2 --method mseminario --output openmm --outdir jobs/mol/b
# 3) 导出最终 IR + 坐标
PYTHONPATH=src python -m ffcompiler.cli convert \
  --ir-in jobs/mol/b/ir.json --coords-in jobs/mol/b/coords.json \
  --output gromacs amber openmm --outdir jobs/mol/final
# 拟合记录的输入（电荷、二面角等）之后被改过 → 退出码 1
PYTHONPATH=src python -m ffcompiler.cli check-staleness --ir-in jobs/mol/b/ir.json
```

`global-fit` 提供 QMDFF 式全局 Hessian 拟合（Hessian 选项相同；另有 `--prior`、`--alpha`）；
有等价类落在 k = 0 下界时写出 `ir_failed_lower_bound.json`，不导出，退出码 1。


在 Python 中使用：

```python
from ffcompiler.importers import from_amber
from ffcompiler.hessian.xtb import run_xtb_hessian
from ffcompiler.bonded import derive_bonded

ir, pos = from_amber("mol.prmtop", "mol.inpcrd", "mol.sdf")
h = run_xtb_hessian([a.atomic_number for a in ir.atoms], pos,
                    method="gfn2", charge=ir.total_charge, workdir="jobs/xtb")
new_ir, geom, report = derive_bonded(ir, h, method="drih", input_positions_nm=pos)
```

## 约定

- **IR 单位**：nm、rad、kJ mol⁻¹、元电荷。
- **谐振项**统一按 E = ½ k Δq² 储存，即 k 就是 Hessian 对角元。GROMACS、OpenMM 直接用 k；
  AMBER 用 K = k/2（E = K Δq²），再换算 kcal 与 Å。
- **DRIH** 内部使用原子单位（Bohr、E_h/a₀²）：stretch 对角元按 E_h/a₀² ≈ 9.37583 × 10⁵
  kJ mol⁻¹ nm⁻² 换算，angle 对角元按 E_h/rad² = 2625.50 kJ mol⁻¹ rad⁻² 换算（不再除 a₀²）。
- 支持的 Hamiltonian：气相、无截断、Lorentz–Berthelot 混合、统一的 AMBER 1-4 缩放。
  其他情况在导入时直接拒绝。

## 验收状态

2026-10-18 在当前源码上（接入可选顺序重拟合之后）重跑了四个 Gate，全部通过，数值与 2026-10-13
那轮一致；每个 Gate 的 pytest 都是 524 + 59。Gate 判据走默认单遍；顺序模式由单元测试和逐位重放
核对覆盖（`experiments/20261018_refinement_verify/`）。

| Gate | 运行目录 | 结果 |
|---|---|---|
| G0：导出等价性 | `jobs/gate_mvp0_20261018_2/` | 通过：11 个分子 × 4 个构型，8 条比较腿，2596 行，0 失败 |
| G1：Hessian 方法 | `jobs/gate_mvp1_20261018_2/` | 通过：三种方法在 12/12 个参考用例上与 Sobtop 一致；全局拟合自洽误差 6.0 × 10⁻⁹ |
| G2：MLP 二面角 | `jobs/gate_mvp2_20261018_2/` | 通过：12/12 根可转键由 MLP 拟合，0 兜底；sander 以外的跨引擎比较全部通过。已知局限（决策 0017）：sander 把接近 π 的拟合相位吸附成 π（nma dihed 2.3 × 10⁻⁴ > 2 × 10⁻⁴ kcal/mol），只报告 |
| G3：电荷 provider | `jobs/gate_mvp3_20261018_2/` | 通过：AM1-BCC 在 12 个固件上复现固件电荷（最大差 4.3 × 10⁻¹¹ e）；3 种 provider × 11 个分子总电荷守恒（≤ 8.3 × 10⁻¹⁷ e），跨引擎比较全部通过；换电荷后重拟合 MLP 二面角：电荷不变时逐项相等，NAGL 电荷下 12/12 根键仍为 MLP |

推导出的力场相对参考 Hessian 的频率 RMSE（11 个分子的中位数，cm⁻¹；只报告、不设门槛）：

| Hessian | GAFF2 | mSeminario | m2Seminario | DRIH | 全局拟合 |
|---|---|---|---|---|---|
| GFN2 | 110 | 109 | 188 | 219 | 118 |
| GFN-FF | 101 | 95 | 76 | 203 | **57** |

自己跑 Gate（每次写一个新的 `jobs/gate_*_<日期>[_N]/` 目录，旧目录保留）：

```bash
python scripts/gate_mvp0.py
python scripts/gate_mvp1.py --jobs 6   # 6 个 xtb 并行 × 每个 4 线程
python scripts/gate_mvp2.py --jobs 8 --threads 4 [--resume 旧G2运行目录]   # MACE 扫描（给了 --resume 就复用其 scans/）
python scripts/gate_mvp3.py --scan-cache jobs/<本次G2运行>/scans          # 必填：本次验收的 G2 的 scans/
PYTHONPATH=src python -m pytest -m "not external" -q   # 单元测试（不调用外部程序）
```

## 已知限制

- **已知缺陷（MVP-2）**：苯甲酸甲酯 1-2 键上，MLP 拟合出的二面角对 MACE 剖面的复现比 GAFF2 还差
  （剖面 RMSE 3.36 vs 2.70 kJ/mol，拟合 RMSE 0.60），12 根可转键里只有这一根。2026-10-17 已诊断
  （`jobs/torsion_diagnosis_20261017/report.md`）：逐键拟合时其余键仍按 GAFF2 扣除，1-2 的拟合项
  吸收了偏硬的 GAFF2 甲基（0-1）项；0-1 换掉后补偿落空。`--refinement-rounds 3` 可降到
  0.44 kJ/mol，合并残差 0.83 → 0.22（样本内，`jobs/torsion_b3_sequential_20261017/report.md`）；
  所有更新顺序收敛到同一不动点，第 3 轮后顺序带来的差 ≤ 1e-4 kJ/mol
  （`jobs/torsion_b3_order_20261018/report.md`）。该模式保持可选：默认仍是单遍，B4 留出验证未做，
  顺序模式不在 Gate 判据里。
- **扫描质量（MVP-2）**：nma 1-3（酰胺 C–N）正反扫滞后最大 57 kJ/mol（36 帧中 14 帧超过 1 kJ/mol，
  集中在 ±60–120° 能垒区）。这是扫描问题而非拟合问题；是否重扫未定，也没有设自动兜底阈值。
- 同一根键的二面角能量全部放在一个四元组上（MVP-2 计划 §7 第 3 条）；其适用范围尚未界定。
- ORCA `.hess` 读取器只用 ORCA 6.1.1 的一个闭壳层样例核对过（频率相对差 5.7 × 10⁻⁸）。
- 线性键角（> 175°）不推导；乙腈只作诊断用例。
- m2Seminario 在 GFN-FF Hessian 上可能给出负的环角力常数（环丙烷），Sobtop 结果相同。
  `derive_bonded` 会抛 `NonPositiveForceConstantError`，不静默回退；这种情况请改用 mSeminario 或 DRIH。
- mSeminario 的 bond 力常数取决于键的书写方向（与 Sobtop 一致）；默认开启的对称化会消除由此带来的不对称。
- 导出只写单个残基；周期体系与其他组合规则不在范围内。

## 目录结构

```text
src/ffcompiler/
  ir/          ForceFieldIR 结构、校验、NumPy 能量、JSON 读写、坐标 sidecar
  importers/   AMBER prmtop 导入
  exporters/   GROMACS、AMBER、OpenMM 导出
  formats/     原生 prmtop / rst7 读写
  hessian/     CartesianHessian、fchk / ORCA .hess / xtb provider、fallback 阶梯
  bonded/      内坐标、mSeminario、m2Seminario、DRIH、对称化、全局拟合、依赖过期检查
  torsion/     可转键、MACE-POLAR 松弛扫描、扣除 MM、Fourier 拟合、MLP provider
  charge/      电荷 provider（载入、AM1-BCC、NAGL）、learned BCC 骨架、写回与二面角重拟合
  validation/  跨引擎比较、GROMACS 适配器、振动工具、pytest 统计
scripts/       Gate、固件生成、Sobtop 参考数据
tests/         单元与跨引擎测试、固件
docs/          计划、决策记录、理论文档
jobs/          正式运行（只跟踪 manifest/acceptance/report）
```

## 文档

- [docs/history/ACCEPTANCE_20261018_CLOSURE.md](docs/history/ACCEPTANCE_20261018_CLOSURE.md)：
  当前版本的结案验收（34 条对照；第 2 节为已知缺陷、待拍板项与未做的验证）
- [docs/theory/hessian_methods.md](docs/theory/hessian_methods.md)：推导、单位与约定
- [docs/decisions/](docs/decisions/README.md)：决策记录 0001–0021
- [docs/superpowers/plans/](docs/superpowers/plans/)：路线图与各期计划
- [HANDOVER.md](HANDOVER.md)：当前已验证状态与下一步
- [CHANGELOG.md](CHANGELOG.md)

## 路线图

MVP-2 MLP 二面角：已完成，Gate G2 通过。MVP-3 电荷 provider：已完成，Gate G3 通过（RESP 往后放，
决策 0018）。MVP-4 学习式 BCC：只有 provider 骨架与立项草案，训练未立项。MVP-5：fallback 框架已就位，
调度研究未开始。按结案报告第 2 节，下一步：B4 留出验证、真实 MACE 参考的 C5 正式组合验收、MVP-4
立项讨论；再往后是 molecule/typing 入口与 preset/单命令 CLI（设计 §1/§17/§22）。待拍板项——默认
bonded 方法、默认重拟合轮数、版本号——都列在那里。仓库以 CRUCIBLE 私有托管。
