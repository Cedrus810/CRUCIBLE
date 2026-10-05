# Crucible

[English](README.md) | **中文** | [日本語](README.ja.md)

引擎无关的小分子力场编译器。已参数化的分子先读进统一的中间表示 `ForceFieldIR`，
可选地用量子化学 Hessian 重新推导参数，再原生写出 **GROMACS**、**AMBER**、**OpenMM** 三种格式。

**版本 0.2.0**：MVP-0（IR + 三后端导出）与 MVP-1（由 Hessian 推导 bond/angle 参数）
已完成，Gate G0 与 Gate G1 均已通过。变更见 [CHANGELOG](CHANGELOG.md)。

## 功能

| 模块 | 内容 |
|---|---|
| IR | `ForceFieldIR`：原子、bond、angle、proper/improper 二面角、LJ 与电荷、排除表与 1-4 对；校验与 NumPy 参考能量 |
| 导入 | AMBER prmtop/inpcrd + SDF（取形式电荷），原生读取，不依赖 ParmEd |
| 导出 | GROMACS top/itp/gro、AMBER prmtop/inpcrd、OpenMM System/XML，全部原生写出 |
| Hessian | Gaussian fchk；xtb 6.7.1（GFN-FF、GFN2），带收敛、电荷与质量检查 |
| 参数推导 | mSeminario、m2Seminario、DRIH、全局 Hessian 拟合；按分子图对称性平均 |
| 验证 | 跨引擎能量/力比较（NumPy、OpenMM ×4、sander ×2、GROMACS）与两道验收 Gate |

## 环境要求

使用现有科学计算环境，本项目不安装任何东西。

- Python ≥ 3.12，需要 `numpy`、`scipy`、`openmm`、`rdkit`
- 外部程序（只在用到对应功能或测试时需要）：
  - `xtb` 6.7.1：默认 `/home/user/software/xtb-dist/bin/xtb`，可用 `FFC_XTB` 替换
  - GROMACS：按 `FFC_GMX`、`PATH`、`/home/user/gmx26.3/bin/gmx` 的顺序查找；mixed 精度即可
  - AmberTools `sander` / `tleap`：跨引擎检查与固件生成
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

- `--hessian-provider {gfnff,gfn2,fchk}`（`fchk` 需要 `--hessian-file`）
- `--method {mseminario,m2seminario,drih}`
- `--no-symmetrize`：不对对称等价项取平均
- `--linear-angles {error,library}`：默认 `error`；`library` 对 > 175° 的键角保留原参数

检查导出文件在每个引擎里都能复现 IR 能量：

```bash
PYTHONPATH=src python -m ffcompiler.cli check-equivalence \
  --prmtop ... --inpcrd ... --sdf ... --workdir jobs/ethanol_check
```

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

| Gate | 运行目录 | 结果 |
|---|---|---|
| G0：导出等价性 | `jobs/gate_mvp0_20261005_3/` | 通过：11 个分子 × 4 个构型，8 条比较腿，2596 行，0 失败 |
| G1：Hessian 方法 | `jobs/gate_mvp1_20261005_4/` | 通过：三种方法在 12/12 个参考用例上与 Sobtop 一致；全局拟合自洽误差 6.0 × 10⁻⁹ |

推导出的力场相对参考 Hessian 的频率 RMSE（11 个分子的中位数，cm⁻¹；只报告、不设门槛）：

| Hessian | GAFF2 | mSeminario | m2Seminario | DRIH | 全局拟合 |
|---|---|---|---|---|---|
| GFN2 | 110 | 109 | 188 | 219 | 118 |
| GFN-FF | 101 | 95 | 76 | 203 | **57** |

自己跑 Gate（每次写一个新的 `jobs/gate_*_<日期>[_N]/` 目录，旧目录保留）：

```bash
python scripts/gate_mvp0.py
python scripts/gate_mvp1.py --jobs 6   # 6 个 xtb 并行 × 每个 4 线程
PYTHONPATH=src python -m pytest -m "not external" -q   # 单元测试（不调用外部程序）
```

## 已知限制

- 线性键角（> 175°）不推导；乙腈只作诊断用例。
- m2Seminario 在 GFN-FF Hessian 上可能给出负的环角力常数（环丙烷），Sobtop 结果相同。
  `derive_bonded` 会抛 `NonPositiveForceConstantError`，不静默回退；这种情况请改用 mSeminario 或 DRIH。
- mSeminario 的 bond 力常数取决于键的书写方向（与 Sobtop 一致）；默认开启的对称化会消除由此带来的不对称。
- 导出只写单个残基；周期体系与其他组合规则不在范围内。

## 目录结构

```text
src/ffcompiler/
  ir/          ForceFieldIR 结构、校验、NumPy 能量、JSON 读写
  importers/   AMBER prmtop 导入
  exporters/   GROMACS、AMBER、OpenMM 导出
  formats/     原生 prmtop / rst7 读写
  hessian/     CartesianHessian、fchk 与 xtb provider
  bonded/      内坐标、mSeminario、m2Seminario、DRIH、对称化、全局拟合
  validation/  跨引擎比较、GROMACS 适配器、振动工具、pytest 统计
scripts/       Gate、固件生成、Sobtop 参考数据
tests/         单元与跨引擎测试、固件
docs/          计划、决策记录、理论文档
jobs/          正式运行（只跟踪 manifest/acceptance/report）
```

## 文档

- [docs/theory/hessian_methods.md](docs/theory/hessian_methods.md)：推导、单位与约定
- [docs/decisions/](docs/decisions/README.md)：决策记录 0001–0015
- [docs/superpowers/plans/](docs/superpowers/plans/)：路线图与各期计划
- [HANDOVER.md](HANDOVER.md)：当前已验证状态与下一步
- [CHANGELOG.md](CHANGELOG.md)

## 路线图

MVP-2 MLP 二面角 · MVP-3 电荷 provider · MVP-4 学习式 BCC · MVP-5 fallback 调度，均未开始。
仓库以 CRUCIBLE 私有托管。
