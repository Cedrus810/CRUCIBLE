# Crucible

[English](README.md) | **中文** | [日本語](README.ja.md)

Crucible 是一个引擎无关的小分子力场编译器。一次参数化结果保存在统一的中间表示
`ForceFieldIR` 中，原生写出 **GROMACS**、**AMBER**、**OpenMM** 三种格式——每次导出都自动做
三引擎的能量/力等价校验。

```text
AMBER prmtop/inpcrd + SDF
        │  导入（原生读取，不依赖 ParmEd）
        ▼
   ForceFieldIR ◄──── 坐标 sidecar（所有命令：--ir-in / --coords-in）
        │
        ├─ bond/angle ← mSeminario · m2Seminario · DRIH · 全局 Hessian 拟合
        │                Hessian：Gaussian fchk · ORCA .hess · CP2K .out · xtb（GFN-FF / GFN2）
        ├─ 二面角     ← MACE-POLAR-1-M 约束松弛扫描 → Fourier 项（GAFF2 兜底，记录原因）
        └─ 电荷       ← 载入 · AM1-BCC（antechamber）· NAGL
        │
        ▼
   ForceFieldIR（完整——每个参数都带 provenance）
        │
   ┌────┼─────────┐
   ▼    ▼         ▼
GROMACS  AMBER   OpenMM   → 跨引擎能量 / 力比较
```

**版本 0.2.0。** MVP-0～MVP-3 已完成并有 Gate（2026-10-18 在当前源码上重跑 G0–G3，全部通过）。
2026-10-13 以来新增：ORCA `.hess` 读取器、CP2K `.out` 读取器、所有命令共享的 IR + 坐标管线、
`global-fit` 与 `check-staleness` 命令、周期盒支持、Hessian fallback 框架（MVP-5，仅 Python API）、
learned BCC 骨架（MVP-4，无训练模型），以及可选的顺序二面角重拟合
（`--refinement-rounds {0,1,2,3}`，默认 0 = 单遍）。对照设计稿的 34 条结案验收见
[ACCEPTANCE_20261018_CLOSURE](docs/history/ACCEPTANCE_20261018_CLOSURE.md)；项目的权威细粒度
状态在 [HANDOVER.md](HANDOVER.md)；变更见 [CHANGELOG](CHANGELOG.md)。

## 能力

### 一份 IR，三个引擎

`ForceFieldIR` 包含原子、bond、angle、proper/improper 二面角、LJ 与电荷、排除表与 1-4 对，
以及可选周期盒；自带校验与 NumPy 参考能量。三个导出器原生写出 GROMACS top/itp/gro、
AMBER prmtop/inpcrd、OpenMM System/XML——不经过 ParmEd，不经过中间格式。同一份 IR 写到三个
引擎必须复现相同的分项能量与原子力；`check-equivalence` 做这个比较（8 条腿：NumPy、
OpenMM ×4、sander ×2、GROMACS），G0 Gate 在 11 个分子 × 4 种几何上运行它。

### 由 Hessian 推导 bonded 参数

四种方法从 Cartesian Hessian 推导 bond/angle 参数，共用一个 provider 抽象和一套单位制
（内部用原子单位）：

- **mSeminario**——局部投影，通用默认候选；
- **m2Seminario**——强调对称性与方向一致的变体，适合刚性/高对称体系；
- **DRIH**——冗余内坐标对角元，适合刚性与特殊成键片段；
- **全局 Hessian 拟合**（QMDFF 式）——扣除固定部分（LJ + Coulomb + 1-4 + 二面角）后，
  全部力常数对整张 Hessian 做拟合，带先验正则，报告条件数、残差与 x₀ 处残余力；
- **localhess**（实验性）——局部柔量方法，仅做了内部一致性验证（决策 0023；不进 Gate）。

图对称等价项默认平均（`--no-symmetrize` 关闭）。非正力常数直接报错，不静默兜底。

### 由 MLP 推导二面角

从分子图识别可转键，对 MACE-POLAR-1-M（经 ASE）做约束松弛扫描，扣除当前力场的 MM 贡献，
剖面拟合成 AMBER 周期项。退回 GAFF2 的键逐个记录原因，重拟合时绝不静默升级原兜底键
（决策 0020）。可选的有界顺序重拟合（`--refinement-rounds {0,1,2,3}`）让每根键对其他键的
当前拟合重拟合；换电荷时按记录的协议从缓存扫描逐位重放。

### 电荷 provider

载入电荷、AM1-BCC（antechamber/sqm）、NAGL（`openff-gnn-am1bcc-1.0.0`）共用一个 provider
接口，总电荷严格守恒。换 provider 时用缓存扫描重拟合 MLP 二面角，不重扫。learned BCC 研究
骨架（按键反对称、电荷守恒的修正）已就位；尚无训练模型。RESP 推迟（决策 0018）。

### 感知层与力场参数库

几何成键感知层（Cordero 2008 半径 ×1.15 加数据文件距离窗口、带线性截断的角/二面角生成、
三斜 PBC 最近镜像、数据驱动的 improper 中心）与 Sobtop 格式 `LJ_param.dat` /
`bonded_param.dat` / `bondcrit.dat` 参数文件的读取/匹配器（X 通配符优先级、IDIVF、负周期、
LJ 换算）支撑参考数据流程；见决策 0022 与 `docs/theory/` 的逆向报告。

### Provenance 与过期检查

每个参数带 `source` 与元数据；各次拟合把输入（哪版电荷、哪版二面角、哪个 Hessian 身份）
记成逐项快照。`check-staleness` 在拟合记录过的任何输入被改动后以退出码 1 报错。IR + 坐标
sidecar 在命令间忠实传递几何——原子顺序、元素、名字永不重排。

## 快速开始

使用现有科学计算环境，本项目不安装任何东西。

- Python ≥ 3.12，需要 `numpy`、`scipy`、`openmm`、`rdkit`；MLP 二面角与 NAGL 电荷功能另需
  ASE + MACE 和 OpenFF NAGL
- 外部程序（只在用到对应功能或测试时需要）：
  - `xtb` 6.7.1：默认 `~/software/xtb-dist/bin/xtb`，可用 `FFC_XTB` 替换
  - GROMACS：按 `FFC_GMX`、`PATH`、`~/gmx26.3/bin/gmx` 的顺序查找；mixed 精度即可
  - AmberTools `sander` / `tleap`：跨引擎检查与固件生成；AM1-BCC 电荷另用
    `antechamber`（`FFC_ANTECHAMBER` 或 `PATH`）
  - MACE-POLAR-1-M 模型（二面角扫描）：默认 `~/MLP/mace/MACE-POLAR-1-M.model`，
    可用 `FFC_MACE_MODEL` 替换
  - Sobtop 2026.1.16：仅用于重新生成参考数据（`scripts/make_sobtop_refs.py`）

```bash
source scripts/env.sh
```

把已参数化的 AMBER 分子转换到三个引擎：

```bash
PYTHONPATH=src python -m ffcompiler.cli convert \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --output gromacs amber openmm --outdir jobs/ethanol_conversion
```

从 GFN2 Hessian 重推 bond/angle（xtb 在 `<outdir>/xtb/` 里运行）：

```bash
PYTHONPATH=src python -m ffcompiler.cli derive-bonded \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --hessian-provider gfn2 --method mseminario \
  --output gromacs amber openmm --outdir jobs/ethanol_mseminario
```

- `--hessian-provider {gfnff,gfn2,fchk,orca-hess,cp2k-out}`（`fchk`、`orca-hess`、`cp2k-out`
  需要 `--hessian-file`）
- `--method {mseminario,m2seminario,drih}`；`--method localhess` 选上面的实验方法
- `--no-symmetrize` 关闭对称等价项平均
- `--linear-angles {error,library}`：默认 `error`；`library` 保留 >175° 角的原参数

从 MACE-POLAR 松弛扫描拟合可转键二面角（缓存在 `<outdir>/scans/`）：

```bash
PYTHONPATH=src python -m ffcompiler.cli derive-torsions \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --output openmm --outdir jobs/ethanol_torsions
```

- `--step-deg`（默认 10°）、`--device`（默认 cpu）、`--scan-cache`（默认 `<outdir>/scans`）
- `--refinement-rounds {0,1,2,3}`：单遍之后的追加顺序轮；默认 0 保持单遍

检查导出文件在各引擎里复现 IR 能量：

```bash
PYTHONPATH=src python -m ffcompiler.cli check-equivalence \
  --prmtop ... --inpcrd ... --sdf ... --workdir jobs/ethanol_check
```

在一份 IR 上串联各步——每步写出 `ir.json` + `coords.json`，下一步读入两者：

```bash
# 1) 电荷
PYTHONPATH=src python -m ffcompiler.cli assign-charges --prmtop P --inpcrd C --sdf S \
  --method am1bcc --output openmm --outdir jobs/mol/q
# 2) 从 Hessian 推 bond/angle，接第 1 步
PYTHONPATH=src python -m ffcompiler.cli derive-bonded \
  --ir-in jobs/mol/q/ir.json --coords-in jobs/mol/q/coords.json \
  --hessian-provider gfn2 --method mseminario --output openmm --outdir jobs/mol/b
# 3) 导出最终 IR + 坐标
PYTHONPATH=src python -m ffcompiler.cli convert \
  --ir-in jobs/mol/b/ir.json --coords-in jobs/mol/b/coords.json \
  --output gromacs amber openmm --outdir jobs/mol/final
# 拟合记录过的输入事后被改过则退出码 1
PYTHONPATH=src python -m ffcompiler.cli check-staleness --ir-in jobs/mol/b/ir.json
```

`global-fit` 暴露 QMDFF 式全局 Hessian 拟合（Hessian 选项相同；`--prior`、`--alpha`）；
任一等价类撞到 k = 0 下界时写出 `ir_failed_lower_bound.json`、不导出、退出码 1。

从 Python：

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
- **谐振项** 存成 E = ½ k Δq²，即 k 就是 Hessian 对角元。GROMACS 与 OpenMM 直接用 k；
  AMBER 用 K = k/2（E = K Δq²），另做 kcal 与 Å 换算。
- **DRIH** 内部用原子单位（Bohr、E_h/a₀²）；拉伸对角元按 E_h/a₀² ≈ 9.37583 × 10⁵
  kJ mol⁻¹ nm⁻² 换算，角对角元按 E_h/rad² = 2625.50 kJ mol⁻¹ rad⁻²（不另乘 a₀²）。
- 支持的 Hamiltonian：气相、无截断、Lorentz–Berthelot 混合、统一 AMBER 1-4 缩放。
  其他在导入时拒绝。

## 验证状态

四道 Gate 已于 2026-10-18（可选顺序重拟合接入之后）在当前源码上重跑并全部通过，数值与
2026-10-13 那轮一致；每个 Gate 的 pytest 都是 524 + 59 通过。Gate 判据是默认单遍；顺序模式
由单元测试与逐位重放核对覆盖（`experiments/20261018_refinement_verify/`）。

| Gate | 运行 | 结果 |
|---|---|---|
| G0 — 导出等价 | `jobs/gate_mvp0_20261018_2/` | passed：11 分子 × 4 几何、8 条比较腿、2596 行、0 失败 |
| G1 — Hessian 方法 | `jobs/gate_mvp1_20261018_2/` | passed：三种方法在 12/12 参考用例上与 Sobtop 一致；全局拟合自洽 6.0 × 10⁻⁹ |
| G2 — MLP 二面角 | `jobs/gate_mvp2_20261018_2/` | passed：12/12 可转键全部由 MLP 拟合、0 兜底；非 sander 跨引擎腿全部通过。已知限制（决策 0017）：sander 在 π 附近咬合拟合相位（nma dihedral 2.3 × 10⁻⁴ > 2 × 10⁻⁴ kcal/mol），只报告 |
| G3 — 电荷 provider | `jobs/gate_mvp3_20261018_2/` | passed：AM1-BCC 复现全部 12 套固件电荷（最大 4.3 × 10⁻¹¹ e）；3 provider × 11 分子总电荷守恒（≤ 8.3 × 10⁻¹⁷ e），跨引擎比较全部通过；换电荷后 MLP 二面角重拟合：电荷不变时逐位相同，NAGL 电荷下 12/12 键保持 MLP |

导出力场对参考 Hessian 的频率 RMSE（11 分子中位数，cm⁻¹；只报告，不进 Gate）：

| Hessian | GAFF2 | mSeminario | m2Seminario | DRIH | 全局拟合 |
|---|---|---|---|---|---|
| GFN2 | 110 | 109 | 188 | 219 | 118 |
| GFN-FF | 101 | 95 | 76 | 203 | **57** |

自己跑 Gate（每个写新的 `jobs/gate_*_<date>[_N]/` 目录；旧目录保留）：

```bash
python scripts/gate_mvp0.py
python scripts/gate_mvp1.py --jobs 6   # 6 并行 xtb × 4 线程
python scripts/gate_mvp2.py --jobs 8 --threads 4 [--resume 旧G2目录]   # MACE 扫描
python scripts/gate_mvp3.py --scan-cache jobs/<该G2运行>/scans          # 必填
PYTHONPATH=src python -m pytest -m "not external" -q   # 单元测试（不需要外部程序）
```

## 已知限制

- **已知缺陷（MVP-2）**：苯甲酸甲酯 1-2 键的 MLP 拟合二面角对 MACE 剖面的复现比 GAFF2 差
  （剖面 RMSE 3.36 对 2.70 kJ/mol；拟合 RMSE 0.60）。12 根可转键里只有这一根如此。2026-10-17
  诊断（`jobs/torsion_diagnosis_20261017/report.md`）：键是逐根对着"其余键为 GAFF2"拟合的，
  1-2 项吸收了硬的 GAFF2 甲基（0-1）二面角；0-1 被替换后这份补偿就错了。`--refinement-rounds 3`
  降到 0.44 kJ/mol、合并残差 0.83 → 0.22（样本内；`jobs/torsion_b3_sequential_20261017/report.md`）；
  所有更新顺序收敛到同一不动点，第 3 轮后顺序差 ≤ 1e-4 kJ/mol
  （`jobs/torsion_b3_order_20261018/report.md`）。保持 opt-in：默认仍是单遍，留出验证（B4）
  未做，顺序模式不进 Gate。
- **扫描质量（MVP-2）**：nma 1-3（酰胺 C–N）正反扫滞后最大 57 kJ/mol（14/36 帧超过 1 kJ/mol，
  集中在 ±60–120° 能垒区）。这是扫描缺陷而非拟合缺陷；是否重扫未定，未设自动兜底阈值。
- 一根键的二面角能量放在单一四元组上（MVP-2 计划 §7 第 3 条）；适用范围未定。
- ORCA `.hess` 读取器只对一个闭壳层 ORCA 6.1.1 样例核对过（频率相对差 5.7 × 10⁻⁸）。
- 线性角（> 175°）不推导；乙腈只作诊断用例。
- m2Seminario 在 GFN-FF Hessian 下可能给出负的环角力常数（环丙烷）；Sobtop 给出同样的值。
  `derive_bonded` 抛 `NonPositiveForceConstantError` 而不是静默兜底——这类体系请用 mSeminario
  或 DRIH。
- mSeminario 的键常数跟随书写方向（与 Sobtop 相同）；对称化（默认开）消除方向不对称。
- 导出写单残基。IR 携带可选周期盒（`box_nm`，决策 0023）：NumPy 能量用最近镜像 + 键图
  解包裹，GROMACS/OpenMM 导出完整三斜盒，AMBER prmtop 仅正交；周期跨引擎腿尚未进入
  check-equivalence。

## 仓库结构

```text
src/ffcompiler/
  ir/          ForceFieldIR schema、校验、NumPy 能量、JSON I/O、坐标 sidecar
  perception/  几何 bond/angle/dihedral/improper 感知、Cordero 半径、PBC 镜像
  fflib/       Sobtop 格式 LJ/bonded/bondcrit 参数库读取与通配匹配
  importers/   AMBER prmtop 导入器
  exporters/   GROMACS、AMBER、OpenMM 写出器
  formats/     原生 prmtop / rst7 读写
  hessian/     CartesianHessian、fchk / ORCA .hess / CP2K .out / xtb provider、fallback 阶梯
  bonded/      内坐标、mSeminario、m2Seminario、DRIH、localhess、对称性、全局拟合、staleness
  torsion/     可转键、MACE-POLAR 松弛扫描、MM 扣除、Fourier 拟合、MLP provider
  charge/      电荷 provider（载入、AM1-BCC、NAGL）、learned-BCC 骨架、写回与二面角重拟合
  validation/  跨引擎比较、GROMACS 适配、振动、pytest 记账
scripts/       Gate、固件生成、Sobtop 参考数据
tests/         单元与跨引擎测试、固件
docs/          计划、决策记录、理论
jobs/          正式运行（只跟踪 manifest/acceptance/report）
```

## 文档

- [forcefield_compiler_design.md](forcefield_compiler_design.md)——设计稿（验收基准；下一阶段
  各节 §3.1、§6.4.4–6.4.6、§10.4–10.5 标注"已设计未实现"）
- [docs/history/ACCEPTANCE_20261018_CLOSURE.md](docs/history/ACCEPTANCE_20261018_CLOSURE.md)——
  当前版本的结案验收（34 条设计要点；§2 为已知缺陷、待拍板项与未做验证）
- [docs/theory/hessian_methods.md](docs/theory/hessian_methods.md)——推导、单位、约定
- [docs/decisions/](docs/decisions/README.md)——决策记录 0001–0024
- [docs/superpowers/plans/](docs/superpowers/plans/)——路线图与 MVP 计划
- [HANDOVER.md](HANDOVER.md)——当前已验证状态与下一步
- [CHANGELOG.md](CHANGELOG.md)

## 路线图

**已完成并有 Gate**：MVP-0（IR + 导出器）、MVP-1/1b（Hessian 方法 + 全局拟合）、MVP-2
（MLP 二面角）、MVP-3（电荷 provider）——G0–G3 在当前源码上通过。

**只有框架，研究未开始/未完成**：MVP-4 learned BCC（骨架 + 立项草案，训练未批准）、
MVP-5 fallback 调度（阶梯已就位，无不确定度调度）。

**已设计、未实现**（2026-10-25 设计增补，决策 0024——设计稿 §3.1、§6.4.4–6.4.6、§10.4–10.5）：

- **MVP-6 分子入口**——SDF 进、骨架 IR 出：antechamber 作外部 GAFF2 typer（只做 typing）、
  GAFF2 库值作先验/兜底、不静默生成构象；解锁单命令 CLI 与 preset（§17/§22）。
- **MVP-7 扫描质量**——多起点松弛扫描 + 逐格点分支诊断；分支只进诊断、不进拟合；nma 1-3
  重扫是第一个应用。
- **MVP-8 平衡几何联合拟合**——q⁰ 进拟合变量 + 梯度匹配（§6.4.3-1 预留接口），opt-in，
  对现行全局拟合逐位退化。
- **MVP-9 多构象联合拟合 + 交替收口**——同一目标函数 C ≥ 1、共享 q⁰、留出构象验证，以及
  charge → torsions → global-fit 交替的 ε 容差完成判据。
- TorsionDrive 保留给二维耦合扫描的 opt-in；是否装进生产环境是待定的显式决策。

**待做验证**：B4 二面角重拟合留出验证、真实 MACE 参考的正式外部组合验收（C5）、Level 3–4
benchmark（设计 §20）。

**待拍板**：默认 bonded 方法、默认重拟合轮数、版本号、TorsionDrive 安装——见 HANDOVER.md
与结案报告 §2。
