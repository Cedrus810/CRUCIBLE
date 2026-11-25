# MVP-6 分子入口与单命令 CLI 施工计划（2026-11-11）

状态：**已执行（2026-11-11）**。`pytest -m "not external"` 629 passed；验收要点 2（12 固件
SDF 重建等价）外部测试 12/12 通过（本机）；要点 1 的完整 Gate（G4）待 Gate 机器资源。
台账：`.superpowers/sdd/2026-11-11-mvp6-molecule-entry/progress.md`。
依据：设计稿 §3.1（入口合同）、
§17（preset）、§21（模块归属）、§22（CLI 草案）、§23（MVP-6 交付物）、决策 0024。
本轮同批补做 A4 快修清单（2026-10-11 计划 §3 的可单元测试项）。**MVP-7/8/9 不在本轮**
（各自是 MVP 级工程，设计已就绪，见 HANDOVER 待办）。

## 1. 开工时二选一的裁定（§3.1 授权，本计划锁定）

**antechamber 用 `-c bcc -at gaff2`**（typing + AM1-BCC 一次取得）。理由：本机
AmberTools 的 antechamber 实测无 `-c zero` 电荷方法（2026-11-11 探测，`Unknown charge
method (zero)`）；§3.1 预留的正是这条退路，且与 12 固件的生成配方完全一致——验收要点 2
的"由 SDF 重建与固件逐项等价"因此最强。后果与边界：

- 骨架电荷 = antechamber AM1-BCC，provenance 记 am1bcc（决策 0024 裁定 1 原文：仍可被
  charge provider 覆盖）。**"占位值/必须覆盖"语义属于 `-c zero` 路径，本轮不触发**。
- typing 结果与电荷方法无关（GAFF2 typing 只看连接关系），验收要点 2 的等价性不受影响。

配方其余部分照搬 `scripts/make_fixtures.py`（已验证 12 固件）：antechamber → 电荷均摊
（Σq = -nc，1e-4 门）→ parmchk2 `-s gaff2` → tleap（`source leaprc.gaff2` / loadmol2 /
loadamberparams / saveamberparm）→ 现有 `importers/amber_prmtop.from_amber`。
**make_fixtures 的 `equalize_mol2_charges`/`ambertools_version` 移入
`typing/gaff2.py` 库函数，make_fixtures 改为导入（数学逐行不变，固件不重新生成）。**

## 2. 几何合同（§3.1，逐条实现于 `molecule/io.py`）

只检查、不优化、不生成：RDKit 读 SDF（`removeHs=False`）；(a) 解析失败→干净报错；
(b) 必须有 3D 构象（无构象或全零坐标→报错，**不做** ETKDG 兜底）；
(c) formal charge 总和 = 声明净电荷（未声明时取 SDF 总和并回显；矛盾→报错）；
(d) 严重原子穿距检查：任意两原子距离 < 0.5 Å→报错（常量进 provenance）；
(e) 记 SDF sha256。入口不写任何坐标文件。

## 3. CLI 形式（§22 逐字 + 一个 entry 子命令）

1. **单命令（§22 原文形式）**：`ffcompiler ligand.sdf --preset fast --output gromacs amber openmm`。
   实现：顶层 argparse 加 `sdf`（nargs=`?`）+ `--preset/--charge/--bonded-method/--hessian/
   --hessian-file/--torsion/--device/--step-deg/--refinement-rounds/--net-charge/--workdir/
   --force/--output`；subparsers 改为非必需，`command is None and sdf` 时走单命令编排。
2. **entry 子命令**（分段使用的入口）：`ffcompiler enter ligand.sdf --outdir D
   [--net-charge N] [--name X] [--output gromacs ...]`——只走配方 + `from_amber` +
   provenance 增强，写 ir.json/coords.json（C1 合同）与配方日志，可选引擎文件。
3. **编排（单命令）**：entry → derive-bonded（preset 的 hessian/method）→
   assign-charges（见 4）→ derive-torsions（preset 有 torsion: mlp 时；--reference mace）→
   最终 IR 一次性导出 `--output` 的引擎文件。每步一个子目录
   `<workdir>/{entry,bonded,charge,torsion}`，中间步 `--output []` 不重复写引擎；
   复用现有 `_do_*` handler（构造 Namespace），保证单命令与分段命令走同一代码路径。
   `workdir` 默认 `./<sdf stem>_ffcompiler`，已存在且无 `--force` 报错；
   写 `run_manifest.json`（各步命令、版本、argv、provenance 摘要）。

## 4. preset 表（§17 逐字冻结于 `presets.py`；不可用项显式报错）

`PRESETS` 字典逐字照抄 §17.1–17.3 + §17.4 的字段。语义裁定（全部写入 run_manifest）：

- `charge: new_charge`（MVP-4，未立项）：**报错**，提示 `--charge entry|am1bcc|nagl|loaded`。
  `--charge entry` = 保留骨架 AM1-BCC 电荷并在 provenance 记 `charge_kept: entry_am1bcc`
  （§3.1 退路语义，非静默沿用）；`am1bcc/nagl/loaded` 走现有 assign-charges。
- `hessian: orca`（high preset）：CLI 的 orca provider 读 `.hess` 文件，preset 不自带
  DFT 计算——缺 `--hessian-file` 时报错并说明。
- **§17.4 rigid 的 `rigid_torsion: method: drih` 是设计表内部不一致**（DRIH 是 bonded
  方法，非二面角 provider）：v1 实现 rigid = typing + bonded(gfn2/m2seminario) +
  charge 覆盖 + **不做 MLP torsion 步**；`rigid_torsion` 字段原样存入 manifest 并标注
  "设计表待澄清，未执行"。此偏离在 HANDOVER 记为待拍板项。

## 5. provenance（§3.1 要求）

`ir.provenance["entry"]`：recipe 字符串、ambertools 版本（conda-meta）、antechamber/
parmchk2/tleap 解析路径与完整命令行、leaprc 名、SDF sha256、几何检查门（0.5 Å）、
电荷均摊日志摘要。各项 `source` 保持 from_amber 的 `gaff2-library` 不变；
`check-staleness` 语义不动（entry 不新增拟合依赖）。

## 6. 测试

- **unit**：SDF 检查（解析失败/缺 3D/净电荷矛盾/原子穿距）；equalize 数学；
  preset 表逐字冻结与编排映射（fast/balanced/high/rigid 的步骤旗标、new_charge 报错、
  rigid 无 torsion 步）；CLI 解析（§22 形式、enter 形式、裸调用退出码 2）；
  provenance 增强后 source 字段不变。
- **external（本机可跑，无 xtb/MACE 依赖）**：12 个固件分子各从其原始 SDF 走 entry，
  与既有 prmtop 导入结果对照——验收要点 2：拓扑（连接/元素/类型）、LJ（σ、ε）、
  exclusions/1-4 逐项一致；bond/angle/torsion 库值逐项一致；电荷 allclose(1e-6)
  （同一配方同一 AmberTools，期望逐位同，不作硬判据）。antechamber 全程 ~1–3 min。
- 单命令端到端（entry→…→export）需要 xtb+MACE，**本机没有**：不在本轮添加不可验证的
  e2e external 测试；编排逻辑以 unit 的映射测试 + `enter` 的 external 测试覆盖，
  端到端留给下次 Gate 资源（HANDOVER 记为待办）。

## 7. A4 快修（同批，2026-10-11 计划 §3 表的可单元测试项）

| 条目 | 处置 |
|---|---|
| 畸形 mol2 报错（charge/loaded.py） | 捕 IndexError/ValueError → 报"mol2 第 N 行列数不足/非数字"，附单元测试 |
| 提前检查 step_deg | **实施修正**：检查放 CLI 层（derive-torsions 与 run 两条入口，扫描开始前拦 `360/step ≥ 13`）；`relaxed_scan` 本身保持纯净——现有测试用 120° 粗网格驱动假计算器，scan 层加检查会破坏它 |
| TorsionRefitError 提示缺参数名（torsion/mlp.py:240） | 消息加 `--scan-cache` |
| rotors 列表共享（torsion/mlp.py:224,269） | provenance 侧存副本 |
| ASE/模型异常分类（torsion/scan.py:85-94） | 包装消息带异常类名，分类逻辑不变 |
| 重复 find_rotors（scripts/gate_mvp3.py:205） | 只算一次 |
| G1 `--providers` help（scripts/gate_mvp1.py:537） | 补 0019 说明 |
| 可转键对称类报告 / SDF 解析提示 / 主进程线程上限时机 / 甲铵正离子回归 | **本轮不做**（前两条例由见表；线程时机动 Gate 语义；MACE 回归需 MACE 机器且不可在本机验证，不加入未经运行的 external 测试） |

## 8. 验收对照（完工后逐条自检）

§3.1 验收要点：1（新分子端到端）→ 本轮以 enter + 编排映射覆盖，全链路待 Gate 资源；
2（12 固件等价）→ external 测试本机执行；3（provenance/staleness 语义）→ 单测断言
source 不变、entry 不新增拟合依赖。Gate：MVP-6 尚无 G4 脚本，**本轮不建 Gate**
（需 xtb/MACE 机器，资源另报——列入 HANDOVER 待办）。

## 9. 明确不做

不自研 typing；不改 importers/amber_prmtop.py 与任何现有 IR 语义；不重新生成固件；
不实现 MVP-7/8/9、MVP-4 训练、MVP-5 调度、RESP、TorsionDrive；不动版本号。
