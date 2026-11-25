# CRUCIBLE（原 forcefield_compiler）交接状态 — 2026-10-28（Asia/Tokyo）

## 2026-11-24 第一批代码层续：xtb_sp 诊断阶段、多起点诊断报告脚本

用户 2026-11-24 确认"就是 22 给你干了一半"，接着补 11-22 批次剩下的代码。计划仍是
`docs/superpowers/plans/2026-11-20-coverage-validation-batch1-plan.md`，§7 待拍板 1–6 仍待定。
- **xtb_sp（计划 §2.2 第 4 项）**：`hessian/xtb.py` 新增 `run_xtb_singlepoint`（`--gfn 2 --grad
  --chrg Q --uhf 0`，入口几何、不优化；运行前清掉旧 gradient/energy；开壳层沿用同一电子数预检）
  与 `read_xtb_gradient`（Turbomole gradient 文件，取最后一个 cycle，Eh 与 Eh/Bohr 不换算）。
  格式按本机 xtb 6.7.1 真实输出核对后写的。`run_xtb_hessian` 未动。驱动
  `scripts/coverage_diagnose.py` 的 xtb_sp 阶段改为真实实现（opt-in，`--xtb-threads` 默认 4）；
  失败按既有 verdict 规则记缺项、不阻塞，规则本身未改。
- **多起点诊断报告（计划 §3.1）**：新脚本 `scripts/torsion_multistart_report.py`，只读多起点
  derive-torsions 的 `torsion_report.json`（可选扫描缓存补逐点角度/x0 能量/逐腿滞后、可选单起点
  基线对照 fit_rmse），产出 report.md / per_point.csv / summary.json。只列数、不判定、不设阈值。
- **验证（本机）**：`test_hessian_xtb.py` + `test_coverage.py` 非 external 30 passed；新 external
  `test_run_singlepoint_ethanol_gfn2` 1 passed（5 s，真调 xtb）；`test_torsion_multistart_report.py`
  4 passed。冒烟 `experiments/20261124_multistart_report_smoke/`（ethanol，`--reference mm
  --step-deg 20 --multistart 2`，1 min 47 s）：报告走通；单起点基线用同一缓存 1.9 s 完成，
  前后 `scans/` 四个 npz 修改时间不变——默认 {x0} 缓存键与多起点 x0 腿一致，实跑确认。
  全量回归仍未跑（11-22 用户指令）。
- **§2.2 第 3、6 项（用户 2026-11-24 "你可以接手给干了"）**：驱动新增两个 opt-in 阶段——
  `fallback`（xtb GFN-FF --ohess → mSeminario，线性角取 library 以便计数；兜底 = 线性角保留
  GAFF2 库值，同时记"fast preset linear_angles=error 会拒绝"；非正 k 记 fail）与 `scan_probe`
  （find_rotors 第一根键，现行 MLPTorsionProvider + MACE、10°，计时；缓存在
  `<name>/scan_probe/scans`；无可转键记 skip；该键兜底记缺项）。新参数 `--mace-threads`、`--device`。
  verdict 规则文字未改，只把 scan_probe 加入判定枚举（属"后续拟合步兜底 → 缺项"）。
- **修了 11-22 驱动的三个潜伏错误**（此前驱动从未在真实分子上跑过，单测没覆盖入口）：
  ① `stage_entry` 返回 ir.json/coords.json 路径但不写文件；② 调用了未定义的 `_write_record`；
  ③ 下游坐标用 `read_rst7`（Å）当 nm，几何放大 10 倍（xtb SCC 不收敛、GFN-FF Hessian 平移残差
  6.9e-3 被质量门拒）——改为读入口落盘的 coords sidecar（nm）。新 external 测试
  `test_driver_light_stages_on_real_ethanol` 覆盖这三处（含 coords×10 与 inpcrd 逐位一致）。
- **冒烟** `experiments/20261124_coverage_stages_smoke/`（ethanol/acetonitrile/benzene，九个阶段全点名，
  16 线程，5 min 24 s；前两轮失败的产物改名保留 `out_2_rst7_angstrom_bug/`、`run_1/2_*.log`）：
  ethanol ok（scan_probe 0-1 键 250 s，fit_rmse 0.022）；benzene ok（无可转键 skip）；
  acetonitrile blocked——gmx 腿 angle 1.1e-3 kJ/mol 超容差 1.0e-4（相对 2.4e-3，近线性角，与决策
  0013 记录的混合精度现象相同），fallback 记 1 个线性角兜底。
- **用户决定 2026-11-24**（原话："可以  就这么定"，针对 Claude 提议"混合精度 GROMACS 腿的超差只记作
  精度缺项、不再判 blocked，其余各腿超差仍算 blocked"）：已实现——超差行全部来自 gmx 腿且
  precision=mixed、无其他腿超差/出错/力超差时，export 记 ok + `gmx_mixed_precision` 缺项（超差行原样
  留在记录里）；verdict 规则 docstring 同步。单测 6 例（1 例转缺项 + 5 例仍 blocked）。实跑：
  acetonitrile 重跑（`out_4_precision_rule/`）判 ok_with_gaps，缺项 = gmx 混合精度（angle,total；
  rel 2.4e-3）+ 1 个线性角库值。`test_coverage.py` 31 passed。
- **覆盖集清单（计划 §2.1，用户已认可 34 个）**：新脚本 `scripts/make_coverage_set.py` →
  `benchmarks/coverage_set/molecules.json` + `sdf/`（22 个新增）。SDF 来源用户未指定，沿用固件配方
  （复用 `make_fixtures.build_sdf`：RDKit ETKDGv3 种子 20261004 → MMFF94）；12 个固件直接引用
  `tests/fixtures/gaff2/` 原文件（计划 §2.4 不动 fixtures）。`tests/unit/test_coverage_set.py` 3 passed。
- **34 分子轻阶段冒烟** `experiments/20261124_coverage_light_smoke/`（entry/rotors/frcmod/export/refcheck，
  47 min 38 s——当时 NFS 慢，进程长时间处于 D 状态 rpc_wait，不是计算量）：0 blocked、31 ok、3
  ok_with_gaps（acetonitrile gmx 混合精度；acetate、methylammonium 离子 MACE 适用性待定）；frcmod 全部
  0 ATTN；19 个分子有可转键、15 个没有。重阶段（mace_smoke/xtb_sp/fallback/scan_probe）未跑，正式跑批
  待用户批准资源。
- **覆盖集正式诊断（用户 2026-11-25 "183走起"）**：183（node-C，40 核）上 8 worker × 5 线程（= 40 物理核；
  计划表原写 181、4×9，用户改到 183 后按空闲机器改为 8×5），34 分子 × 9 阶段，输出
  `jobs/coverage_20261124/`，启动脚本 `jobs/launch_coverage_20261124.sh`，PID 在该目录 `pids.txt`
  （worker0–7：79515–79543，02:48:58 起）。183 上的运行账号无个人工具配置 → 用 `FFC_CONFIG` 指向同目录
  `tools.toml`（原作者账号的配置原样复制；含真实路径，被 jobs/** 规则挡在库外）；183 上 check-env 全 ok。
  续跑 = 重跑启动脚本。已全部结束并汇总（见下条）。
  误启动一次（02:48:24）：脚本把数组命名为 `GROUPS`（bash 内置只读，赋值被静默忽略），11 个进程拿到
  组 ID 当 `--molecules`，参数检查报错自行退出、无产物；日志存 `misfire_GROUPS_builtin/`，已改名
  `MOL_GROUPS` 并在 183 上空跑核对后重启。
- **覆盖集正式诊断结果**（183，02:48:58–03:32:42，44 min；汇总在 183 上执行，8 个 worker 日志 0 个
  UNEXPECTED/Traceback）：34 分子 **26 ok / 8 ok_with_gaps / 0 blocked**（`jobs/coverage_20261124/report.md`）。
  - export、mace_smoke、xtb_sp：34/34 ok；scan_probe：19 ok（全部 status=mlp，无兜底）、15 skip（无可转键）。
  - **fallback（fast preset bonded：GFN-FF --ohess → mSeminario）5 个失败** → 单命令 `run --preset fast`
    对它们会失败：acetone（2 虚频，最低 −56 cm⁻¹）、methyl_acetate（−63）、styrene（−41）、
    glycine_zwitterion（−222）——GFN-FF vtight 优化落在鞍点、被 Hessian 质量门拒；oxirane——mSeminario
    两个三元环角 k ≈ −92（与决策 0014 环丙烷同类方法局限）。
  - **glycine_zwitterion 的 scan_probe 病态但判定未报**：fit_rmse 1473 kJ/mol、滞后 38.8 kJ/mol，status
    仍为 mlp（扫描质量不设自动阈值，§10.4 / 用户 10-17 裁定）。【未核实推断】气相两性离子在约束扫描中
    发生质子转移。只报告，未加门。
  - 扫描滞后偏大：aniline 17.5、anisole 6.35 kJ/mol（其余 ≤ 0.06）——与 nma 1-3 同类扫描分支问题，
    可作 MVP-7 多起点重扫候选。
  - scan_probe 耗时 140–1742 s/键（5 线程），19 键合计 12180 s；大致随原子数增长。
  - **4 个虚频分子已诊断（用户 2026-11-25 "可以"）** `experiments/20261125_saddle_diagnosis/`：从 xtb 自己写的
    `xtbhess.xyz`（沿虚频扭曲）用同一 `run_xtb_hessian`（GFN-FF --ohess vtight）重启一次，四个全部无虚频、
    过 `check_hessian`、能量更低（−0.18 ~ −9.4 kJ/mol）→ GFN-FF 优化停在鞍点、管线不重启（代码里不处理
    xtbhess.xyz），不是势能面或 Hessian 问题。acetone/methyl_acetate 虚频主要是甲基扭转（重叠 0.71/0.97），
    styrene/glycine 为混合模式。停在鞍点的原因（起点近对称）是推断，未验证。
  - **已修（决策 0027，用户 2026-11-25 "可以"）**：`run_xtb_hessian` 在最低频率 < −20 cm⁻¹ 且有
    xtbhess.xyz 时从它重启（`imag_restart_<k>/`，最多 2 次），只在发生时记 meta/provenance；不触发时逐位
    不变（ethanol 实测 hessian 逐字节相同）。验证：`test_hessian_xtb.py` 23 passed（新 5 单测 + acetone
    external）；test_assign/pipeline_chain/fallback/coverage 64 passed；5 个失败分子重跑 fallback
    （`experiments/20261125_imag_restart_verify/`）：4 个转通过，oxirane 仍失败（负环角 k，另一原因）。
    G0–G3 未重跑（Gate 分子都过质量门，不触发重启）。
- **覆盖集正式结果（第 2 轮，0027 之后；用户 2026-11-25 "181继续跑"）**：181 上 8 worker × 4 线程
  （181 有他人负载，未用满 40 核），`--force` 重算，08:57–09:06；扫描走缓存（缓存文件仍是第 1 轮 03:08）。
  **30 ok / 4 ok_with_gaps / 0 blocked**（`jobs/coverage_20261124/report.md` 已替换）：fallback 33 ok / 1 fail
  （oxirane 负环角 k）；剩余缺项 = acetate、methylammonium（离子 MACE 待拍板 2）、acetonitrile（gmx 混合
  精度 + 线性角库值）、oxirane。**glycine_zwitterion 现判 ok，但其 scan_probe 仍病态（fit_rmse 1473 kJ/mol，
  判定不看扫描质量）**。修复前结果留底在 `jobs/coverage_20261124/prefix_0027/`（34 record + report/summary/manifest）；
  启动脚本加 `NT`/`TAG` 与参数透传，第 2 轮 PID/日志带 `_r2`。
  - 另见：NonPositiveForceConstantError 的提示文字固定写"可改用 mseminario 或 drih"，当前方法就是
    mseminario 时会误导（老代码，未改）。
- **离子/两性离子参考势改用 xtb GFN2（决策 0028，用户 2026-11-25 "可以"，计划 §7 待拍板 2）**：新模块
  `torsion/xtb_ref.py`（判定 + ASE 计算器，每步调 `xtb --gfn 2 --grad`）；CLI `--reference` 加 `auto`（默认）
  与 `xtb-gfn2`，`run` 用 auto；覆盖集 refcheck/scan_probe 同规则。Gate G2/G3 直接建 MACE，不受影响。
  顺带修 11-22 潜伏错误：`run` 给 derive-torsions 的 Namespace 缺 `multistart`/`multistart_threshold`
  （run 走到二面角步必 AttributeError）；新测试从 handler 源码扫 `args.X` 逐项核对，确认能抓到原错误。
  测试：test_torsion_xtb_ref（含 external：力 = −dE/dx 有限差分、能量 = 直接单点）+ test_coverage 38 passed；
  test_cli_entry 10 passed。冒烟 `experiments/20261125_gfn2_reference_smoke/`：CLI methylammonium auto→gfn2，
  35 s 扫完、fit_rmse 0.018；随后 assign-charges am1bcc 按缓存 refit 成功（cache_only + `_no_reference`，
  身份不符会兜底 no_cached_scan，实为 mlp）；覆盖集 acetate、methylammonium 用 GFN2 → **ok**（各 28 s）；
  **glycine_zwitterion 在 GFN2 下扫描同样病态**（fit_rmse 1205 kJ/mol、滞后 38.2）——不是 MACE 特有，
  【未核实推断】气相两性离子约束扫描中质子转移。正式覆盖集（jobs/coverage_20261124）的两个离子仍是 MACE
  结果（ok_with_gaps），未按 0028 重跑。
- **正式覆盖集第 3 次更新（0028 之后；用户 2026-11-25 "可以"）**：181 上只对 acetate、methylammonium、
  glycine_zwitterion `--force` 重跑 9 个阶段（日志 `jobs/coverage_20261124/rerun_0028.log`），再汇总全量。
  **现行结果：32 ok / 2 ok_with_gaps / 0 blocked**（剩余缺项：acetonitrile gmx 混合精度 + 线性角库值；
  oxirane 负环角 k）。三者 scan_probe 均为 xtb-gfn2：acetate 28 s / rmse 0.0035；methylammonium 27 s / 0.018；
  glycine_zwitterion 117 s / **rmse 1205、滞后 38.2（病态，判定不看扫描质量）**。重跑前的 3 份 record、
  它们的 MACE 扫描缓存（同名 0-1.npz 会被 GFN2 覆盖）与汇总留底在 `jobs/coverage_20261124/prefix_0028/`。
- **发现（未改，待用户定）**：CLI 帮助称 `--multistart N` 为"额外的 ETKDG 起点数"，实际 N 是
  嵌入的构象数，与 x0/彼此重复的被去重（冒烟中 N=2 只得到 1 条 etkdg 腿）。
- **未做**：§20 状态标注（改设计稿，计划写明批准后随 §C 提交，计划未批）；一切正式计算。
- **用户决定 2026-11-24**：计划 §7 待拍板 1 认可 34 个分子（原话："可以"）；SDF 来源未定。
- 提交草稿：`.commit_drafts/20261125_1113_coverage-run-saddle-restart-ion-gfn2.md`（24 号，2026-11-25 用户确认整理）。

## 2026-11-22 第一批代码层：覆盖集诊断、MVP-7 多起点、独立验证装配（2026-11-20 计划 §A/§B/§C）

用户 2026-11-22 指示"先补代码"后动工。只写代码与单测，未跑任何计算/MD/外部工具，
未动环境、无 git 操作。计划：`docs/superpowers/plans/2026-11-20-coverage-validation-batch1-plan.md`（状态行已更新）。
- **§A 覆盖集**：`validation/coverage.py`（记录 schema、frcmod ATTN 解析、verdict 判定
  【blocked = entry/typing 失败或导出不等价；ok_with_gaps = frcmod 缺项/参考势适用性/拟合兜底；
  规则施工前锁定】、清单加载与汇总）+ 驱动 `scripts/coverage_diagnose.py`（默认轻阶段
  entry/rotors/frcmod/export/refcheck；mace_smoke opt-in；**xtb_sp 本轮未实现，点名即报错**；
  `<name>/record.json` 存在即跳过，续跑同松弛检查惯例）。**覆盖集清单 JSON 尚未建**（34 分子
  草案待拍板 1），跑批需要 antechamber 等外部工具，属另批资源。
- **§B MVP-7（§10.4 落地；默认 {x0} 逐位不变——有单测核对：默认缓存键无 start 块、
  collect_candidates 不改主数组、默认路径 provenance 无 multistart 块）**：
  - `clustering.py`：Kabsch 对齐（无反射）+ 对齐重原子 RMSD + 能量升序贪心聚类（§B 分支
    与 §C1 盆地共用一套实现）；
  - `torsion/scan.py`：新增 `collect_candidates`（逐方向候选块）与 `grid_origin_rad`
    （非 x0 起点网格对齐 x0 的 φ0）；默认调用路径不变；
  - `torsion/multistart.py`：`etkdg_starts`（ETKDG 嵌入 → 参考势优化 → 去重；身份含
    seed/参数/RDKit 版本，进缓存键）+ `merge_start_scans`（逐格点合并出包络 + 分支诊断：
    逐点分支数/分支间能差/包络相对单起点增益；**分支只进诊断、永不进拟合目标**）；
  - `torsion/mlp.py`：`MLPTorsionProvider(multistart=...)`（fresh = 调优项；refit = 按
    provenance 存档身份只读缓存，不调 RDKit/参考势）；缓存键加 start 块（x0 = None 与
    现行键逐位一致，旧 MACE 缓存不失效）；候选块以 cand_* 落盘，旧缓存无此块仍可读；
    `refit_torsions` 从 provenance 恢复多起点身份；CLI `--multistart N
    --multistart-threshold nm`（derive-torsions 与 mm-self 参考势两路都接）；
  - 分支聚类阈值默认 0.05 nm（施工裁定，可参数覆盖；§7 待拍板 4 未涵盖此项，列为可调）。
- **§C**：`validation/holdout.py`（B4 落地：盆地整进整出划分 + holdout 上限
  max_holdout_frac=0.4；单盆地/放不下 → 留出为空并写 note，不硬造）；`validation/variants.py`
  （`substitute_terms`/`comparison_matrix`：bonded/torsion/charge/vdw 单替换 + full，
  每个变体过 validate）；`validation/liquid.py`（`pack_box` 立方格点 + 抖动装箱，确定性
  seed，格点放不下即报错不静默降密度；`to_liquid_system` = PME + LJ 截断 + 色散修正，
  exceptions 语义同导出器，**导出器 NoCutoff 路径不动（决策 0023 语义）**；`liquid_dynamics`
  = NPT Simulation（LangevinMiddle + MonteCarloBarostat，CPU 平台））。
- **验证（本机、逐文件单跑）**：test_coverage 16 passed；test_torsion_multistart 16 passed
  （142 s，GAFF2 MM 作参考势，不用 MACE）；test_holdout_variants + test_liquid 17 passed。
  **全量 `pytest -m "not external"` 未跑**（用户 2026-11-22 指令"回归停掉"，命令未执行）；
  external / G0–G3 / G4 均未动。开工前 G0–G3 第 3 轮（10-28）仍为最新通过记录。
- **未做**：一切计算（覆盖集跑批、nma 1-3 多起点重扫、二聚体/液体基准）；CLI 不加覆盖集
  子命令（脚本驱动）；xtb_sp 诊断阶段；多起点诊断报告脚本；设计 §20 状态标注；§7 待拍板
  1–6 原样待定。
- 提交草稿：`.commit_drafts/20261122_batch1-code.md`（23 号）。

## 2026-11-11 MVP-6 分子入口与单命令 CLI（决策 0024 裁定 8 的首个实施）+ A4 快修

按施工计划 `docs/superpowers/plans/2026-11-11-mvp6-molecule-entry.md` 实施（开工裁定
`-c bcc` 在其 §1 锁定：本机 AmberTools 无 `-c zero`，取 §3.1 预留退路，与 12 固件配方一致）。
- 新增：`molecule/io.py`（SDF 几何合同：解析、必须有 3D、formal charge=净电荷、0.5 Å 穿距门、
  sha256；不做任何构象生成）、`typing/gaff2.py`（antechamber → 电荷均摊 → parmchk2 → tleap →
  现有 from_amber；provenance 增 `entry` 块，各项 source 保持 gaff2-library）、`presets.py`
  （设计 §17 逐字冻结）；CLI 新增 `enter` 与 `run` 子命令，§22 字面形式
  （`ffcompiler ligand.sdf --preset fast ...`）经 argv 垫片支持；编排复用各分段 handler。
  `equalize_mol2_charges`/`ambertools_version` 自 make_fixtures 移入库（数学逐行不变，固件未重新生成）。
- preset 语义：new_charge（=MVP-4 未立项）显式报错；`--charge entry` 保留入口 AM1-BCC
  （决策 0024 裁定 1）；设计名 orca → CLI orca-hess 翻译；rigid 的 rigid_torsion:drih 是
  设计表内部不一致（DRIH 是 bonded 方法）→ v1 不执行、run_manifest 记待澄清（**待用户拍板**）。
- 验证：`pytest -m "not external"` **629 passed**（原 603 + 新 26；yayoi，79 s）；新外部测试
  **12 固件 SDF 重建等价 12/12（5.3 s，本机，只调 antechamber/parmchk2/tleap）**——设计 §3.1
  验收要点 2/3 完成；真机冒烟 `run ethanol.sdf --preset fast --charge entry` **完整走通并导出
  gromacs 三件套，exit 0（约 40 min，时间几乎全在 MACE 扫描）**——验收要点 1 以冒烟形式达成，
  细节见台账 `.superpowers/sdd/2026-11-11-mvp6-molecule-entry/progress.md`。
- A4 快修六项（畸形 mol2 报错、step_deg 提前检查【CLI 层，scan.py 保持纯净】、
  TorsionRefitError 消息、异常分类、gate_mvp3 重复 find_rotors、gate_mvp1 help）；
  rotors 列表共享复核为已被 deepcopy 覆盖；SDF 提示确认不触发。
- **未做**：G4 Gate 脚本与正式端到端验收（需 Gate 机器资源另报）；MVP-7/8/9；MVP-4 立项；
  MVP-5 调度；perception/fflib 接线；localhess oracle；Level 3–4（见"代码侧未做的事"清单）。
提交草稿：`.commit_drafts/20261111_mvp6-molecule-entry.md`。

## 2026-11-11 受约束松弛检查 阶段 1（2026-11-01 计划获用户批准后执行）

按 `docs/superpowers/plans/2026-11-01-relaxed-check-plan.md`（用户 11-11 批准阶段 1）与施工计划
`docs/superpowers/plans/2026-11-11-relaxed-check-phase1-build.md` 执行完毕。新脚本
`scripts/torsion_relaxed_check.py`，报告 `jobs/torsion_relaxed_check_20261111/report.md`。
只读 10-13 G2 缓存、只算 MM（无 MACE/xtb/重扫）；12 键 × {only, all}，P0/X/Y/P1 按锁定判据计算。

- **实施偏差**：dmso 的 bg 松弛要同时固定 4 个同轴二面角，ASE FixInternals 投影不收敛
  （maxstep/epsilon 调参均不可救）。脚本内实现 `MultiDihedral`（对全部约束同时解 Gauss-Newton
  投影），可行集合、fmax、LBFGS 不变，统一用于全部 12 键。验证：methyl_benzoate 三键 bg_relax
  与 11-01（node-A、FixInternals）最大差 4.6e-6 kJ/mol（跨机复现）。首版投影符号错在写帧前
  由插桩发现并修复；汇总层两处 bug（RMSE 双减参考能、numpy 列未减 E_ref）修正后从帧 CSV 重建。
- **验收**（施工计划 §5）：P0 对 B2 十二键最大偏差 2.7e-11；自检 24 行全过（约束漂移最差
  9.7e-13 rad）；机器 yayoi（8 核），并行相共约 75 s。
- **结果**：三条停止条件均未触发。Y > σ=0.41 在 only 基 5/11（ethanol 0-1/1-2、mb 1-2、
  mb 2-4、nma 3-4），其中 ethanol 两键与 mb 2-4 的固定几何指标近乎完美（0.02–0.04）——
  笔记 §3.2 的横向自由度机制存在于无已知缺陷的分子上；X 分类 0 假阳、真阳多于 P0/P1（两基）；
  但逐帧比值 0.95–1.00，X 只在量级层面分类、不逐帧预测；mb 1-2 的 all 基上 X 失明（跨扭转折补
  住在被剥掉的拟合项里）。nma 1-3 单列不计入。
- **下一步**：阶段 2（新增 5 分子、需算 MACE 扫描 1–2.5 h）**待用户批准**；本结果不进 pipeline、
  不加门（2026-11-01 计划 §6）。提交草稿见 `.commit_drafts/20261111_relaxed-check-phase1.md`。

## 2026-10-28 外部工具位置配置（决策 0026）

按用户要求（"给别人交付的时候应该给人家一个 config""允许用户缺项，缺了就少对应的功能"）实现了以下内容：
- `ffcompiler/toolconfig.py`：工具查找顺序为 `FFC_*` → 用户自己的 `~/.config/crucible/tools.toml`（不入库）→ PATH → 旧默认位置。
- `init-config`、`check-env` 两个子命令。
- `crucible.example.toml` 模板；`env.sh` 的路径可以覆盖。
- 缺 sander 时，CLI 跳过 sander 比对；Gate 仍要求 sander 必须在。

计划与 ledger：`docs/superpowers/plans/2026-10-28-tool-config.md`、`.superpowers/sdd/2026-10-28-tool-config/progress.md`。

验证：
- `pytest -m "not external"` 603 passed。
- 用非原作者账号、不改 HOME、不设 `FFC_*`，只靠配置文件：`check-env` 8/8 ok；原来失败的两个 external 测试（gfnff CLI、sander 周期 prmtop）都通过。

**只靠配置文件的全量验证（10-28 11:17–11:43，node-A，非原作者账号，不改 HOME、不设 `FFC_<工具>`）：全部通过。**
- 全量 external：60 passed（`experiments/20261028_tool_config_verify/`）。
- G0–G3 第 3 轮（`jobs/gate_mvp{0,1,2,3}_20261028_3/`）四个都 passed；pytest 每个都是 unit 603 + external 60，0 失败；四个 report.md 与第 2 轮逐字相同；G2 metrics.csv 相同。
提交草稿：`.commit_drafts/20261028_tool-config.md`（排在 `20261028_review-fixes.md` 之后执行）。

## 2026-10-28 审阅 10-24 代码并修复（决策 0025）

用户要求审阅，随后指示“该修全修”。修了 8 处，详见决策 0025：
- 周期盒只作模拟盒元数据，NumPy 参考能量不再取镜像。原来的语义和三引擎都不一致，实测 −8.45 vs 8.8e6 kJ/mol。
- AMBER 周期 prmtop 补齐 SOLVENT_POINTERS / ATOMS_PER_MOLECULE。修之前 sander 一读就报错。
- GROMACS 导出检查盒子是否为下三角形式。
- CLI `cp2k-out` 改用文件里的几何。原来用输入几何覆盖，几何检查等于没做。
- bondcrit 窗口改为替换语义，依据是反编译的 setbondcrit。
- 负 pn 只取绝对值，phase 不变。
- improper 匹配改为按 findffimp 的规则。
- MACE 默认路径测试去掉写死的 `/home/user`。

验证：`pytest -m "not external"` 579 passed，0 failed（104 s）。相关的 external 测试已跑：`test_ir_box.py`、`test_torsion_mace.py`、`test_cli.py` 共 34 passed，其中包括新增的 sander 读周期 prmtop。还复现了原来的问题：sander 在 `ntb=1` 下现在能正常运行。证据在 `experiments/20261028_review_pytest/`。
**G0–G3 第 2 轮（10-28 07:57–08:21，node-A，非原作者账号）：四个 Gate 全部通过。** 四个 Gate 的 report.md 与 10-18 第 2 轮逐字相同，G2 metrics.csv 也相同；pytest 四次都是 unit 579 + external 60，0 失败。这一轮之前修了 `test_cli_derive_bonded_ethanol_gfnff`：CLI 子进程现在会继承 HOME 和 FFC_XTB。运行目录 `jobs/gate_mvp{0,1,2,3}_20261028_2/`。

第 1 轮（03:03–03:29，已被第 2 轮取代，结果保留）：四个 Gate 都判为未通过，唯一失败项是环境问题。
- 各 Gate 自己的判据与 10-18 第 2 轮逐位一致：G0 的 report.md 逐字相同；G2 的 metrics.csv 逐字相同；G1、G2、G3 的 report.md 只差 pytest 判据那一行（G1.7、G2.5、G3.5）。
- pytest：unit 部分四次都是 579 passed；external 部分四次都是 59 passed、1 failed，失败的都是 `test_assign.py::test_cli_derive_bonded_ethanol_gfnff`。
- 失败原因：这个测试启动 CLI 子进程时只传 `PYTHONPATH` 和 `PATH`，不传 HOME。在非原作者账号下，`~` 被解析成该账号自己的家目录，xtb 找不到。
- 区分原因的检查（`experiments/20261028_review_pytest/gfnff_cli_home_check.log`）：同一条命令，不传 HOME 时 exit=1，补上 `HOME` points at the author's home 后 exit=0。0025 没有改动这条路径。
- 用户同意后改了这个测试，然后跑了第 2 轮（见上）。
- 运行目录 `jobs/gate_mvp{0,1,2,3}_20261028/`。
提交草稿：`.commit_drafts/20261028_review-fixes.md`。


## 2026-10-25 README 重写（仅文档）

用户确认大改 README：三份（EN/zh-CN/ja）全量重写为"工具是什么"的结构——管线图前置、
能力分节（一份 IR 三引擎 / Hessian bonded / MLP 二面角 / 电荷 provider / 感知层与参数库 /
provenance）、真实 quickstart、验证状态与已知限制原样保留、路线图重写（已完成有 Gate /
只有框架 / **已设计未实现** MVP-6～9 / 待验证 / 待拍板）。分工约定：README 是门面，
**当前状态的权威记录仍在本文件**；Roadmap 里未实现项一律标注"designed, not implemented"，
与设计稿（下一阶段）标注对应。提交草稿：`.commit_drafts/20261025_readme-overhaul.md`
（README.md 的 10-24 未提交改动随之并入该笔；16 号 update 草稿的 add 列表已把 README 移出）。

## 2026-10-25 设计稿增补：下一阶段四项（仅文档，无代码无计算）

按用户指示把四项下一阶段方向写入设计稿，**全部为新增小节、既有章节编号未动**：
§3.1 分子入口（antechamber 作外部 GAFF2 typer，骨架 IR，不静默生成构象）、
§6.4.4–§6.4.6（平衡几何联合拟合＝§6.4.3-1 预留接口落地，多构象联合拟合与其共享
目标函数、C=1 退化，charge→torsions→global-fit 交替协议与 ε 容差组合完成判据＝
§6.4.2/#14 遗留的收口）、§10.4–§10.5（一维多起点扫描原生做＋分支只进诊断；
TorsionDrive 仅二维 opt-in，装环境属用户显式决策）、§23 MVP-6～9、§26、§27。
关键裁定见**决策 0024**；验收对照时新节按“未实现（已设计）”归类，待用户确认后
并入现行基准。nma 1-3 多起点重扫、MVP-9 的留出验证分别对应待办 5 与 B4 的落地
形式。提交草稿：`.commit_drafts/20261025_next-phase-design.md`。

## 2026-10-18 本轮无计算交付

用户授权“自己干除了计算之外的事儿”，随后要求继续准备结案。已把顺序重拟合接成可选能力：
`derive-torsions --refinement-rounds {0,1,2,3}`，默认 0；初次拟合与换电荷重拟合共用同一 provider，
保存原 replaced_terms、原兜底集合、轮数与顺序、逐轮指标和固定 proper 依赖。
测试代码及交接已准备（施工方未运行测试）。（2026-10-18 接收方验证：相关 6 个测试文件 66 passed；非 external 全量 524 passed；3 轮与 B3 实验逐位相同；默认 0 轮与改动前单遍逐位相同。见 `experiments/20261018_refinement_verify/`。external 与 Gate 未跑。）
细节与接收方检查顺序：[无计算交接](docs/superpowers/plans/2026-10-18-sequential-refinement-handover.md)。
提交草稿：`.commit_drafts/20261018_optional-sequential-refinement.md`。
下文早间 G0–G3/561 项通过是本轮修改前的源码快照；本轮改动之后，已于 08:50–09:15 在新源码上重跑 G0–G3（第 2 轮），全部通过，见“正式运行”。

## 必读（开工前全部读完）

- `仓库代理规程`：设计稿路径（`forcefield_compiler_design.md`，现行，验收/开工前全文读完；`.v0.md` 已作废）、环境（openmm_dev 是生产环境）、分工。
- 本文件：当前状态、待办、未提交的草稿。
- `docs/superpowers/plans/2026-10-15-crucible-bc-detailed-next-steps.md`：B（二面角质量诊断）与 C（完整管线）的任务、依赖与完成标准——下一步工作的规格。
- `docs/superpowers/plans/2026-10-17-unbuilt-modules.md`（**先读 §0**）：10-17 搭建的 ORCA 读取器 / C1–C5 / MVP-5 框架 / MVP-4 骨架，及当晚验收修复清单与仍未做的事。
- **`docs/history/ACCEPTANCE_20261018_CLOSURE.md`：当前版本的结案验收（设计稿逐条对照 34 条），第 2 节是结案时的已知缺陷、待拍板项和没做的验证——下一阶段从这里开始。**（取代 A3 报告 §1 作为现行偏离清单；`ACCEPTANCE_20261013_A3.md` 保留作历史。）
- `docs/decisions/0017`–`0020`：仍在生效的验收与运行规则（OpenMM 为准、MVP-3 范围、子集不能通过、重拟合不升级兜底键）。
- `docs/decisions/0022`、`0023`、`0025`：10-24 吸收的 Sobtop 功能（感知层、fflib、localhess、CP2K、周期盒）以及 10-28 审阅对它们的修订。
- `docs/superpowers/plans/2026-10-08-ffc-mvp2-mlp-torsion.md` §7：MVP-2 还没定的问题（与 B 相关）。
- `docs/superpowers/plans/2026-10-17-mvp4-learned-bcc-proposal.md`：MVP-4 立项草案（S1–S6 待定，批准前禁止训练）。

参考（做到相关部分再读）：MVP-3 计划与台账（`2026-10-10-ffc-mvp3-charge-providers.md`、`.superpowers/sdd/2026-10-10-…/progress.md`）、A1/A2 计划与台账（`2026-10-11-ffc-a1a2-gate-scope-provenance.md`、`.superpowers/sdd/2026-10-11-ffc-a1a2/progress.md`）、`2026-10-11-crucible-next-steps.md`（A–F 总顺序，A4 小项表）、10-17 施工台账 `.superpowers/sdd/2026-10-17-unbuilt-modules/progress.md`。
逐日工作记录：`docs/history/HANDOVER_LOG_202610.md`。

## 当前状态

**版本 0.2.0（代码已超出 MVP-3，版本号未改，等用户定）。**

| 部分 | 状态 | 证据 |
|---|---|---|
| MVP-0～MVP-3 | **G0–G3 在当前源码上正式通过（2026-10-28 第 3 轮，决策 0026 之后，只靠配置文件）**，report.md 与 10-18 第 2 轮逐字相同 | 下文“正式运行”；本轮无计算交接 |
| G1.8 理论文档 | 已收口（四处修正，状态"已审阅"） | `docs/theory/hessian_methods.md` |
| ORCA .hess 读取器（路线图 D3） | 完成，按 ORCA 6.1.1 真实文件实现；频率与 ORCA 相对差 5.7e-8 | `experiments/20261017_orca_hess_format/`；`tests/fixtures/hessian/orca/` |
| C1 坐标合同 | 完成：`ffcoords/0.2` sidecar（带原子名），各命令 `--ir-in/--coords-in`，derive-bonded 输出 ir/coords/report | `ir/coords.py`、`cli.py` |
| C2 局部路线串联 | 单元级完成（2026-10-18 补上 torsion）：bonded → 电荷 → 二面角（`--reference mm` 替身）→ 三后端，参数与 Hessian 几何全程保留，`check-staleness` 不过期；再换电荷重拟合后依赖记录刷新 | `tests/unit/test_pipeline_chain.py::test_local_route_chain_with_torsions` |
| C3/C4 依赖记录与过期 | 完成：H_fixed / MM 扣除快照、`check-staleness`、`global-fit` CLI（下界命中不导出） | `bonded/staleness.py`、`cli.py` |
| C5 组合验收 | 单元/组合级完成；**正式外部组合验收未做**（10-15 §5，需报资源） | 同上 |
| MVP-5 fallback | 框架完成（gfnff → gfn2 → DFT 文件，质量门）；CLI 不暴露；调度研究未开始 | `hessian/fallback.py` |
| MVP-4 learned BCC | 只有骨架 + 立项草案；**研究未立项** | `charge/learned_bcc.py`、立项草案 |
| 感知层 / fflib（10-24，决策 0022） | 有单测，**还没接进 CLI 和管线**；10-28 修正了 bondcrit 窗口、improper 匹配、负 pn 三处语义（0025） | `perception/`、`fflib/` |
| MVP-6 分子入口 + 单命令 CLI（2026-11-11） | **完成**：enter/run CLI、preset（§17）；12 固件 SDF 重建等价 12/12（本机）；**G4 未建、端到端正式验收待 Gate 机** | `molecule/`、`typing/`、`presets.py`、`tests/unit/test_molecule_entry.py` |
| localhess（10-24，决策 0023） | experimental，CLI 可选；没有 oracle，只做了内部一致性检查 | `bonded/localhess.py` |
| CP2K provider（10-24，决策 0023） | 和官方 FeCl2 示例逐位一致；10-28 起 CLI 用文件几何（0025） | `hessian/cp2k.py` |
| IR 周期盒（10-24，决策 0023/0025） | 只作模拟盒元数据；三个导出器都写盒子，AMBER 只支持正交盒；sander 可读（10-28 起） | `ir/`、`exporters/` |
| B 二面角质量 | **B1/B2 完成（2026-10-17）**：1-2 缺陷由能量账确认为单遍拟合的跨键耦合（主要是 0-1 甲基的 GAFF2 项被 1-2 吸收）；nma 1-3 扫描滞后 57 kJ/mol 另列。**B3 顺序重拟合机制实验完成**：1 轮后 1-2 由 3.36 降到 0.44（GAFF2 2.70），全局残差 0.83 → 0.22，无键明显变差；样本内结果，未采用为默认。**顺序依赖已查（2026-10-18）**：全部顺序趋向同一不动点，第 1 轮拟合项差 ≤2.5e-2，第 2 轮 ≤1.2e-3，第 3 轮 ≤1e-4 kJ/mol | `jobs/torsion_diagnosis_20261017/report.md`、`jobs/torsion_b3_sequential_20261017/report.md` |

最近一次测试（2026-10-28，node-A，0025 修复后的源码）：`pytest -m "not external"` 579 passed（104 s，`experiments/20261028_review_pytest/pytest_not_external_after_fix.log`）；external 只跑了改动涉及的部分（见上）。

更早一次（2026-10-17 夜，node-A，去掉 Git 探测之后的源码）：**全量 561 个全部通过**——`pytest -m "not external"` 502 passed（85 s），`pytest -m external` 59 passed（546 s，其中 MACE 的 `test_derive_torsions_mace_ethanol` 415 s）。日志 `experiments/20261017_full_pytest/`。G0–G3 已于 2026-10-18 在当前源码上重跑并通过（见“正式运行”）。

## 正式运行

**最新：2026-10-28 第 3 轮，决策 0026 之后的源码，只靠用户配置文件，G0–G3 全部 passed**（11:17–11:43 JST；G0 894 s、G1 920 s、G2 838 s、G3 639 s；`jobs/gate_mvp{0,1,2,3}_20261028_3/`；report.md 与第 2 轮逐字相同；pytest 603 + 60，0 失败）。

前一轮：2026-10-28 第 2 轮，决策 0025 修复后的源码，G0–G3 全部 passed（node-A，非原作者账号，07:57–08:21 JST；G0 791 s、G1 921 s、G2 784 s、G3 641 s；`jobs/gate_mvp{0,1,2,3}_20261028_2/`，G2 `--resume jobs/gate_mvp2_20261018_2`，G3 `--scan-cache jobs/gate_mvp2_20261028_2/scans`）。四个 Gate 的 report.md 与 10-18 第 2 轮逐字相同；pytest 579 + 60，0 失败。同日第 1 轮（`_20261028`）因测试的 HOME 问题判为未通过，已被取代。

上一轮：2026-10-18 第 2 轮，可选顺序重拟合接入后的源码，G0–G3 全部 passed（node-A，三路并发，08:50–09:15 JST）。四次运行的 `source_sha256` 彼此一致（119 个文件），与运行后的工作区 0/119 不同；相对凌晨那轮改动的是 `torsion/mlp.py`、`torsion/mm.py`、`bonded/staleness.py`、`cli.py` 和相关测试，新增 `test_torsion_refinement.py`、`scripts/b3_order_dependence.py`。`git_commit` 为 null（决策 0021）。启动脚本 `jobs/launch_gates_20261018_2.sh`、`jobs/chain_g2_g3_20261018_2.sh`。

| Gate | 运行目录 | 结果 |
|---|---|---|
| G0 | `jobs/gate_mvp0_20261018_2/` | passed（837 s）：11 固件 × 4 构型，8 条腿 2596 行 0 失败 |
| G1 | `jobs/gate_mvp1_20261018_2/` | passed（841 s）：G1.1–G1.7，scope=full；Sobtop 三种方法 12/12；全局拟合自洽 6.03e-9 |
| G2 | `jobs/gate_mvp2_20261018_2/` | passed（853 s，`--resume jobs/gate_mvp2_20261013`）：12 根键全部 cached、mlp，0 兜底；metrics.csv 与 10-13 逐字相同；sander nma dihed 2.3e-4 只报告 |
| G3 | `jobs/gate_mvp3_20261018_2/` | passed（622 s，`--scan-cache jobs/gate_mvp2_20261018_2/scans`）：G3.1 ≤4.29e-11 e、compare 33/33、G3.4 逐项相同 |

四个 Gate 的 pytest 都是 524 + 59，0 失败。各项数值与凌晨那轮、10-13 那轮一致。Gate 走的是默认单遍；三轮顺序模式只由单元测试和 `experiments/20261018_refinement_verify/` 覆盖，不在 Gate 判据里。

上一轮：2026-10-18 凌晨（`jobs/gate_mvp{0,1,2,3}_20261018/`，117 个文件，可选模式接入前）；2026-10-13 A3（`jobs/gate_mvp{0,1,2,3}_20261013/`，验收报告 `docs/history/ACCEPTANCE_20261013_A3.md`）；更早：`gate_mvp0_20261005_3`、`gate_mvp1_20261005_4`、`gate_mvp2_20261010_2`、`gate_mvp3_20261010_2`（保留）。

其他运行目录全部保留不覆盖：
- `gate_mvp2_20261008_3/`：G2 首跑。它的 `scans/` 是原始 MACE 扫描缓存（10-13 的 G2 由它复制而来），**不要删**。
- `gate_mvp2_20261010/`、`gate_mvp3_20261010/`：调试冒烟（部分分子、跳过 pytest）。
- 更早的 G0/G1 运行见 CHANGELOG。

PID 记录：
- 2026-10-28 第 3 轮（决策 0026 之后，只靠配置文件；均已结束）：external pytest 2162489（`experiments/20261028_tool_config_verify/pids.txt`），`gate_mvp0.pids` 2162531，`gate_mvp1.pids` 2162568，`gate_mvp2.pids` 2162569（串联 bash）/ 2162766（G2）；启动脚本 `jobs/launch_gates_20261028_3.sh`、`jobs/chain_g2_g3_20261028_3.sh`、`experiments/20261028_tool_config_verify/run_external.sh`。
- 2026-10-28 第 2 轮（修好测试的 HOME 传递之后，均已结束）：`gate_mvp0.pids` 1656442，`gate_mvp1.pids` 1656443，`gate_mvp2.pids` 1656444（串联 bash）/ 1656487（G2）；启动脚本 `jobs/launch_gates_20261028_2.sh`、`jobs/chain_g2_g3_20261028_2.sh`。
- 2026-10-28（0025 修复后，node-A，均已结束）：`jobs/gate_mvp0.pids` 1293456，`gate_mvp1.pids` 1293459，`gate_mvp2.pids` 1293467（串联 bash）/ 1293517（G2），G3 的 PID 由串联脚本追加到 `gate_mvp3.pids`。启动脚本 `jobs/launch_gates_20261028.sh`、`jobs/chain_g2_g3_20261028.sh`（只对这些进程把 HOME 指向原作者家目录，原因见脚本注释）。
- 以下均已结束：
- 2026-10-18 第 2 轮：`jobs/gate_mvp0.pids` 1936619，`gate_mvp1.pids` 1936639，`gate_mvp2.pids` 1936640（串联 bash）/ 1936695（G2），`gate_mvp3.pids` 2169319。
- 2026-10-18：`jobs/gate_mvp0.pids` 1436597，`gate_mvp1.pids` 1436600，`gate_mvp2.pids` 1436609（串联 bash）/ 1436637（G2），`gate_mvp3.pids` 1663387。
- 2026-10-13：`jobs/gate_mvp0.pids` 3046691，`gate_mvp1.pids` 3046692，`gate_mvp2.pids` 84807 / 84820（node-B），`gate_mvp3.pids` 3325239；启动脚本 `jobs/launch_a3_20261013.sh`、`jobs/launch_g3_20261013.sh`（当时没用 `scripts/run_logged.sh`，因为它调用 `git -c safe.directory`；2026-10-17 起已去掉 git 探测，决策 0021）。
- `jobs/gate_mvp2.pids`：1389493 / 1389505 是 10-08 首跑，68591 / 68603 是 10-10 在 node-B 上的重跑。
- `jobs/gate_mvp3.pids`：4099417 / 4099447。

每对里前一个是外层 bash，后一个是 Gate 的 python。

## 待办

**下一步工作（按 10-15 计划，执行前报资源计划）：**
0. ~~0025 修复后重跑 G0–G3~~：10-28 第 2 轮全部通过（见顶部）。（当时非原作者账号要靠改 HOME 才能跑；决策 0026 之后改为各自写 `~/.config/crucible/tools.toml`。）
1. B3：机制实验与顺序依赖检查已做；可选 0–3 轮接口及换电荷同协议重放已于本轮准备，默认仍为单遍。下一步由计算接收方先跑新增回归与旧实验一致性检查，再做 B4 留出验证。有限轮按本键扫描更新不等价于联合最小二乘；全局拟合对照与默认策略不自动决定。见本轮无计算交接。
2. C 收尾：~~C2 串联测试补 torsion 一步~~（2026-10-18 已补）；C5 正式外部组合验收（10-15 §5）。
3. MVP-4 立项讨论（立项草案 S1–S6）。

**等用户定：**
4. MVP-2 计划 §7 第 2–4 条：扫描几何变形是否作兜底条件；同一键能量并到一个四元组的近似；`spin=1` 的含义。
5. nma 1-3 扫描滞后：用户 2026-10-17 定为先看分布、不设自动兜底阈值。分布已看（B3 报告末节）：432 帧中位数 0.001、90% 分位 0.01 kJ/mol；10 根键 ≤0.01，methyl_benzoate 1-2 为 0.41，nma 1-3 为 57.4（14 帧 > 1，集中在 ±60–120° 能垒区）。nma 1-3 单列为扫描质量问题，其"拟合改善"不算科学收益；是否重扫另定。
5b. **已修（2026-10-17，决策 0021）**：RunContext 删除 Git 探测，git_commit / git_dirty 记 null；run_logged.sh 同步移除探测。回归 8 passed，未重跑 Gate。
6. `derive-bonded --method` 是否默认 mseminario（A3 偏离 #8，设计 §6.1；10-15 计划把默认值列为接口建议）。
7. 版本号（0.2.0）要不要改。（后半"10-17 改动后是否重跑 G0–G3"已完成：10-28 第 3 轮已在最新源码上重跑并通过。）
8. 受约束松弛检查**阶段 2**（2026-11-01 计划 §5 已锁分子与判据）：新增 5 分子、需算 MACE 扫描 1–2.5 h，σ 按新扫描重取。阶段 1 已执行（见 2026-11-11 条目：停止条件均未触发，X 保留为候选诊断），阶段 2 做不做由用户定。

### 代码侧未做的事（2026-11-11 盘点；依据：结案验收 §1、README Roadmap、上表）

**已设计、代码未写**（决策 0024 的 MVP-7～9；MVP-6 已于 2026-11-11 完成，见顶部条目）：
1. ~~MVP-6 分子入口 + 单命令 CLI + preset~~：**已完成（2026-11-11）**；遗留：G4 Gate 脚本与正式端到端验收（需 Gate 机器，资源另报）、rigid preset 的 rigid_torsion 语义待拍板。
2. MVP-7 扫描质量：多起点松弛扫描 + 逐格点分支诊断（分支只诊断、永不拟合）；首个应用 nma 1-3 重扫。
3. MVP-8 平衡几何联合拟合：q⁰ 进拟合 + 梯度匹配（§6.4.3-1 预留接口落地），opt-in，对现行 global-fit 逐位退化。
4. MVP-9 多构象联合拟合 + charge→torsions→global-fit 交替协议与 ε 完成判据（§6.4.4–6.4.6）。
4. MVP-8 平衡几何联合拟合：q⁰ 进拟合 + 梯度匹配（§6.4.3-1 预留接口落地），opt-in，对现行 global-fit 逐位退化。
5. MVP-9 多构象联合拟合 + charge→torsions→global-fit 交替协议与 ε 完成判据（§6.4.4–6.4.6）。

**只有骨架/框架，等研究立项（缺研究决定，不缺代码）**：
6. MVP-4 learned BCC：`charge/learned_bcc.py` 骨架，无训练模型、未接 CLI；立项草案 S1–S6 未定，批准前禁训练。
7. MVP-5 fallback 调度：`hessian/fallback.py` 阶梯已搭（gfnff→gfn2→DFT + 质量门），CLI 不暴露，§16 第 1 级"可信库参数"未接，不确定度调度未研究。

**代码在、未接进管线（集成债）**：
8. 感知层 `perception/` 与 `fflib/`（决策 0022）：有单测，未接 CLI 和管线（0025 修了 bondcrit 窗口、improper 匹配、负 pn 三处语义）。
9. `localhess`（0023）：experimental、CLI 可选，无 oracle，只做过内部一致性检查。
10. Level 3–4 验证基准（设计 §20）：完全没有（结案验收 #31）。

**小项**：A4 六项已完成（2026-11-11：重复 `find_rotors`、畸形 mol2 报错、`step_deg` 提前检查、TorsionRefitError 消息、扫描异常分类、G1 help）；A4 剩余：可转键对称类报告（需小规格）、主进程线程上限时机（动 Gate 语义）、甲铵正离子 MACE 回归（需 MACE 机器，不加未经运行的 external 测试）；RESP provider 推迟（0018）；OpenMM 原生 CustomTorsionForce 模式（计划不做，MVP-2 计划 §1 第 7 条）；TorsionDrive 2-D opt-in（装环境是显式决定）。

**验证/计算类（收尾用，不是写代码）**：B4 留出验证；C5 正式外部组合验收（需报资源）；G4（MVP-6 的 Gate，需 Gate 机器）；松弛检查阶段 2（等上面第 8 条拍板）。

**已知小问题（不阻塞）：**
- A4 未做的小项（重复 `find_rotors`、畸形 mol2 报错、`step_deg` 提前检查等），表见 `2026-10-11-crucible-next-steps.md` §3 A4。
- MVP-2 其余 Minor 见 `.superpowers/sdd/2026-10-08-ffc-mvp2-mlp-torsion/progress.md`。
- 非 external 的 pytest 在 node-A 高负载时要 72–131 s，超过 MVP-0 §1.5 的"1 分钟"约定；空闲时 63–78 s。

## 待用户提交的草稿（`.commit_drafts/`）

按 `.git/logs/HEAD`，1–18 号草稿都已在 2026-10-25 17:25 之前提交；作废的那几份本来就不执行。
**状态更新（2026-11-25 整理执行）：19–21 号（review-fixes / tool-config / relaxed-check-phase1）已并入 2026-11-11 03:02 的提交（`f20f4cd`：sobtop 吸收 + 决策 0022–0025；`88e370d`：tool config 决策 0026 + relaxed-check 脚本；`c4216ce`：实验记录与 gate manifests）。22–24 号于 2026-11-25 按序执行，以下清单保留作历史记录：**

19. `20261028_review-fixes.md`：10-24 代码审阅修复（决策 0025）及相应文档。**入库前要先脱敏**（运行记录里有真实路径，见 CHANGELOG 和 0025）。
20. `20261028_tool-config.md`：外部工具位置配置（决策 0026）。按 19 → 20 的顺序执行。
21. `20261111_relaxed-check-phase1.md`：受约束松弛检查阶段 1（脚本 + 报告 + 计划状态 + HANDOVER/CHANGELOG）。
22. `20261111_mvp6-molecule-entry.md`：MVP-6 分子入口 + 单命令 CLI + preset + A4 快修（新模块/CLI/测试/计划文档/HANDOVER/CHANGELOG）。21 → 22 顺序执行。
23. `20261122_batch1-code.md`：第一批代码层（覆盖集诊断 / MVP-7 多起点 / 独立验证装配 + 52 个新单测；全量回归按用户指令未跑）。22 → 23 顺序执行。
24. `20261125_1113_coverage-run-saddle-restart-ion-gfn2.md`：覆盖集诊断收口（xtb_sp/fallback/scan_probe、驱动三个潜伏错误、34 分子清单、多起点报告）+ 决策 0027 鞍点重启 + 决策 0028 离子 GFN2 参考势与 `run` 参数修复 + 覆盖集结果（32/2/0）。已脱敏，相关测试 129 passed。23 → 24 顺序执行。

## 结论备忘

- **MVP-1**：
  - Sobtop oracle 12/12：mSeminario（决策 0011）、m2Seminario（C2，0011）、DRIH（原子单位，0012）。
  - 全局拟合自洽误差 6.0e-9。
  - 已知局限：环丙烷 / GFN-FF / m2Seminario 的环角 k ≈ −33.6（0014）。
  - 频率 RMSE 中位数：GFN-FF 下全局拟合最好（57 cm⁻¹）。
- **MVP-2**：
  - 用 MACE-POLAR-1-M 经 ASE（不经 openmm-ml，0016），10° 约束松弛扫描，Fourier 拟合 n=1..6，相位任意。
  - 扫描耗时：每根键 222–2610 s，12 根键 8 进程 × 4 线程共 52 min。
- **MVP-3**：
  - AM1-BCC 与固件做法完全一致（`antechamber -c bcc -at gaff2 -pf y`，再均摊残差）。
  - NAGL 与 AM1-BCC 逐原子最多差 0.040 e。
  - 换电荷后用缓存的扫描重拟合二面角，不重扫。
  - 拓扑检查比对原子序数和键的连接关系（审阅后加强，偏离计划 §2.3 原文）。
- **10-17 管线（C）**：
  - 依赖用逐项快照记录、不用 hash；MLP 的 MM 扣除快照**不含其他键的 proper**（单遍拟合的已知近似，归 B2）。
  - global_fit 排在 charge/torsion 之后，但它改 bond/angle 会让 MLP 的 MM 扣除过期：按设计顺序走完，`check-staleness` 一定报 mlp 过期（按代码推断，未实跑这条完整链），这是 C4 规定的"标记过期"，不是 bug；耦合策略待 B2/B3。
- **验收约定**：
  - 0013：mixed GROMACS 用 fp32 容差，乙腈不进 Gate；
  - 0015：AMBER 原生读写；
  - 0017：MLP 二面角 IR 以 OpenMM 为准，sander 只报告。

## 环境

- node-A（本机，40 核，与他人共用）；node-B（8 核，同一 NFS，`source scripts/env.sh` 后 gmx/sander/MACE 都能用）。注意：node-B 的 NFS 属性缓存很长，在 node-A 上看 node-B 写的 log 可能是旧的，要到 node-B 上直接读。
- `source scripts/env.sh`（openmm_dev，Python 3.12，**生产环境**，不装不改）；ffcompiler 用 `PYTHONPATH=src`。
- GROMACS 2026.3 mixed（`FFC_GMX`）、xtb 6.7.1（`FFC_XTB`）、antechamber（`FFC_ANTECHAMBER`，默认 PATH）、Sobtop（`tools/sobtop/`，`OMP_STACKSIZE=64M OMP_NUM_THREADS=4`）。
- ORCA 6.1.1（node-A）：`/home/user/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2/orca`（10-17 运行时设了 `LD_LIBRARY_PATH` 含该目录与其 `lib/`，是否必需没试）。`~/.local/bin/orca` 是 orca-ide，**不是**量化程序。本项目只用它生成过 .hess 样例。
- MACE 模型 `/home/user/MLP/mace/MACE-POLAR-1-M.model`；NAGL 模型 `openff-gnn-am1bcc-1.0.0.pt`（openff-nagl-models 2025.9.0）。

## 本地 Git

工作区 `/home/user/CRUCIBLE`，git 元数据在工作区 `.git/`，远端 `origin` = `https://github.com/Cedrus810/CRUCIBLE.git`。
工作区在 NFS 上，agent 不执行任何 git 命令；提交写成 `.commit_drafts/` 草稿，由用户在实际机器上提交。
