# CRUCIBLE（原 forcefield_compiler）交接状态 — 2026-10-25（Asia/Tokyo）

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
- `docs/superpowers/plans/2026-10-08-ffc-mvp2-mlp-torsion.md` §7：MVP-2 还没定的问题（与 B 相关）。
- `docs/superpowers/plans/2026-10-17-mvp4-learned-bcc-proposal.md`：MVP-4 立项草案（S1–S6 待定，批准前禁止训练）。

参考（做到相关部分再读）：MVP-3 计划与台账（`2026-10-10-ffc-mvp3-charge-providers.md`、`.superpowers/sdd/2026-10-10-…/progress.md`）、A1/A2 计划与台账（`2026-10-11-ffc-a1a2-gate-scope-provenance.md`、`.superpowers/sdd/2026-10-11-ffc-a1a2/progress.md`）、`2026-10-11-crucible-next-steps.md`（A–F 总顺序，A4 小项表）、10-17 施工台账 `.superpowers/sdd/2026-10-17-unbuilt-modules/progress.md`。
逐日工作记录：`docs/history/HANDOVER_LOG_202610.md`。

## 当前状态

**版本 0.2.0（代码已超出 MVP-3，版本号未改，等用户定）。**

| 部分 | 状态 | 证据 |
|---|---|---|
| MVP-0～MVP-3 | **G0–G3 在当前源码上正式通过（2026-10-18 第 2 轮，含可选顺序重拟合）**，数值与之前各轮一致 | 下文“正式运行”；本轮无计算交接 |
| G1.8 理论文档 | 已收口（四处修正，状态"已审阅"） | `docs/theory/hessian_methods.md` |
| ORCA .hess 读取器（路线图 D3） | 完成，按 ORCA 6.1.1 真实文件实现；频率与 ORCA 相对差 5.7e-8 | `experiments/20261017_orca_hess_format/`；`tests/fixtures/hessian/orca/` |
| C1 坐标合同 | 完成：`ffcoords/0.2` sidecar（带原子名），各命令 `--ir-in/--coords-in`，derive-bonded 输出 ir/coords/report | `ir/coords.py`、`cli.py` |
| C2 局部路线串联 | 单元级完成（2026-10-18 补上 torsion）：bonded → 电荷 → 二面角（`--reference mm` 替身）→ 三后端，参数与 Hessian 几何全程保留，`check-staleness` 不过期；再换电荷重拟合后依赖记录刷新 | `tests/unit/test_pipeline_chain.py::test_local_route_chain_with_torsions` |
| C3/C4 依赖记录与过期 | 完成：H_fixed / MM 扣除快照、`check-staleness`、`global-fit` CLI（下界命中不导出） | `bonded/staleness.py`、`cli.py` |
| C5 组合验收 | 单元/组合级完成；**正式外部组合验收未做**（10-15 §5，需报资源） | 同上 |
| MVP-5 fallback | 框架完成（gfnff → gfn2 → DFT 文件，质量门）；CLI 不暴露；调度研究未开始 | `hessian/fallback.py` |
| MVP-4 learned BCC | 只有骨架 + 立项草案；**研究未立项** | `charge/learned_bcc.py`、立项草案 |
| B 二面角质量 | **B1/B2 完成（2026-10-17）**：1-2 缺陷由能量账确认为单遍拟合的跨键耦合（主要是 0-1 甲基的 GAFF2 项被 1-2 吸收）；nma 1-3 扫描滞后 57 kJ/mol 另列。**B3 顺序重拟合机制实验完成**：1 轮后 1-2 由 3.36 降到 0.44（GAFF2 2.70），全局残差 0.83 → 0.22，无键明显变差；样本内结果，未采用为默认。**顺序依赖已查（2026-10-18）**：全部顺序趋向同一不动点，第 1 轮拟合项差 ≤2.5e-2，第 2 轮 ≤1.2e-3，第 3 轮 ≤1e-4 kJ/mol | `jobs/torsion_diagnosis_20261017/report.md`、`jobs/torsion_b3_sequential_20261017/report.md` |

最近一次测试（2026-10-17 夜，node-A，去掉 Git 探测之后的源码）：**全量 561 个全部通过**——`pytest -m "not external"` 502 passed（85 s），`pytest -m external` 59 passed（546 s，其中 MACE 的 `test_derive_torsions_mace_ethanol` 415 s）。日志 `experiments/20261017_full_pytest/`。G0–G3 已于 2026-10-18 在当前源码上重跑并通过（见“正式运行”）。

## 正式运行

**最新：2026-10-18 第 2 轮，可选顺序重拟合接入后的源码，G0–G3 全部 passed**（node-A，三路并发，08:50–09:15 JST）。四次运行的 `source_sha256` 彼此一致（119 个文件），与运行后的工作区 0/119 不同；相对凌晨那轮改动的是 `torsion/mlp.py`、`torsion/mm.py`、`bonded/staleness.py`、`cli.py` 和相关测试，新增 `test_torsion_refinement.py`、`scripts/b3_order_dependence.py`。`git_commit` 为 null（决策 0021）。启动脚本 `jobs/launch_gates_20261018_2.sh`、`jobs/chain_g2_g3_20261018_2.sh`。

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

PID 记录（均已结束）：
- 2026-10-18 第 2 轮：`jobs/gate_mvp0.pids` 1936619，`gate_mvp1.pids` 1936639，`gate_mvp2.pids` 1936640（串联 bash）/ 1936695（G2），`gate_mvp3.pids` 2169319。
- 2026-10-18：`jobs/gate_mvp0.pids` 1436597，`gate_mvp1.pids` 1436600，`gate_mvp2.pids` 1436609（串联 bash）/ 1436637（G2），`gate_mvp3.pids` 1663387。
- 2026-10-13：`jobs/gate_mvp0.pids` 3046691，`gate_mvp1.pids` 3046692，`gate_mvp2.pids` 84807 / 84820（node-B），`gate_mvp3.pids` 3325239；启动脚本 `jobs/launch_a3_20261013.sh`、`jobs/launch_g3_20261013.sh`（当时没用 `scripts/run_logged.sh`，因为它调用 `git -c safe.directory`；2026-10-17 起已去掉 git 探测，决策 0021）。
- `jobs/gate_mvp2.pids`：1389493 / 1389505 是 10-08 首跑，68591 / 68603 是 10-10 在 node-B 上的重跑。
- `jobs/gate_mvp3.pids`：4099417 / 4099447。

每对里前一个是外层 bash，后一个是 Gate 的 python。

## 待办

**下一步工作（按 10-15 计划，执行前报资源计划）：**
1. B3：机制实验与顺序依赖检查已做；可选 0–3 轮接口及换电荷同协议重放已于本轮准备，默认仍为单遍。下一步由计算接收方先跑新增回归与旧实验一致性检查，再做 B4 留出验证。有限轮按本键扫描更新不等价于联合最小二乘；全局拟合对照与默认策略不自动决定。见本轮无计算交接。
2. C 收尾：~~C2 串联测试补 torsion 一步~~（2026-10-18 已补）；C5 正式外部组合验收（10-15 §5）。
3. MVP-4 立项讨论（立项草案 S1–S6）。

**等用户定：**
4. MVP-2 计划 §7 第 2–4 条：扫描几何变形是否作兜底条件；同一键能量并到一个四元组的近似；`spin=1` 的含义。
5. nma 1-3 扫描滞后：用户 2026-10-17 定为先看分布、不设自动兜底阈值。分布已看（B3 报告末节）：432 帧中位数 0.001、90% 分位 0.01 kJ/mol；10 根键 ≤0.01，methyl_benzoate 1-2 为 0.41，nma 1-3 为 57.4（14 帧 > 1，集中在 ±60–120° 能垒区）。nma 1-3 单列为扫描质量问题，其"拟合改善"不算科学收益；是否重扫另定。
5b. **已修（2026-10-17，决策 0021）**：RunContext 删除 Git 探测，git_commit / git_dirty 记 null；run_logged.sh 同步移除探测。回归 8 passed，未重跑 Gate。
6. `derive-bonded --method` 是否默认 mseminario（A3 偏离 #8，设计 §6.1；10-15 计划把默认值列为接口建议）。
7. 版本号（0.2.0）要不要改；10-17 改动后是否重跑 G0–G3（同一份源码身份）。

**已知小问题（不阻塞）：**
- A4 未做的小项（重复 `find_rotors`、畸形 mol2 报错、`step_deg` 提前检查等），表见 `2026-10-11-crucible-next-steps.md` §3 A4。
- MVP-2 其余 Minor 见 `.superpowers/sdd/2026-10-08-ffc-mvp2-mlp-torsion/progress.md`。
- 非 external 的 pytest 在 node-A 高负载时要 72–131 s，超过 MVP-0 §1.5 的"1 分钟"约定；空闲时 63–78 s。

## 待用户提交的草稿（`.commit_drafts/`）

按时间顺序执行；是否已提交过由用户核对（agent 不跑 git）。

1. `20261008_*`（10 份）：10-08 的改名、T24–T30、评审修复、G2 首跑。
2. `20261010_1200_docs-tidy.md`、`20261010_1715_g2-openmm-reference.md`、`20261010_1800_mvp3-plan.md`、`20261010_1945_mvp3-charge-providers.md`。
3. `20261011_2330_a1a2-gate-scope-provenance.md`。
4. `20261013_1340_refit-fallback-cache-paths.md`、`20261013_1530_a3-acceptance.md`。
5. `20261017_1210_g18-doc-closure.md`：G1.8 文档收口（与 6 不重叠）。
6. `20261017_2000_unbuilt-modules-accepted-and-fixed.md`：10-17 施工 + 当晚验收修复，一个提交。
7. `20261017_2100_docs-tidy.md`：本次文档整理。
8. `20261017_runtime-no-git.md`：移除 RunContext/日志包装器的 Git 探测，8 条回归通过；B3 联合目标等价条件澄清。
9. `20261017_2300_b1b2b3-torsion-diagnosis.md`：B1/B2 诊断脚本与报告、B3 顺序重拟合机制实验、滞后分布，及相应文档。
10. `20261018_0215_gates-rerun-current-source.md`：全量 pytest 与 G0–G3 在当前源码上重跑（全部通过）及相应文档。
11. `20261018_1000_b3-order-and-c2-torsion-chain.md`：B3 顺序依赖检查；C2 含 torsion 的串联测试。
12. `20261018_optional-sequential-refinement.md`：可选顺序更新接口、换电荷协议重放、报告/依赖与新增测试；2026-10-18 接收方验证结果已补进草稿。
13. `20261018_0920_gates-rerun-after-refinement.md`：可选顺序重拟合接入后的 G0–G3 第 2 轮重跑（全部通过）及相应文档。
14. `20261018_1000_closure-acceptance.md`：结案验收报告（放在最后执行；提交后把哈希补记进报告）。
15. `20261008_allen-crosscheck.md`：决策 0011 的 Allen 2018 公式对照补记与原文 PDF（10-08 内容的遗留认领；PDF 入库与否执行时定）。
16. `20261024_update-round2.md`：感知层、参数数据读取、CP2K 读取器、IR 周期盒、localhess 方法及决策 0022/0023（**提交信息按用户指示只写 "update"**；会连带 CHANGELOG 的 10-25 条目）。
17. `20261025_next-phase-design.md`：设计稿下一阶段增补（§3.1、§6.4.4–6.4.6、§10.4–10.5、MVP-6～9）与决策 0024；仅文档。
18. `20261025_readme-overhaul.md`：三份 README 全量重写（EN/zh-CN/ja；含 10-24 的 README 改动）；仅文档。
    执行顺序 15 → 16 → 17 → 18。
- **作废、不要执行**：`20261017_1240`、`_1500`、`_1540`、`_1620`、`_1650`、`_1720`（由 6 取代，文件开头已标注）。

几份草稿涉及同一批文件（README、HANDOVER、CHANGELOG、roadmap）的工作区状态：先执行的那条会把后面的改动一起带上，后面的命令仍可执行。

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
