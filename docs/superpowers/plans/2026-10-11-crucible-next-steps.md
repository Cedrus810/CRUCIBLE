# CRUCIBLE 后续工作整理（2026-10-11，Asia/Tokyo）

> **2026-10-17：** C 的工程部分与 MVP-4 骨架、MVP-5 框架已搭（`2026-10-17-unbuilt-modules.md`）；当前状态看 `HANDOVER.md`。
>
> **2026-10-15 状态入口：** A1/A2、决策 0020 的运行细节与 A3 四道 Gate 重跑均已完成；本次复核四份同源 manifest 的 104 个已列文件与当前内容一致。G1.8 还差已审出的 4 处文档修改。下文 §1 是 10-11 快照，不再代表当前源码对应状态；历史建议中“先做 A1/A2”“A3 未做”等措辞也已过时。后续 B/C 的具体任务、依赖及完成标准见 [10-15 细化清单](2026-10-15-crucible-bc-detailed-next-steps.md)，验收事实见 [A3 报告](../../history/ACCEPTANCE_20261013_A3.md)。

本文是建议顺序与工作清单，不是新的验收标准或已批准的施工计划。依据现行设计稿、HANDOVER、决策 0016–0018、MVP-2/3 执行台账，以及本次读取的 Gate 产物和关键源码。此次只整理文档，未启动测试、扫描或 benchmark，未修改实现及已有状态记录。

## 1. 当前基线

MVP-0 到 MVP-3 均有正式 Gate 通过记录；MVP-4 learned BCC、MVP-5 自动 fallback/不确定度调度尚未开始。版本号仍为 0.2.0。

本次直接核对四份 acceptance.json，均为 passed=true；另按各 manifest 已记录的 source_sha256 清单逐文件比较当前内容：

| Gate | 正式运行 | 与当前内容不同的清单条目 | 版本身份 |
|---|---|---:|---|
| G0 | jobs/gate_mvp0_20261005_3 | 9 / 75 | 旧提交 c20d6a2，git_dirty=true |
| G1 | jobs/gate_mvp1_20261005_4 | 9 / 75 | 旧提交 c20d6a2，git_dirty=true |
| G2 | jobs/gate_mvp2_20261010_2 | 5 / 91 | git_commit / git_dirty 均为 null |
| G3 | jobs/gate_mvp3_20261010_2 | 0 / 104 | git_commit / git_dirty 均为 null |

比较限定于 manifest 已列出的文件；不是 Git 工作区检查，也不证明清单覆盖了所有新文件。G3 清单与当前内容一致，是有效的内容对应证据；干净提交身份仍待补齐。G0/G1/G2 历史通过记录保留，但不能据此宣布当前整套源码已全部重新验收。

已确认的工程能力包括：IR、三后端导出、Hessian bond/angle 推导与全局拟合函数、MLP 二面角、loaded/AM1-BCC/NAGL 电荷、换电荷后的缓存重拟合。跨引擎等价和重拟合自洽尚不足以说明整体力场比 GAFF2 更准确。

## 2. 建议推进顺序

| 顺序 | 工作包 | 完成时应交付 |
|---|---|---|
| A | 收尾验收与来源追踪 | Gate 覆盖范围明确，参数来源可追溯，历史扫描缓存不被改写 |
| B | 诊断二面角质量 | 苯甲酸甲酯反例的原因分析；质量/几何/近似处理的决策草案 |
| C | 串联完整参数化管线 | 同一最终 IR 保留全部模块结果；全局拟合依赖和更新策略明确 |
| D | 建立独立科学 benchmark | 冻结的分子、参考数据、拆分、指标及消融方案 |
| E | MVP-4 learned BCC 立项 | 基础电荷、模型适用域、训练目标、数据来源和 G4 草案 |
| F | MVP-5 调度 | 经校准的质量/不确定度规则及带原因记录的 fallback graph |

A 的小修与 B 的缓存诊断可并行安排；C 的依赖策略需在组合 benchmark 和 MVP-4 正式训练/评价前确定。D 的方案、E 的立项讨论可以提前，实际科学结论依赖可追溯的工程基线。不要因先做收尾而把所有低影响条目都变成研究阻塞项。

## 3. A：近期优先收尾

### A1. Gate 覆盖与历史产物

- G2/G3 使用 --molecules 子集时，正式通过判定应仍核对完整规定矩阵；区分冒烟通过与正式 Gate 通过。manifest.debug_flags 本身不能防止误读 passed=true。**【2026-10-11 已做：T38，决策 0019；G1 `--providers` 一并处理】**
- MVP-2 --resume 目前可能往旧 scans/ 写入。建议读取既有缓存，新扫描写入新运行目录，并记录来源与文件身份。**【已做：T39，复制到本次目录，manifest 记 scan_cache_seeded_from】**
- 不覆盖旧 manifest、acceptance 或 report；不删除 jobs/gate_mvp2_20261008_3/scans，该缓存被 G2/G3 共用。

建议回归用例：子集不能取得正式验收；缺比较腿/缺测试不能通过；复用缓存后旧产物内容不变。保留现有数值容差和决策 0017。

### A2. 电荷 provider 与 provenance

- assign_charges 的中央接口检查电荷形状、原子数和有限性；当前 zip 非 strict，不能把防护只交给各 provider，尤其后续还要接 learned BCC。**【已做：T40】**
- LoadedCharges() 保留原始电荷计算来源，并另外记录本次载入/残差处理，而不是将历史来源覆盖成 loaded。**【已做：T40，provenance["charge"]["previous"]】**
- 重拟合报告保留原兜底原因；区分原扫描失败与本次缓存缺失。**【已做：T41】**
- --reference mm 的缓存需包含参考 IR/Hamiltonian 身份；换电荷等改变 MM 参考势的操作不得命中旧参考能缓存。MACE 参考扫描与最终 MM 参数分开标识。**【已做：T42，只记电荷（字面值，不用 hash）；bonded/LJ 改动留给 C】**

建议回归用例：自定义 provider 多/少返回电荷时明确失败；重复 loaded 后仍可追溯到原 provider；换参考 MM 后缓存失效；重拟合断言确实检查更新后的 provenance。

### A3. 可追溯验收与发布收尾

**【2026-10-13：G0–G3 已在同一份源码上重跑并全部通过；设计稿逐条验收与 G1.8 审阅见 `docs/history/ACCEPTANCE_20261013_A3.md`。G1.8 待修 4 处文档措辞；版本号待用户定】**

- 完成 G1.8 理论文档人工审阅；已有公式编号补充不能代替审阅记录。
- 先完成计划内修复并固定源码，再由用户在实际机器上提交。agent 不在本 NFS 工作区运行 Git，提交说明沿用 .commit_drafts/。
- 在满足项目 Git/运行规则的执行方式下重跑完整 G0/G1/G2/G3，记录源码、输入、模型及缓存身份；若 Git 信息不可获得，明确记录替代的源码身份，不能将 null 当成干净提交。
- 验收对应后再统一版本号、CHANGELOG、三语 README 与状态摘要；当前 README 功能表仍写“两道验收 Gate”，可纳入文档收尾。

### A4. 可合并处理的低影响遗留

MVP-3：重拟合测试的弱断言、--scan-cache 错误提示、重复 find_rotors、畸形 mol2 的错误信息。台账中“AM1-BCC 检查位置过晚”已被 Important 1 的修复覆盖，无需重复排工。

MVP-2：提前检查 step_deg、主进程线程上限的设置时机、ASE/模型异常分类、甲铵正离子专门回归、SDF 解析失败提示、可转键对称类报告。测试耗时超一分钟先确认负载与瓶颈，避免为时限削减必要覆盖。

**2026-10-13 逐条核实（读现行代码，未改）：**

| 条目 | 现状 | 位置 |
|---|---|---|
| 重拟合测试弱断言 | 已修（T40）：重拟合后断言 `charge_method == "shifted-test"` 且项已变 | `tests/unit/test_torsion_refit.py:94-95` |
| --scan-cache 错误提示 | 一半：G3 已做（决策 0020）；`TorsionRefitError` 的"needs scan_cache"一句没写 CLI 参数名 | `src/ffcompiler/torsion/mlp.py:240` |
| 重复 find_rotors | 仍在，只在报错分支多算一次 | `scripts/gate_mvp3.py:205-206` |
| 畸形 mol2 报错 | 仍在：列数不足 / 非数字时原样抛 IndexError / ValueError | `src/ffcompiler/charge/loaded.py:22-23` |
| 提前检查 step_deg | 一半：整除 360 在扫描开头查；点数 < 13 到 `fit_fourier` 才报（整根键已扫完）；默认 10° 不触发 | `torsion/scan.py:75-77`、`torsion/fourier.py:42` |
| 主进程线程上限时机 | 仍在：模块顶部已 import numpy，`_init_worker` 在 main 里才设环境变量；子进程有效 | `scripts/gate_mvp2.py:35` |
| ASE/模型异常分类 | 仍在：LBFGS 与取能量/力里任何异常都包成"参考势求值失败"→ `mlp_error`（原文保留） | `torsion/scan.py:85-94` |
| 甲铵正离子回归 | 仍在：MACE 带电测试只有醋酸根 | `tests/unit/test_torsion_mace.py:64` |
| SDF 解析失败提示 | 实际不触发：`from_amber` 先读 SDF，失败报 `cannot read SDF`；后面二次解析只在同一文件已读成功后走到 | `importers/amber_prmtop.py:53-55,142` |
| 可转键对称类报告 | 仍在：按扫描四元组算类，只影响报告 | `torsion/rotors.py:91` |
| （A1/A2 Minor）rotors 列表共享 | 仍在：report 与 IR provenance 是同一 list | `torsion/mlp.py:224,269-270` |
| （A1/A2 Minor）G1 `--providers` help | 仍在：只写了 0019 注释 | `scripts/gate_mvp1.py:537` |

以上都不改变 Gate 结果或参数值；随 B（二面角相关）/ C（CLI、电荷）的改动一起处理。

A1/A2 施工记录：`2026-10-11-ffc-a1a2-gate-scope-provenance.md` 与 `.superpowers/sdd/2026-10-11-ffc-a1a2/progress.md`。A3（G1.8 审阅、提交后重跑 G0–G3）、A4 未做。

2026-10-13：A1/A2 最终审阅的三条推迟项按决策 0020 收掉（重拟合不升级原兜底键；G2 `--resume` 路径写错报错；G3 `--scan-cache` 必填且须是目录，A4 的“--scan-cache 错误提示”一并完成）。A3 的 G3 用同一轮 G2 的 `scans/`。

## 4. B：二面角质量诊断

苯甲酸甲酯 1-2 键的 G2 实测：拟合残差 0.604426 kJ/mol，最终剖面 RMSE 3.35511，GAFF2 为 2.69745；最大滞后 0.410472 kJ/mol、键长变化 2.88695 pm。当前用户决定是记为已知缺陷，修法尚未定。

建议先用既有缓存做可解释诊断：

1. 在相同扫描构型上分解目标键项、其他二面角项、bond/angle 和非键贡献；检查逐键拟合后替换其他键是否改变最终剖面。
2. 同时报告原有“各减自身最小值”指标、最佳常数平移后的误差、能垒及极小值位置。新指标作为补充，保留历史指标与结果。
3. 比较逐键单遍与有限轮重拟合；如观察到收益，先锁定收敛条件、轮数上限与失败处理，再决定采用。
4. 在其他键和未参与拟合的构型上检查，防止为一个反例制定特殊修法。

尚待定的三项：扫描几何变形是否触发兜底及阈值；同一中心键能量集中到一个四元组的适用范围；MACE 接口 spin=1 的确切约定。spin 语义应核对实际安装接口/权威说明，目前本整理不替它作结论。

NMA 几何偏差是此前审查的候选诊断项；旧记录为约 45–63 pm，未在本次复核，不应当作当前源版本已确认的问题。若继续诊断，先核对构象、扭转角和原子映射。

## 5. C：把模块串成完整管线

现行设计的入口是分子结构，经 typing/base FF 到最终 IR。目前公共 CLI 主要从已参数化的 AMBER prmtop/inpcrd + SDF 开始；任意 SDF 的 GAFF2 初始化仍需形成明确工作流或单列后续任务。

当前 derive-bonded 从原始 AMBER 输入起步，没有 --ir-in/IR JSON 输出选项；global_fit 已有 Python 函数，但不在 derive-bonded 的 --method 选项里。单独模块通过 Gate，尚不能代替整条命令链保持所有已生成参数的验证。

建议先支持已有 AMBER 输入的完整流程，再扩展分子结构入口与 presets：

- 局部 Hessian 路线：基础 IR → bond/angle → 最终电荷 → 二面角 → 同一最终 IR → 三后端比较。
- 全局拟合路线：根据设计 §6.4，先固定 charge/torsion，再扣 H_fixed 拟合 bond/angle。
- 依赖更新：换电荷已会重拟合 MLP 二面角；全局拟合的 bond/angle 也依赖 charge/torsion。反过来改 bond/angle 又可能改变二面角残差，需要事先确定分阶段冻结、失效重算或受控迭代策略，不能默认一次排序就满足所有依赖。
- 参数来源记录实际输入的 charge、torsion、Hessian 和 IR 身份；失效状态须可识别。

组合验收至少覆盖一条局部路线、一条 global_fit 路线，以及参数更新后的依赖处理；输出中同时保留拟合后的 bond/angle、电荷和二面角。

## 6. D：独立科学验证

先冻结一个覆盖普通中性、离子、共轭/耦合扭转、含硫及环系的小分子矩阵，规模由计算预算决定。参考数据来源与资源安排另行确定。

第一层比较独立参考的优化几何、频率、留出构象能量差、力和扭转能垒。把 GAFF2、仅换 bonded、仅换 torsion、仅换 charge、完整组合分开比较；固定 LJ 与 1-4 约定，记录各比较所需的二面角重拟合。

拟合所用 MACE/Hessian 的自洽结果继续报告；它们不能单独承担独立准确度结论。HFE/溶剂转移、片段相互作用，再到 RBFE/ABFE，按设计 §20 分层推进。科学 benchmark 的阈值先锁定，不能用工程 Gate 的 passed 代替。

## 7. E/F：下一期研究与长期路线

MVP-4 单独立项，首先定下：

1. 基础电荷：AM1-BCC、EEQ 或其他 prior；NAGL 的角色是基线还是初始化。
2. 第一版适用域：元素、总电荷、自旋、构象/环境依赖和最终固定电荷的定义。
3. 数据与标签：来源、QM 层级、训练/验证/测试拆分，避免近重复/相同骨架泄漏；AM1-BCC/NAGL 相互接近不等于物理真值。
4. 第一批训练目标、独立评价目标、固定 GAFF2 LJ 的边界及 G4 草案。遵循设计的多目标方向，可分批纳入目标，不必第一版同时跑完所有昂贵性质。

设计建议的核心形式是反对称键电荷转移修正，天然保持总电荷。实现与研究验收分开：守恒、对称性和适用域是工程要求，独立物理性能是研究要求。

MVP-5 保持路线图中的最后阶段：先通过数据校准不确定度/质量信号，再制定 GFN-FF → GFN2 → 高等级参考等升级策略；已有 MLP → GAFF2 的局部兜底继续保留。

RESP 继续后置（0018）；sander/pmemd 的相位问题维持报告范围（0017）；ORCA Hessian reader 等真实样例；LJ 联合学习、Hessian 二面角和 enhanced OpenMM 模式不自动加入近期工作包。

## 8. 本整理的资料入口

- [当前交接](../../../HANDOVER.md)、[现行设计](../../../forcefield_compiler_design.md)
- [既有路线图](2026-10-04-ffc-roadmap.md)
- [MVP-2 计划](2026-10-08-ffc-mvp2-mlp-torsion.md)、[MVP-3 计划](2026-10-10-ffc-mvp3-charge-providers.md)
- [MVP-2 台账](../../../.superpowers/sdd/2026-10-08-ffc-mvp2-mlp-torsion/progress.md)、[MVP-3 台账](../../../.superpowers/sdd/2026-10-10-ffc-mvp3-charge-providers/progress.md)
- [决策索引](../../decisions/README.md)

建议下一次具体施工从 A1/A2 的小范围修复开始；B 先做缓存诊断方案。本文不改变既有决策，也不授权启动科学计算或改生产环境。
