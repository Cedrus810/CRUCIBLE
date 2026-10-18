验收基准：forcefield_compiler_design.md（已全文读完）

# 结案验收：当前版本（2026-10-18）

范围：10-13 A3 验收（`ACCEPTANCE_20261013_A3.md`）之后的全部新增，即 10-17 的 ORCA 读取器、C1–C5 管线合同、MVP-5 fallback 框架、MVP-4 骨架（当晚验收修复见 `docs/superpowers/plans/2026-10-17-unbuilt-modules.md` §0），10-17/18 的 B1–B3 诊断与实验，以及 10-18 的可选顺序重拟合；在 A3 的逐条对照表基础上更新。

源码身份：最新一轮 Gate `jobs/gate_mvp{0,1,2,3}_20261018_2/` 的四份 manifest 共用同一组 `source_sha256`，共 119 个文件，与验收时的工作区 0/119 不同。`git_commit` 为 null（决策 0021：NFS 上不探测 git）。提交还在 `.commit_drafts/` 里等用户在实际机器上执行，执行后把哈希补记到本报告。

补记（2026-10-18）：草稿已全部执行完毕（用户委托 代理 代执行；草稿中的共同作者尾注按仓库规矩剥除，`仓库代理规程` 与 `.superpowers/` 等本地工具文件按用户指示不入库，另有一笔 gitignore 把它们挡在库外）。六份作废草稿（`20261017_1240/_1500/_1540/_1620/_1650/_1720`）未执行；`docs/decisions/0011` 的一处工作区改动（Allen 2018 公式对照补记）无草稿认领，留在本地未提交。本补记提交之前的全部提交：

- ed5fc25 chore: keep local agent tooling and commit drafts out of the repository
- 3bc58ff docs(handover): workspace renamed to CRUCIBLE; update local git notes
- 885aa27 feat(torsion): Fourier fit and AMBER periodic-term conversion (MVP-2 T24)
- 3230a95 feat(torsion): rotatable-bond perception and dihedral symmetry classes (MVP-2 T25)
- 4cf80ab feat(torsion): ASE constrained relaxed dihedral scan and MM helpers (MVP-2 T26)
- 36d24c3 feat(torsion): MACE-POLAR-1-M reference calculator via ASE (MVP-2 T27)
- c3cc1a0 feat(torsion): MLPTorsionProvider with MM subtraction, fit and GAFF2 fallback (MVP-2 T28)
- f6616e2 feat(cli): derive-torsions subcommand (MVP-2 T29)
- 25981b6 feat(gate): Gate G2 script for MLP torsions (MVP-2 T30)
- e8000b3 fix(torsion): model provenance, cache identity, cache-only gate fitting (MVP-2 review)
- 5c19450 docs(mvp2): record first Gate G2 run (failed on G2.3) and nma diagnosis
- 2f31d34 docs: sync MVP-2 status, annotate roadmap, retire design v0
- b28e2ca gate_mvp2: judge G2.3 with OpenMM as reference, sander legs report-only
- b1f4ca3 docs: MVP-3 charge-provider plan and decision 0018
- 67d06f1 feat(charge): MVP-3 charge providers with torsion refit (Gate G3)
- ba63436 fix(gates,charge,torsion): A1/A2 gate scope and provenance close-out
- 22ae02e Refit keeps original fallback bonds; strict resume/scan-cache paths (decision 0020)
- c40d8ee A3 acceptance: G0-G3 re-run on one source tree, all passed
- 0648979 Close G1.8 review: four wording fixes in hessian_methods.md
- 4a805f8 Add ORCA reader, C1-C5 pipeline contract, MVP-5 fallback, MVP-4 skeleton; fix on acceptance
- 987f086 Docs: sync status to 2026-10-17; slim HANDOVER, move daily logs to history
- 113e7a5 fix(runtime): remove automatic Git probes from run recording
- ba762a3 B1-B3: diagnose MLP torsion profile defect; sequential refit mechanism experiment
- 2d0bd7f Re-run full pytest and Gates G0-G3 on the current source: all pass
- 11b188b B3 order-dependence check; C2 chain test now includes the torsion step
- 36f0f40 feat(torsion): add optional bounded sequential refinement and replay protocol
- 2db460b Re-run Gates G0-G3 after adding optional sequential refinement: all pass
- 62837fc Closure acceptance of the current version against the design document

验收方式：设计要点逐条对照代码（只读），最后才看测试和 Gate。本报告不含新的计算；引用的运行结果见第 3 节。

## 1. 设计要点逐条对照

结论：**符合**、**偏离（已定）**＝有决策或计划记录、**偏离（待定）**＝要用户拍板、**部分**、**未实现**。

| # | 设计要点（章节） | 结论 | 依据 |
|---|---|---|---|
| 1 | 参数生成与输出格式解耦；算法不直接写 .itp / prmtop / XML（§1、§3） | 符合 | 写引擎文件只在 `exporters/{gromacs,amber,openmm}.py`、`formats/amber.py`。其余写文件的地方是外部程序的输入（`hessian/xtb.py`、`charge/am1bcc.py`）、IR 与坐标 JSON（`ir/io.py`、`ir/coords.py`）和验证层 |
| 2 | 完整管线从分子出发：化学感知 → atom typing（GAFF/GAFF2）（§1、§3、§21 `molecule/`、`typing/`） | **未实现** | 没有 `molecule/`、`typing/` 模块。输入必须是已经参数化好的 AMBER prmtop/inpcrd + SDF（`importers/amber_prmtop.py`）。任意 SDF 的 GAFF2 初始化不在本期（10-15 计划 §5 C2 明确排除） |
| 3 | AMBER/GAFF 型 Hamiltonian；LJ 用 GAFF2；1-4 按 AMBER 约定（§2） | 符合 | `ir/schema.py:84-98` |
| 4 | IR 字段（§4） | 符合 | `ir/schema.py:15-98` |
| 5 | 每个参数带 provenance（§4.1） | 符合（比 A3 时更完整） | 各项都有 `source` + `meta`。IR 级记录新增：`bonded.h_fixed_dependency`（`bonded/global_fit.py:245`）、`torsion.mm_dependency` / `fit_protocol`（`torsion/mlp.py:299-306`）、`charge.previous`（`charge/apply.py:36`）；坐标来源记在 sidecar（`ir/coords.py`） |
| 6 | Hessian provider 与参数算法分开，统一单位与原子顺序（§5.1） | 符合 | `hessian/base.py` 的 `CartesianHessian`；各 provider 只产出它 |
| 7 | provider 包括 Gaussian / ORCA / xTB / GFN-FF / CP2K（§5.1） | 部分 | 有 Gaussian fchk、xtb（GFN-FF、GFN2）、**ORCA .hess（新增，`hessian/orca.py:108`）**。ORCA 读取器按 ORCA 6.1.1 的真实输出实现，频率与 ORCA 自己的结果相对差 5.7e-8，但只核对过这一个样例（闭壳层、每块 5 列）。CP2K 没做 |
| 8 | 核心四种方法；原始 Seminario 不对用户暴露（§6.0.1） | 符合 | `bonded/{mseminario,m2seminario,drih,global_fit}.py`。CLI：`derive-bonded --method {mseminario,m2seminario,drih}`，另有 `global-fit` 子命令（`cli.py:533`） |
| 9 | mSeminario 是默认方法（§6.1、§7） | **偏离（待定）** | `derive-bonded --method` 仍是必填参数，没有默认值（`cli.py:498`）。A3 #8 至今未定 |
| 10 | 全局拟合：扣 H_fixed、线性最小二乘、以先验为正则（§6.4.1） | 符合 | `bonded/global_fit.py`（与 A3 时相同） |
| 11 | 报告条件数、残差、x₀ 处的力（§6.4.3） | 符合 | `GlobalFitReport`；`global-fit` 把它写进 `global_fit_report.json` |
| 12 | 对称等价项共享 k（§6.4.3-3） | 符合 | `term_classes` |
| 13 | provenance 记录 H_fixed 用的是哪一版电荷 / 二面角（§6.4.2） | **符合（A3 #12 已关闭）** | 逐项快照（不用 hash）加 Hessian 身份（`bonded/global_fit.py:245-256`）；拟合后输入被改过，`check-staleness` 报过期（`bonded/staleness.py:119`） |
| 14 | 全局拟合排在 charge / torsion 之后，换模型要重拟合（§6.4.2） | 部分 | 依赖记录和过期检查已有（C3/C4），不会静默沿用旧拟合。但 bond/angle 与 torsion 的双向依赖怎样交替更新、什么时候算"组合完成"，**还没定**（10-15 计划 C4） |
| 15 | GFN-FF 当 Hessian provider；参数叫 GFN-FF-derived effective（§8、§8.2） | 符合 | provenance 里 Hessian 来源与方法分开记 |
| 16 | GFN-FF native topology 不做主路径（§9、§26） | 符合 | 没有这条代码路径 |
| 17 | flexible torsion 不从 Hessian 来（§10.1、§26） | 符合 | `bonded/assign.py` 只推 bond/angle |
| 18 | TorsionProvider / MLPTorsionProvider（§10） | 符合 | `torsion/base.py`、`torsion/mlp.py:160` |
| 19 | MLP 失败要有明确兜底（§6.0） | 符合 | 逐键兜底原因；重拟合不升级原兜底键（决策 0020） |
| 20 | 输出 Fourier 经典项（§10.3 A） | 符合（新增可选项） | `torsion/fourier.py`；**可选顺序重拟合** `--refinement-rounds {0..3}`，默认 0 仍是单遍（`torsion/mlp.py:264`、`cli.py:514`）；换电荷时按记录的协议重放（`torsion/mlp.py:316`） |
| 21 | OpenMM 原生 CustomTorsionForce 模式（§10.3 B、§18.3） | 未实现（计划不做） | MVP-2 计划 §1 第 7 条 |
| 22 | ChargeProvider：载入 / RESP / AM1-BCC / 新模型（§11） | 偏离（已定） | 载入、AM1-BCC、NAGL；RESP 推迟（决策 0018）；新模型见 #23 |
| 23 | 新电荷模型 learned BCC：q = q⁽⁰⁾ + Σ Δq，按键反对称、守恒（§13） | 部分（只有骨架） | `charge/learned_bcc.py:62`：形式和符号约定已锁定，工件缺失就报错，base 电荷来源必须显式给出。**没有训练好的模型，没有接进 CLI**；立项草案 S1–S6 都还没定 |
| 24 | 多目标训练、charge 与 LJ 的边界（§14、§15） | 未实现（下一阶段） | MVP-4 立项草案 |
| 25 | fallback graph（§16） | 部分（只有框架） | `hessian/fallback.py:68,125`：GFN-FF → GFN2 → DFT 文件，用 `check_hessian` 作质量门。§16 第 1 级"可信库参数"没有接；CLI 不暴露；不确定度调度研究没开始 |
| 26 | preset（§17）、单命令 CLI `ffcompiler ligand.sdf --preset`（§22） | **未实现** | 现在只有子命令。preset 依赖 #2 和 #14 |
| 27 | GROMACS 输出 itp/top/gro（§18.1） | 符合 | `exporters/gromacs.py:42` |
| 28 | AMBER 输出 mol2 + frcmod，可另生成 prmtop（§18.2） | 偏离（已定） | 只写 prmtop + inpcrd（决策 0015；frcmod 按 YAGNI 不做） |
| 29 | OpenMM 标准力（§18.3 portable） | 符合 | `exporters/openmm.py:38-59` |
| 30 | 跨引擎等价：分项能量 + 力（§19） | 符合 | `validation/engines.py`；G0 共 8 条腿 |
| 31 | 参数验证层级（§20） | 部分 | Level 1（力常数、剖面、电荷）与 Level 2（几何、频率）有；Level 3–4 没有（属 D）。B1–B3 是工程/方法诊断，**全部是样本内结果** |
| 32 | MVP 顺序（§23） | 符合 | MVP-0～3 与 MVP-1b 完成且有 Gate；MVP-4 只有骨架、MVP-5 只有框架，符合"最后才做" |
| 33 | 当前明确不做的事（§26） | 符合 | 没有重训完整 MLFF、没有重拟合 LJ、没有混用力场、没有用 Hessian 拟合 flexible torsion、没有做 residue/蛋白模板 |
| 34 | 第四种方法叫 QMDFF-inspired global Hessian fit，不写成"实现了 QMDFF"（§6.0） | 符合 | CLI 帮助写的是"QMDFF 式全局 Hessian 拟合"；代码注释同样 |

## 2. 结案时的未完成项与已知缺陷

**已知缺陷（当前版本按默认设置就会遇到）：**
1. **苯甲酸甲酯 1-2 键**：默认单遍下剖面 RMSE 3.36 kJ/mol，比 GAFF2（2.70）差。B2 已用能量账确认原因是单遍拟合的跨键耦合（`jobs/torsion_diagnosis_20261017/report.md`）。`--refinement-rounds 3` 降到 0.44，但只是样本内结果；所有更新顺序都收敛到同一不动点，第 3 轮后顺序带来的差 ≤ 1e-4 kJ/mol（`jobs/torsion_b3_order_20261018/report.md`）。默认值没有改。
2. **nma 1-3（酰胺 C–N）扫描**：正反扫滞后最大 57.4 kJ/mol，14/36 帧超过 1 kJ/mol，是能垒区的扫描分支问题。它的"拟合改善"不算科学收益；是否重扫没有定。没有设自动兜底阈值（用户决定先看分布）。
3. **单一四元组近似**：同一根键的能量全部放在一个四元组上（MVP-2 §7 第 3 条）。B2 已排除它是苯甲酸甲酯 1-2 的主因，但它在什么范围内适用没有定。
4. **ORCA 读取器**只用 ORCA 6.1.1 的一个闭壳层样例核对过。

**要用户拍板的：**
5. `derive-bonded` 是否默认 mseminario（上表 #9）。
6. 版本号（现在是 0.2.0，代码已超出 MVP-3）。
7. 顺序重拟合要不要改成默认，取几轮。依据要等 B4。
8. MVP-2 §7 剩下的几条：`spin=1` 的语义、几何变形（键长变化）阈值。

**没做的验证：**
9. **C5 真实组合验收**：串联测试（`tests/unit/test_pipeline_chain.py::test_local_route_chain_with_torsions`）用的是 `--reference mm` 替身，**真实 MACE 参考势的完整组合验收完全没做**。
10. **B4 留出验证**：B1–B3 都在拟合用的帧上评估；要有独立构象或扫描，才能判断能否推广。
11. **D 独立科学 benchmark**（§20 Level 3–4）。

**下一阶段（不在本期）：** MVP-4 learned BCC 的数据、模型与训练（S1–S6）；MVP-5 调度的校准；§1 的 typing 入口与 §17 preset；global-fit 与 torsion 的交替策略（#14）。

**10-13 以来已关闭的 A3 偏离：** #6 ORCA（已实现，见 #7）；#12 H_fixed provenance（见 #13）；#13 管线串联（局部路线的单元级完成，全局路线的策略待定，见 #14）。

## 3. 验证证据（最后才看）

| 项 | 运行 | 结果 |
|---|---|---|
| Gate G0–G3（当前源码） | `jobs/gate_mvp{0,1,2,3}_20261018_2/`（node-A，10-18 08:50–09:15） | 全部 passed，数值与 10-13、10-18 凌晨两轮一致；每个 Gate 的 pytest 都是 524 + 59，0 失败 |
| 单元测试 | 同上（Gate 内） | 非 external 524 passed；external 59 passed（含 MACE CLI 用例） |
| 顺序重拟合数值一致性 | `experiments/20261018_refinement_verify/` | 3 轮与 B3 实验逐位相同；默认 0 轮与改动前的单遍逐位相同 |
| 10-17 施工的验收修复 | `2026-10-17-unbuilt-modules.md` §0 | 源码 bug 2 个、代码问题 7 个、测试错误 8 个，均已修 |
| ORCA 格式核对 | `experiments/20261017_orca_hess_format/` | 频率与 ORCA 相对差 5.7e-8 |
| B1–B3 | `jobs/torsion_diagnosis_20261017/`、`jobs/torsion_b3_sequential_20261017/`、`jobs/torsion_b3_order_20261018/` | 见第 2 节第 1 条 |

Gate 走的是默认单遍；三轮顺序模式没有对应的 Gate 判据，只由单元测试和数值一致性核对覆盖。

## 4. 结论

当前版本可以结案，范围如下：MVP-0～MVP-3（含 MVP-1b）有 Gate 证据；C1–C4 的工程合同（IR + 坐标传递、依赖记录、过期检查、`global-fit` CLI）和 C2 局部路线的单元级串联已完成；ORCA 读取器、MVP-5 fallback 框架、MVP-4 骨架、可选顺序重拟合已交付，并分别注明了边界。

第 2 节的已知缺陷（1–4）、待拍板项（5–8）和没做的验证（9–11），作为明确的未完成项留给下一阶段。本报告不把它们算作"通过"。结案还差两步：用户执行积压的提交草稿，并把提交哈希补记到本报告。
