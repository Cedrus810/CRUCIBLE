验收基准：forcefield_compiler_design.md（已全文读完）

# A3 验收：MVP-0～MVP-3 当前版本（2026-10-13）

范围：设计稿 §23 的 MVP-0、MVP-1（含 MVP-1b）、MVP-2、MVP-3。MVP-4（新电荷模型）、MVP-5（fallback /
不确定度调度）还没开始，下表标"未实现（不在本期）"的条目只是记录，不算本期验收的失败。
源码身份：NFS 上不跑 git，四个 manifest 的 `git_commit` 都是 null；四次运行记录的 `source_sha256`
覆盖同一组 104 个文件（scripts/、src/、tests/），彼此一致，和 2026-10-13 15:20 的工作区也一致（0/104 变化）。

## 1. 设计要点逐条对照

| # | 设计要点（章节） | 结论 | 依据 |
|---|---|---|---|
| 1 | 参数生成与输出格式解耦；算法不直接写 .itp / prmtop / XML（§1、§3） | 符合 | 写文件只出现在 `exporters/{gromacs,amber,openmm}.py`、`formats/amber.py`；`validation/engines.py:316` 序列化 OpenMM XML 是验证层，不是参数算法 |
| 2 | AMBER/GAFF 型 Hamiltonian，LJ 用 GAFF2，1-4 按 AMBER 约定（§2） | 符合 | `ir/schema.py:84-98`（`forcefield_family`、`combination_rule`、`coulomb_14_scale`、`vdw_14_scale`、`exclusions`、`pairs14`） |
| 3 | IR 字段（§4） | 符合 | `ir/schema.py`：Atom 有 element / atom_type / formal_charge / partial_charge / sigma / epsilon，另加 name、atomic_number、mass |
| 4 | 每个参数带 provenance（§4.1） | 符合（另有 1 处偏离，见 #12） | Bond/Angle/ProperTorsion/Improper 都有 `source` + `meta`（`ir/schema.py:36,48,62,78`）；IR 级 `provenance`（torsion、charge 两节） |
| 5 | Hessian provider 和参数算法分开，统一单位、原子顺序（§5.1） | 符合 | `hessian/base.py:37` `CartesianHessian`（nm、kJ/mol/nm²）；`hessian/xtb.py`、`hessian/gaussian.py` 只产出它 |
| 6 | provider 包括 Gaussian / ORCA / xTB / GFN-FF / CP2K（§5.1、§23 MVP-1） | 偏离（已定） | 实现了 Gaussian fchk、xtb（GFN-FF、GFN2）；ORCA `.hess` 要等真实样例（路线图 D3），CP2K 没做 |
| 7 | 核心四种方法：mSeminario、m2Seminario、DRIH、全局拟合；原始 Seminario 不对用户暴露（§6.0.1） | 符合 | `bonded/{mseminario,m2seminario,drih,global_fit}.py`；CLI 只有 `mseminario/m2seminario/drih`（`cli.py:245-246`） |
| 8 | mSeminario 是默认方法（§6.1、§7） | 偏离（小） | `derive-bonded --method` 是必填，没有默认值（`cli.py:245`）。CLI/preset 留给 C |
| 9 | 全局拟合：扣掉 H_fixed、线性最小二乘、以先验为正则（§6.4.1） | 符合（计划 §2.9 细化） | `bonded/global_fit.py:138-175`：k = s·k_prior，λ‖s−1‖²，s ≥ 0；理论推导见 `docs/theory/hessian_methods.md` §5 |
| 10 | 报告条件数、残差、几何偏差（§6.4.3） | 符合 | `GlobalFitReport`（cond、两个相对残差、x₀ 处最大力，`global_fit.py:233-244`）；G1 的 `metrics.csv` 有 FF 优化后的 kabsch_rmsd、键长/键角 MAE |
| 11 | 对称等价项共享 k（§6.4.3-3） | 符合 | `global_fit.py:144-160` 用 `term_classes` |
| 12 | provenance 要记录 H_fixed 用的是哪一版电荷 / 二面角（§6.4.2） | **偏离（未记录）** | bond/angle 的 meta 只有 prior_method、alpha、symmetrized、subtract_fixed（`global_fit.py:200-205`），没有记电荷和二面角来源。全局拟合也还没接进 CLI。归到 C（重拟合依赖） |
| 13 | 全局拟合排在 charge、torsion 之后（§6.4.2） | 未实现（属 C） | 管线还没串起来；MVP-1b 按 §23 用 GAFF2 torsion + 载入电荷做 H_fixed |
| 14 | GFN-FF 当 Hessian provider；参数叫 GFN-FF-derived effective（§8、§8.2） | 符合 | provenance 里 hessian 来源与方法分开记 |
| 15 | GFN-FF native topology 不做主路径（§9、§26） | 符合 | 没有这条代码路径 |
| 16 | flexible torsion 不从 Hessian 来（§10.1、§26） | 符合 | `bonded/assign.py:159` 只推 bond/angle；二面角只给 DRIH 当冗余坐标 |
| 17 | TorsionProvider / MLPTorsionProvider 接口（§10） | 符合 | `torsion/base.py:13`、`torsion/mlp.py:126` |
| 18 | MLP 失败要有明确兜底，不能静默留空（§6.0 最后一行） | 符合 | 兜底原因有 unsupported_element、odd_electrons、mlp_error、scan_not_converged、no_cached_scan（`torsion/mlp.py:120-206`），写进 provenance |
| 19 | 输出 Fourier 经典项（§10.3 A） | 符合 | `torsion/fourier.py`，n = 1..6，相位任意 |
| 20 | OpenMM 原生 CustomTorsionForce 模式（§10.3 B、§18.3） | 未实现（计划不做） | MVP-2 计划 §1 第 7 条 |
| 21 | ChargeProvider：载入 / RESP / AM1-BCC / 新模型（§11、§23 MVP-3） | 偏离（已定） | 载入、AM1-BCC、NAGL（`charge/`）；RESP 推迟（决策 0018）；新模型属 MVP-4 |
| 22 | 新电荷模型 learned BCC（§13–§15） | 未实现（不在本期） | MVP-4 |
| 23 | fallback graph（§16）、preset（§17） | 未实现（不在本期） | MVP-5 / C |
| 24 | GROMACS 输出 itp/top/gro（§18.1） | 符合 | `exporters/gromacs.py:42`（topol.top、molecule.itp、conf.gro） |
| 25 | AMBER 输出 mol2 + frcmod，可另生成 prmtop（§18.2） | 偏离（已定） | 只写 prmtop + inpcrd（`exporters/amber.py:215-223`）；frcmod 不做（路线图 §3.4，决策 0015） |
| 26 | OpenMM 标准力（§18.3） | 符合 | `exporters/openmm.py:38`（HarmonicBond/Angle、PeriodicTorsion、Nonbonded） |
| 27 | 跨引擎等价：分项能量 + 力（§19） | 符合 | `validation/engines.py`；G0 8 条腿 |
| 28 | 参数验证层级 Level 1–2（§20） | 部分 | Level 1（力常数、二面角剖面、电荷）与 Level 2（几何、频率）已有；Level 3–4 属 D |
| 29 | 单命令 CLI `ffcompiler ligand.sdf --preset`（§22） | 未实现（属 C） | 目前是子命令：convert、check-equivalence、derive-bonded、derive-torsions、assign-charges |
| 30 | 不做的事（§26） | 符合 | 没有重训完整 MLFF、没有重拟合 LJ、没有混用力场、没有用 Hessian 拟合 flexible torsion |

需要处理的偏离：#12（H_fixed 来源没进 provenance）与 #8（没有默认方法），两者都并入 C。
其余偏离都有决策或计划记录（#6 路线图 D3，#21 决策 0018，#25 决策 0015，#20 MVP-2 计划 §1）。

## 2. G1.8 理论文档审阅（`docs/theory/hessian_methods.md`）

对照 MVP-1 计划 §2（锁定定义）和 T15 的内容要求：

- T15 要求的 5 项都写了：mSeminario 公式与原文编号、m2Seminario C1/C2、DRIH 推导与梯度前提、全局拟合推导、单位换算。
- 和计划 §2 不一样的两处都有决策：mSeminario 的键方向和角平均（0011）；DRIH 的 B 用 Bohr 构造（0012，用户已确认）。
- 代码常数和文档一致：DRIH 伪逆阈值 1e-8（`bonded/drih.py:33`），B 用 Bohr（`drih.py:31,35`）；全局拟合 α = 1e-3、差分步长 1e-5 nm（`global_fit.py:97,101`）；线性角 175°（`bonded/internal.py:18`）；Hessian 检查阈值（`hessian/base.py:20-25`）。
- 换算核对：Ha/Bohr² → kJ/mol/nm² ×9.375829e5；频率公式里的 10²⁴；SF6 2.928942e5 ÷ 836.8 = 350.02，都对。

发现的问题（文档措辞，代码不受影响）：

| 位置 | 问题 | 改法 |
|---|---|---|
| §4.3 | 有冗余时的对角元写成 \(\sum_u k_u(\sum_i A_{ti}B_{ui})\)，少了平方。正确式是 \((\mathbf H_q)_{tt}=\sum_u k_u P_{tu}^2\)，\(\mathbf P=\mathbf A\mathbf B^\mathsf T\)。随机算例数值核对：对角元 0.86872 / 2.752792 / 3.683551，与带平方的式子逐位一致，不带平方的给出 0.14739 / 4.431767 / 4.173792。代码直接算 \(\mathbf A\mathbf H\mathbf A^\mathsf T\)，不用这个式子 | 补平方 |
| §6 第一段 | "角 kJ mol⁻¹ rad²"，应为 rad⁻² | 改笔误 |
| §5.3 "为什么是线性最小二乘" | "（几何上是平衡值 \(q^0\)，已由几何直接给定，不是拟合变量）"一句不通 | 改成"平衡值 \(q^0\) 由几何直接给定，不是拟合变量"；同时注明带下界 s ≥ 0 后是凸问题，仍是全局解 |
| §2 angle 段 | \(\mathbf u_{AB}\) 在 bond 段是 A→B，在 angle 段改成 B→A，同一个记号两种方向 | 加一句：投影力常数对 \(\mathbf u\) 的符号不变，方向只影响记号 |

结论：实质内容通过。上面 4 处改好以后，把文档开头的"待审阅（G1.8）"改成"已审阅"，G1.8 就算完成。

## 3. Gate 运行（2026-10-13）

| Gate | 运行目录 | 节点 / 耗时 | 结果 |
|---|---|---|---|
| G0 | `jobs/gate_mvp0_20261013/` | node-A / 2421 s | passed：11 固件 × 4 构型，8 条腿合计 2596 行，0 失败 |
| G1 | `jobs/gate_mvp1_20261013/` | node-A / 2610 s | passed，scope=full：G1.1–G1.7 全部为真；Sobtop 三种方法都是 12/12；全局拟合自洽 6.03e-9；1 条已知方法局限（环丙烷 m2Seminario，决策 0014） |
| G2 | `jobs/gate_mvp2_20261013/` | node-B / 2815 s | passed，scope=full：12 根键全部命中缓存（`--resume jobs/gate_mvp2_20261008_3`，缓存键重算 12/12 一致），0 兜底；sander nma dihed 2.3e-4 只报告（决策 0017） |
| G3 | `jobs/gate_mvp3_20261013/` | node-A / 1215 s | passed，scope=full：G3.1 ≤4.29e-11 e；G3.2 ≤8.3e-17 e；compare 33/33；G3.4 六个分子重拟合逐项相同；`--scan-cache jobs/gate_mvp2_20261013/scans` |

- 四个 Gate 的 pytest 都是 not external 453 passed、external 59 passed。
- 耗时比 10-05 / 10-10 长 1.5–3.5 倍：同一时段 NFS 服务器响应慢。node-B 上 `cat` 一个小文件要 0.29 s；G2 从 START 到开始读缓存花了约 28 min，10-10 只用了 56 s。
- 和旧运行对比：
  - G2 `metrics.csv` 前 20 列与 `gate_mvp2_20261010_2` 逐字相同，苯甲酸甲酯 1-2 键仍是 3.355 对 GAFF2 的 2.697，属已知缺陷，交给 B。
  - G1 `metrics.csv` 的频率列与 `gate_mvp1_20261005_4` 逐字相同；只有 benzene/gfn2 两行的几何列末位不同，差约 1e-6 pm、1e-8 deg，是 FF 极小化的数值抖动。
- 启动方式：没用 `scripts/run_logged.sh`，因为它会调用 `git -c safe.directory`。改用 `jobs/launch_a3_20261013.sh`、`jobs/launch_g3_20261013.sh`；PID 记在 `jobs/gate_mvp{0,1,2,3}.pids`。

## 4. 结论

G0～G3 在同一份源码上全部通过（正式运行，scope=full）。设计稿逐条对照没有发现与"已定决策"冲突的地方。
新发现的偏离只有 1 条（#12，全局拟合的 H_fixed 来源没进 provenance）。它和 #8 一起交给 C。
G1.8 等 4 处文档措辞修好以后收口。
