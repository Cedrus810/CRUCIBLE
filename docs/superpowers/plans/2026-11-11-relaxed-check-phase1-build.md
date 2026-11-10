# 受约束松弛检查 阶段 1 施工计划（2026-11-11）

状态：**已执行（2026-11-11）**，报告见 `jobs/torsion_relaxed_check_20261111/report.md`。
执行偏差一处：§2 第 1 条"X 的松弛直接调 `_relax`"改为脚本内 `MultiDihedral` 约束——ASE
FixInternals 在 dmso 的 4-同轴约束上不收敛（详见报告 §2），定义、判据、验收门均未动。
阶段 2 仍待用户批准。规格与判据以
`docs/superpowers/plans/2026-11-01-relaxed-check-plan.md` 为准（其 §2 的量定义、§3 的判据
**照抄锁定，不得看结果后再改**）；本文只补充代码接口、验收数值、运行步骤与交付物。
**阶段 2（新增分子、算 MACE 扫描）不在本计划内，需另行批准。**

## 必读（开工前全部读完）

1. `docs/superpowers/plans/2026-11-01-relaxed-check-plan.md`：全文，尤其是 §2 定义、§3 判据与停止条件、§6 明确不做。
2. `docs/theory/innovation_screen_20261029.md` §6（研究门槛）、§8（2026-11-01 冒烟结果与用户结论原文，含两条边界）。
3. `jobs/torsion_force_bound_20261101/report.md`：X 量（relax）的已有实现与苯甲酸甲酯对照值；用户已否掉欧氏力下界，本轮不做力。
4. `jobs/torsion_diagnosis_20261017/report.md`：B2 能量账；本计划 §5 验收数值的来源。
5. `scripts/torsion_relax_energy.py`：`_relax`（X 的松弛照用）、`independent_quartets`。
6. `scripts/torsion_force_bound.py`：`_load_fixture`、`quartets_on_bond`。
7. `scripts/diagnose_torsion_profiles.py`：only/all IR 的构造方式（`diagnose_molecule` 内）、`_bestshift_rmse` 口径。
8. `HANDOVER.md`「环境」节与「本地 Git」节。

## 0. 分工与边界

- 施工方：写脚本 → 冒烟 → 全量运行 → 写 `report.md` → 更新 HANDOVER → 写提交草稿与台账。
  **不运行任何 git 命令**（NFS 工作区，决策 0021 与仓库规程）；**不动环境**（只
  `source scripts/env.sh` 进 openmm_dev，不 install/upgrade 任何包）；不加载 MACE、不调 xtb、
  不重扫（`MLPTorsionProvider(..., cache_only=True)`）。
- 只读输入：`jobs/gate_mvp2_20261013/`（`scans/`、`metrics.csv`、`ir/`）与 `tests/fixtures/gaff2/`。
- 新产物只有三处：`scripts/torsion_relaxed_check.py`（新）、`jobs/torsion_relaxed_check_<执行日期>/`
  （新）、台账 `.superpowers/sdd/2026-11-11-relaxed-check-phase1/progress.md`（新）。`src/` 一行不动。
- 任何验收数值（§5）不符、任何松弛不收敛、任何约束漂移超界 → **停**，把现象写清楚等用户，不自行解释过去。

## 1. 规格锁定（摘自 2026-11-01 计划 §2/§3，原文口径）

对键 r、扫描帧 s（参考几何 x_s、参考能量 E_ref(s)），IR ∈ {only, all}（同 B2：only 只把 r 换成
拟合项，all 为 G2 最终 IR）：

- e_fix(s) = E_IR(x_s) − E_ref(s)
- e_rel(s) = min{E_IR(y) : φ_q(y) = φ_q(x_s)} − E_ref(s)。只约束扫描四元组 q，**逐帧从 x_s 出发**
  （不做顺序扫描，避免分支），LBFGS fmax 0.001 eV/Å；e_rel 的松弛用**完整 IR**（保留 r 上的拟合
  扭转项），即 OpenFF 式重新松弛。r 上非扫描四元组可以动——那是机制本身，不是错误。
- relax(s) = E_bg(x_s) − min{E_bg(y) : Φ_r(y) = Φ_r(x_s)}，Φ_r = r 上全部四元组（取线性无关子集，
  被去掉的由恒等式随之固定并逐帧核对）；E_bg = 该 IR 去掉 r 上的 proper（only → GAFF2 背景 =
  `strip_bond_torsions(base, bond)`，all → 最终 IR 背景 = `strip_bond_torsions(redo, bond)`）。与 θ_r 无关。
- e_rel_pred(s) = e_fix(s) − [relax(s) − mean(relax)]（mean 对该键 36 帧取）。
- P0 = RMSE_bs(e_fix)（G2 现有指标）；P1 = 扫描最大滞后；Y = RMSE_bs(e_rel) − P0；
  X = RMSE_bs(e_rel_pred) − P0。RMSE_bs = 最佳常数平移 RMSE（kJ/mol），参考值都是同一条 MACE
  扫描 E_ref。
- σ = 0.41 kJ/mol（B2 中扫描质量正常的 11 根键最大滞后；只适用本项目扫描参数）。
- 逐帧机制：报告 RMS(e_rel − e_rel_pred) 与 RMS(e_rel − e_fix) 之比，只报告不设阈值。
- 判据只用 11 根键（nma 1-3 滞后 57 kJ/mol，按用户 2026-10-17 决定单列，**照算但不进混淆表**）。
- 停止条件（任一成立即停，不把 X 当诊断指标推进）：没有任何键 Y > σ；或 X 的分类不比 P0、P1 好；
  或只在已知的苯甲酸甲酯 1-2 上成立。

## 2. 代码任务：`scripts/torsion_relaxed_check.py`（新）

一个独立实验脚本，风格照 `torsion_relax_energy.py`（docstring 写清量定义与命令；`ROOT` 加
sys.path；manifest 记身份）。要点：

1. **复用**：`torsion_force_bound._load_fixture`、`quartets_on_bond`；
   `torsion_relax_energy._relax`、`independent_quartets`。X 的松弛直接调 `_relax`
   （E_bg、线性无关子集、fmax、max_steps 同默认）；**不需要**再重复 11-01 的鞍点扰动检查。
2. **新增** e_rel 的松弛函数（只约束 q）：ASE `Atoms(numbers, positions=x/ANGSTROM_NM)` +
   `OpenMMCalculator(ir, validate_ir=False)`；**e0 在 set_constraint 之前**对同一计算器求值
   （保证 e_fix 与 e_rel 同一能量函数，差值干净）；target 取 x 处的 φ_q；LBFGS(fmax=0.001,
   steps≤5000)。返回 e0、e1、收敛标志、φ_q 漂移（应为 ~0）、Φ_r∖{q} 的最大漂移（数据列）、
   Kabsch 对齐后最大位移。
3. **IR 构造**：`base, pos, mol = _load_fixture(name)`；
   `redo, report = MLPTorsionProvider(..., cache_only=True).parameterize(base, pos, mol)`；
   only 的构造照 `diagnose_torsion_profiles.py`（base 去掉 r 的 GAFF2 proper，加上 redo 里 r 的
   拟合项）。只处理 `status == "mlp"` 的转子（12 根键）。断言 scan 四元组 ∈ quartets_on_bond
   （允许反序）。
4. **输出**：
   - `frames/<分子>_<键>.csv`：逐帧一行。必含 frame、phi_deg、converged、hysteresis、E_ref，
     以及每个变体（only/all）的 e_fix、e_rel、relax_drop(=e_fix−e_rel)、relax_conv、dphi_q_rad、
     dphi_others_rad、max_disp_nm、bg_relax、bg_conv、bg_dphi_rad、bg_max_disp_nm、
     e_fix_numpy（NumPy `evaluate` 交叉核对）。列名可调，字段必须齐。
   - `summary.csv`：每（分子、键、变体）一行，含 n_frames、scan_all_converged、p1_max_hysteresis、
     p0_rmse_bs、rmse_bs_e_rel、Y、rmse_bs_e_rel_pred、X、rms_e_rel_minus_pred、
     rms_e_rel_minus_fix、rms 比值、bg_relax 的 median/min/max/range、relax_drop 的 min/median、
     收敛标志、自检列（见 §5）。
   - `manifest.json`：task、argv、host、python、wall_s、fmax、max_steps、`sigma_kJ_mol: 0.41`、
     g2_run、molecules、inputs_sha256（metrics.csv、ir/*.json、scans/*.npz）、source_sha256
     （本脚本 + src/ffcompiler 全部 .py；照 diagnose_tvp 的写法，只作身份不作复用判断）。
5. **续跑**：frames CSV 已存在的键跳过（按文件存在判断，不用 hash）；summary/manifest 从
   CSV 重建，覆盖写。
6. `RMSE_bs` 口径照 `diagnose_torsion_profiles.py::_bestshift_rmse`。

## 3. 分子与规模

G2 的 6 个分子 12 根键：acetate(0-1)、dmso(0-1, 1-2)、ethanol(0-1, 1-2)、
methyl_benzoate(0-1, 1-2, 2-4)、methylammonium(0-1)、nma(0-1, 1-3, 3-4)。
每键 36 帧 × 4 组松弛（X、Y 各 × only/all）= 1728 次 LBFGS。预计每分子 ≤ 10 min（冒烟实测后修正）。

## 4. 运行步骤

1. **冒烟（先做，单分子）**：确认负载（`uptime`；40 物理核与他人共用）后：
   ```bash
   cd /home/ruigengji/CRUCIBLE && source scripts/env.sh
   OMP_NUM_THREADS=1 python scripts/torsion_relaxed_check.py \
       --g2-run jobs/gate_mvp2_20261013 --outdir jobs/torsion_relaxed_check_<日期> \
       --molecules methyl_benzoate
   ```
   stdout 落 `jobs/torsion_relaxed_check_<日期>/run.log`，PID 记 `jobs/<…>.pids`（沿用现有惯例）。
   冒烟完立即对 §5 验收数值；**对不上先停**。实测耗时回填全量估计。
2. **全量**：其余 5 个分子（acetate dmso ethanol methylammonium nma）各一个进程，参数同上，
   同一个 outdir（methyl_benzoate 已有的键自动跳过）。6 进程 × 1 线程远低于机器容量。
3. 汇总 `summary.csv` → 写 `report.md`（§6）。

## 5. 验收数值（施工后逐条对照，不符即停）

**(a) P0 与 B2 bestshift 表逐键一致**（容差 1e-3 kJ/mol；同机同口径应远小于此）：

| 分子 键 | P0 only | P0 all | | 分子 键 | P0 only | P0 all |
|---|---|---|---|---|---|---|
| acetate 0-1 | 0.0074 | 0.0074 | | methyl_benzoate 1-2 | 0.60 | 2.70 |
| dmso 0-1 | 0.023 | 0.099 | | methyl_benzoate 2-4 | 0.037 | 0.045 |
| dmso 1-2 | 0.025 | 0.062 | | methylammonium 0-1 | 0.0061 | 0.0061 |
| ethanol 0-1 | 0.022 | 0.038 | | nma 0-1 | 0.030 | 0.031 |
| ethanol 1-2 | 0.017 | 0.103 | | nma 1-3 | 0.68 | 0.95 |
| methyl_benzoate 0-1 | 0.060 | 0.060 | | nma 3-4 | 0.021 | 0.061 |

**(b) only 变体 bg_relax 与 11-01 relax 报告一致**（methyl_benzoate，同一流程同一量；容差 0.01）：
0-1 中位 12.79（12.58–13.65）、1-2 中位 10.10（8.73–14.75）、2-4 中位 11.37（10.66–12.81）。

**(c) 自检全部通过**：所有 LBFGS 收敛；relax_drop 与 bg_relax 均 ≥ −1e-9（松弛只降不升）；
dphi_q、bg_dphi ≤ 1e-6 rad（11-01 实测 1.4e-9，超 1e-6 说明约束没固定住）；
e_fix 的 NumPy vs OpenMM 差 ≤ 1e-5（B2 实测 2.7e-6）；扫描帧 36/36 收敛（照抄 B2 状态）。

## 6. report.md 必须包含（写到 `jobs/torsion_relaxed_check_<日期>/report.md`）

1. 命令、机器、墙钟、输入身份（manifest 摘要）、与 §5 验收数值的逐条对照结果。
2. 逐键表：P0、X、Y（only、all 各一组）、P1、各自 >σ 标志；nma 1-3 单列并注明不进判据。
3. 混淆表（11 键，不算单一分数）：X>σ、P0>σ、P1>σ 各自对"Y>σ"分类的 2×2。
4. 逐帧机制：每键 RMS(e_rel−e_rel_pred) / RMS(e_rel−e_fix) 比值，只报告。
5. 停止条件三条逐条判定，结论一句话写明 X 是否值得推进（推进也只意味着"可考虑阶段 2"，
   阶段 2 仍需用户批准）。
6. 不写"已证明/首个"；结果只在样本内、机制验证层面（11 键，Spearman 标准误 ≈ 0.32）。

## 7. 明确不做（照 2026-11-01 计划 §6）

不改 G2 拟合流程；不实现 OpenFF 式重新松弛拟合；判据不进 pipeline、不加门；不做自由能、MD、
独立 QM 核查；不做阶段 2；不做 11-01 已否掉的欧氏力下界。

## 8. 收尾

- `HANDOVER.md` 顶部追加本轮记录（做法、验收对照、结果摘要、PID/脚本路径、产物路径）。
- 台账 `.superpowers/sdd/2026-11-11-relaxed-check-phase1/progress.md`。
- 提交草稿 `.commit_drafts/20261111_relaxed-check-phase1.md`：格式照 `20261017_2300_b1b2b3-torsion-diagnosis.md`
  （commit message + 改动文件 + 命令）。入库：脚本、`jobs/.../report.md`、`manifest.json`、
  HANDOVER、2026-11-01 计划的状态行；`frames/`、`summary.csv` 等已被 `.gitignore` 的 jobs 产物只留磁盘。
  运行记录里的真实路径先脱敏再写进任何入库文档。
