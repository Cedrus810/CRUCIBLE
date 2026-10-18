# ffcompiler 路线图（MVP-0 → MVP-5）

> 依据：`forcefield_compiler_design.md`（2026-10-04 版，已含 QMDFF §6.0 / §6.4）。
> 分工：代理写施工计划 + 每期 Gate 验收；施工方按计划施工。计划只给规格、接口、锁定的数学定义、必测用例和验收判据，不给实现代码。
>
> **2026-10-17 状态补注**：§1 状态列更新到 2026-10-17；D3（ORCA）已结案；本机有 ORCA 6.1.1（见 §0 补注）。
>
> **2026-10-10 状态补注**：分工已改为 代理 完全接手（含施工，决策 0016）。各期状态见 §1 状态列，D1–D6 的结论见 §2 状态列；§0、§4 中已过时的事实就地加了补注，原文保留作历史。当前状态以 `HANDOVER.md` 为准。

## 0. 环境事实（2026-10-04 查过）

| 东西 | 在哪 | 用途 |
|---|---|---|
| Python 环境 | mamba `openmm_dev`：openmm 8.5.2、rdkit 2025.09.6、parmed 4.3.1、openff-toolkit 0.19.0、openff-interchange 0.5.4、openff-nagl 0.6.1、numpy 2.4.3、scipy 1.17.1、torch 2.12.0、pydantic 2.13.5、networkx 3.6.1、h5py、pytest 9.1.1 | 全部 Python 工作 |
| AmberTools | `openmm_dev/bin`：antechamber、parmchk2、tleap、sqm、sander、cpptraj、resp、espgen | GAFF2 固件、AMBER 能量参照 |
| xtb 6.7.1 | `/home/user/software/xtb-dist/bin/xtb` | GFN-FF / GFN2-xTB Hessian |
| Sobtop 2026.1.16 | `/home/user/software/sobtop_2026.1.16/`（Linux 版，无执行权限，需复制到本地再 chmod） | mSeminario / m2Seminario / DRIH 的黑箱参考 |
| **没有** | GROMACS（任何环境都没有）；ORCA / Gaussian / psi4 / pyscf（`~/.local/bin/orca` 是 IDE 启动器，不是量化程序） | 见下面 D1、D3（2026-10-17：ORCA 实际有，见 §0 补注） |

不得创建、安装、升级任何环境或包。缺东西就停下报告。

> 2026-10-10 补注（上表是 2026-10-04 的快照）：GROMACS 2026.3 mixed 已有，`/home/user/gmx26.3/bin/gmx`（`FFC_GMX` 可替换，决策 0013）；Sobtop 已复制到 `tools/sobtop/` 并可运行；ParmEd 不再是依赖（决策 0015）；MVP-2 用的 MACE 模型在 `/home/user/MLP/mace/`（决策 0016）。

> 2026-10-17 补注：node-A 上有 ORCA 6.1.1：`/home/user/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2/orca`（2026-10-04 只查了 `~/.local/bin/orca` 那个 IDE 启动器）。用它生成了 .hess 样例并据此实现读取器（`experiments/20261017_orca_hess_format/`）。

## 1. 分期

| 期 | 内容 | 计划文件 | Gate | 状态（2026-10-18） |
|---|---|---|---|---|
| MVP-0 | ForceFieldIR + AMBER prmtop 导入 + OpenMM / AMBER / GROMACS 三个导出器 + 跨引擎等价测试 | `2026-10-04-ffc-mvp0-ir-exporters.md` | G0 | 通过（`jobs/gate_mvp0_20261018_2/`，当前源码） |
| MVP-1 | Hessian 容器与读取器（Gaussian fchk、xtb）+ xtb 运行器 + mSeminario / m2Seminario / DRIH + 对称化 + 全局 Hessian 拟合（MVP-1b） | `2026-10-04-ffc-mvp1-hessian-bonded.md` | G1 | 通过（`jobs/gate_mvp1_20261018_2/`）；G1.8 审阅已收口（2026-10-17）；ORCA `.hess` 读取器 2026-10-17 补上（D3） |
| MVP-2 | MLP 二面角接入（TorsionProvider）。原始材料：仓库根目录 `orb_torsion_scan_from_abfe_ibs.py`（MACE-OFF24 刚性扫描）；实际实现改为 MACE-POLAR-1-M 经 ASE 的约束松弛扫描；2026-10-08 起由代理施工（决策 0016） | `2026-10-08-ffc-mvp2-mlp-torsion.md` | G2 | 通过（`jobs/gate_mvp2_20261018_2/`，决策 0017：以 OpenMM 为准，sander 腿只报告）；methyl_benzoate 1-2 已知缺陷归 B |
| MVP-3 | 电荷 provider 框架：载入电荷、AM1-BCC（antechamber/sqm）、NAGL（openff-nagl）；RESP 往后放（决策 0018） | `2026-10-10-ffc-mvp3-charge-providers.md` | G3 | 通过（`jobs/gate_mvp3_20261018_2/`）。之后的管线串联（C1–C5）见 `2026-10-15-crucible-bc-detailed-next-steps.md`、`2026-10-17-unbuilt-modules.md` |
| MVP-4 | 新电荷模型（learned BCC） | 研究性质，单独立项（草案 `2026-10-17-mvp4-learned-bcc-proposal.md`） | G4 | 只有 provider 骨架；研究未立项 |
| MVP-5 | 自动 fallback / 不确定度调度 | 最后做 | G5 | 框架已有（`hessian/fallback.py`，CLI 不暴露）；调度研究未开始 |

MVP-0 和 MVP-1 可以并行开工：MVP-1 只依赖 MVP-0 的 T1–T3（IR schema、拓扑工具、校验）和 T5（OpenMM 导出器）。

## 2. 待用户确认的决策

| # | 决策 | 建议（计划按建议写，用户否决再改） | 状态（2026-10-10） |
|---|---|---|---|
| D1 | 本机没有 GROMACS，GROMACS 那一路怎么验证 | 用 OpenMM 的 `app.GromacsTopFile` 读我们写出的 `.top`/`.itp`，它是对 GROMACS 拓扑语义的一套独立实现；真 `gmx` 的测试写好，没有 `gmx` 时自动 skip。以后哪台机器有 GROMACS 再补跑 | 已过时：本机已有 GROMACS 2026.3，`gmx` 腿实跑（决策 0013） |
| D2 | m2Seminario 是卢天自己提出的方法，**没有发表的公式** | 先按计划里的候选定义实现，再和 Sobtop 的 m2Seminario 输出比对。全部固件都对得上，才叫 `m2seminario`；对不上就改名 `symseminario`（我们自己的对称投影法），Sobtop 的结果只作对照 | 结案：C2 与 Sobtop 12/12 一致，方法名 `m2seminario`（决策 0011） |
| D3 | DFT 级 Hessian 从哪来：openmm_dev 里没有任何 QM 程序 | MVP-1 只解析现成文件：Gaussian fchk（Sobtop 样例里有 SF6、苯甲酸甲酯）。ORCA `.hess` 读取器要等用户提供一个真实样例文件再做。是否有别的机器能跑 ORCA/Gaussian，请告知 | 结案（2026-10-17）：本机找到 ORCA 6.1.1，生成真实 .hess 样例后实现读取器，频率与 ORCA 相对差 5.7e-8 |
| D4 | 全局拟合的权重 W | 默认不加权的 Cartesian Hessian（力常数与质量无关）；质量加权作为 benchmark 选项 | 按建议执行（MVP-1 计划 §2.9，权重为 1） |
| D5 | 对称等价项是否默认取平均 | 默认取平均（按图自同构等价类），并记录类内离散度。理由：Sobtop 实测，SF6 上即使 m2Seminario/DRIH 也给出不相等的等价键角力常数（90° 角 608–692 kJ/mol/rad²），不平均会在 MD 中破坏对称性 | 按建议执行（`symmetrize=True` 默认） |
| D6 | mSeminario 的公式来源 | 施工方需要 Allen et al., *J. Chem. Theory Comput.* 2018, 14, 274 原文。请用户提供 PDF，或允许施工方在 T14 联网查阅原文（测试本身不联网） | 已解决：公式与原文编号已誊入 `docs/theory/hessian_methods.md` §2（2026-10-08） |

## 3. 设计文档里需要后续回答的问题（不阻塞 MVP-0/1）

> 2026-10-10 补注：第 1–3 条已由决策 0016 回答（AMBER 周期形式、不加谐振二面角；MLP 失败退回 GAFF2；暂用 GAFF2 电荷与 1-4 缩放）。第 3 条的"换电荷时重拟合二面角"留到 MVP-3/4。

1. **刚性二面角的谐振形式**：设计 §10.2 允许刚性二面角由 Hessian 推出（Sobtop 用 GROMACS 的 funct 2 谐振二面角）。AMBER prmtop 没有谐振二面角，只能近似成周期项。MVP-1 不推二面角；进入 MVP-2 前要定：IR 是否增加谐振二面角项、AMBER 导出如何处理（报错 / 近似 / 只允许 OpenMM 和 GROMACS）。
2. **MLP 失败兜底**：参照 QMDFF（TB 失败的 torsion 退回 Hessian 拟合），MVP-2 必须定义 MLP 推理失败或不确定度过高时的 fallback，不能静默留空。
3. **二面角与电荷的耦合**：MLP 二面角是扣掉 1-4 静电/LJ 之后的残差，绑定了训练时的电荷模型和 1-4 缩放。MVP-2 要查清 MLP 的训练基线；MVP-4 换电荷时要评估是否需要重拟合二面角。
4. **mol2 + frcmod 路线**：Hessian 推出的参数是按原子实例给的，frcmod 按原子类型存。AMBER 导出以 prmtop 为主；frcmod 路线不做（YAGNI），除非用户需要。

## 3a. MVP-2 对接说明（给对接人员）

> **2026-10-08 已作废（施工方部分）**：用户决定由代理完全接手，MVP-2 由代理写计划并施工（决策 0016）。下面的材料说明仍有效。

~~MVP-2 由对接人员负责，代理不写施工计划。~~下面是现有材料的情况：

- **代码**：仓库根目录 `orb_torsion_scan_from_abfe_ibs.py`，取自 ABFE_IBS commit 56c2a92 的 `abfe_core.py`（`OrbScanner`、`OrbMMHybridFactory.get_rotatable_torsions_rdkit`）。
- **模型**：MACE-OFF24 medium，放在 `/home/user/MLP/mace/`，经 openmm-ml 的 `MLPotential("mace", modelPath=...)` 加载。工程里一律用这个路径，**不要**沿用代码里写死的 `~/.cache/mace/MACE-OFF24_medium.model`（那是另一个账号家目录下的缓存，换账号就找不到）。
- **环境**：openmm_dev 里已有 openmm-ml、mace 0.3.16、e3nn 0.4.4、openmmtorch，不需要装任何东西。
- **这不是一个预测二面角参数的模型**：它用 MACE-OFF24 对单个二面角做 1D 扫描，返回原始能量曲线 E(φ)（kJ/mol），没有扣掉 MM 部分，也没有拟合成参数。ABFE_IBS 那边是把修正项拟合成 Fourier 级数 offset + Σ cₙcos(nφ) + sₙsin(nφ)，用 `CustomTorsionForce` 加回 System，不写回拓扑。接进 IR 需要：
  1. 扣掉 MM 部分：在同一批构型上算 \(E_\text{MLP}(\phi)-E_\text{MM}^{\setminus\text{torsion}}(\phi)\)，其中 MM 部分要用**最终的**电荷、LJ、1-4 缩放和 bond/angle（见上面 §3.3）；
  2. 拟合成 IR 的 \(\sum k_n[1+\cos(n\phi-\delta_n)]\) 项（设计 §10.3 A）；
  3. 定义失败兜底（上面 §3.2）。
- **扫描是刚性的**（RDKit `SetDihedralDeg` 直接转，不做约束优化）：能垒容易偏高，原子可能撞在一起。建议改成约束松弛扫描。
- **可转键的选择规则** `[!#1;!D1]-[!#1;!D1]` 没有排除环内键（RDKit 转环内键的二面角会直接报错），而且对同一根键的每个四原子组合各扫一遍，有重复。应改为每根非环可转键扫一次，再拟合这根键上的全部二面角项。
- `device` 默认写死 `"cuda"`；本机 GPU 是 RTX 2080 Ti（11 GB）。

## 4. 全期通用约束

- 仓库：`/home/user/forcefield_compiler`，包名 `ffcompiler`（src 布局）。NFS 上 git 写对象会出权限错，git 目录放本地盘 `/home/user/.local/share/git-sep/ffcompiler.git`（与 Amber-X 相同做法）。
  - **2026-10-10 补注（上一条已过时）**：工作区 2026-10-08 改名为 `/home/user/CRUCIBLE`，git 元数据在工作区 `.git/`；agent 不在这里执行 git，提交写成 `.commit_drafts/` 草稿由用户提交（见 `HANDOVER.md`「本地 Git」）。
- 日志：每行带 `[YYYY-mm-dd HH:MM:SS +elapsed]`，开头 START（主机、argv、git commit），结尾 END（退出码、总 wall time）。
- 正式运行（Gate、benchmark）输出到 `jobs/<name>_<YYYYMMDD>/`，带 `manifest.json` 和 `acceptance.json`。产物放共享盘，不放 `/tmp`。
- 单元测试不联网；单元测试里不调用 xtb、sander、Sobtop，这些放到标了 marker 的集成测试里（`-m external`）。
- 正确性测试按核数分线程并行跑；只有性能基准才需要独占机器。
