# 决策 0016：MVP-2（MLP 二面角）的前置决定

- 2026-10-08 用户决定（对应 roadmap §3 的三个问题与 §3a）：
  1. **二面角形式**：用户原话"我是按照 amber 的模式写的"、"hessian 算不了二面角"。
     二面角一律用 IR 已有的 AMBER 周期形式 \(\sum k_n[1+\cos(n\phi-\delta_n)]\)；
     不新增谐振二面角（GROMACS funct 2）项，Hessian 不推二面角。因此 AMBER 导出不受影响，
     roadmap §3.1 关于"谐振二面角在 AMBER 里只能近似"的问题不再存在。
  2. **MLP 失败兜底**：退回 GAFF2 的二面角参数，并在报告里逐条记录哪根键退回、原因是什么。
  3. **电荷耦合**：暂时用现有的（GAFF2 prmtop）电荷与 1-4 缩放扣除 MM 部分；换电荷（MVP-3/4）时再评估。
  4. **模型路径**：用户确认用 `/home/user/MLP/mace/`，不用 `~/.cache/mace/`。
     路径是**大写** `MLP`。2026-10-08 代理把 roadmap §3a 里的小写 `mlp` 当真，误建了 `/home/user/mlp/mace/`
     并拷了一份模型；用户指出后已删除（删前 `cmp` 与 `MLP/mace/` 里的同名文件一致，冒烟结果不受影响）。
     `/home/user/MLP/mace/` 里另有 MACE-POLAR-1-S/M/L、MACE-omol-0 等模型。
- 模型覆盖的元素（2026-10-08 读模型文件 `atomic_numbers` 核实）：H C N O F P S Cl Br I，r_max = 6 Å。
  模型没有总电荷输入。
- **修订（2026-10-08，同日）**：改用 `/home/user/MLP/mace/MACE-POLAR-1-M.model`，**直接经 ASE**
  （mace 的 `MACECalculator(model_type="PolarMACE")`），不经 openmm-ml。用户原话"对的直接ase"。
  - 原因：POLAR 有总电荷输入（Fukui 均衡），离子可以照常扫；openmm-ml 1.6 不传 `fermi_level` /
    `external_field`，加载后求能报 `'fermi_level'`（`experiments/20261008_mace_smoke/polar_m_cpu.log`）。
  - 第 2 条兜底里原计划的"带电分子整分子退回 GAFF2"作废；改为奇电子数整分子兜底（只做闭壳层，spin=1）。
  - 冒烟：ASE `FixInternals` + LBFGS 10° 一步扫描，乙醇、醋酸根（−1）每步 2.4–4.2 s，键长最大变化 0.54 pm
    （`polar_m_ase_scan.log`）。
- **施工方**：2026-10-08 用户原话"你继续做了。。。你完全接手这个项目了"。MVP-2 起由代理写计划并施工；
  roadmap §3a 的"由对接人员负责，代理不写施工计划"作废（材料说明部分仍可参考）。
