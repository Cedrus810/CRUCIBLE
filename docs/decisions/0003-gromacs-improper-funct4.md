# 决策 0003（勘误）：GROMACS funct 4 是周期 improper

- 日期：2026-10-04
- 勘误：最初记录把 funct 4 说成谐振 improper，并据此推断 OpenMM 与真 gmx 语义不一致。这是错误解释，原始文字保留在基线 commit 02acbd9 中。
- 正确语义：funct 2 是谐振 improper；funct 4 是周期 improper，形式 k[1+cos(nφ−δ)]，参数为 phase_deg、k、multiplicity，与本 IR 一致。
- 依据：[GROMACS 官方拓扑格式表](https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html)。
- 实测（2026-10-04，OpenMM 8.5.2 `openmm/app/gromacstopfile.py`）：OpenMM 解析 funct 4 时走 `'1','4','9'` 同一分支（**周期扭转**），参数读作 fields[5]=phase(deg)、fields[6]=k、fields[7]=periodicity；因此写 `i j k l 4 phase k n` 在 OpenMM 读回时与 IR 严格一致（T8 用例 1 验证）。
- 同日补验：用户提供 `/home/user/gmx26.3`，加载既有模块后已实际运行 2026.3。benzene/NMA 的 improper 分项通过当前文件容差；首轮 PATH 缺少 gmx 不能作为未安装的依据。
- 决定：
  1. 保持 funct 4 和原子存储顺序，不做谐振近似。
  2. 真实引擎检查接入生产 compare，报告构建精度并按严格文件容差比较分项和总能量。运行或解析失败不能跳过。
- 影响：修正文档事实，不改变 IR 函数形式。混合精度超差另外记录，不应归因于 funct 4 语义。
