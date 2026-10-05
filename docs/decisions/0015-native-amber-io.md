# 决策 0015：AMBER prmtop / rst7 原生读写，去掉 ParmEd 依赖

- 2026-10-05 用户要求：openmm / gmx / amber 导出原生支持，不用第三方库；prmtop 导入一并改为原生。
- 现状：GROMACS 导出本来就是原生文本；OpenMM 导出用 OpenMM 自身 API（目标引擎，不算第三方库）；
  ParmEd 用在 AMBER 导出、prmtop 导入、`validation/engines.py`（给 tleap prmtop 配 rst7）、
  `scripts/make_fixtures.py` 的核对和三个测试文件。
- **实现**：
  - `ffcompiler/formats/amber.py`：格式层。按 `%FLAG`/`%FORMAT` 与 Fortran 定宽格式
    （`20a4`、`5E16.8`、`10I8`、`1a80`、`1I8`）解析/写出原始数组；按 POINTERS 核对每个数组长度，
    缺段、截断、重复段、坐标不足都报 `PrmtopError`。ASCII rst7 读写（6F12.7）。
  - 导入器（`importers/amber_prmtop.py`）建在格式层上，语义与 ParmEd 相同：项的顺序先不含氢、
    后含氢；原子下标 = 数组值 // 3；二面角第 3 个原子为负 = 不计 1-4，第 4 个为负 = improper；
    LJ 取 A/B 对角元（Rmin/2 = (2A/B)^{1/6}/2，ε = B²/4A）。**新增**：核对非对角 A/B 符合
    Lorentz–Berthelot（相对 1e-6），不符合就拒绝（IR 只支持 LB，不再默默假设）；元素只按
    ATOMIC_NUMBER 定，不按质量猜。
  - 导出器（`exporters/amber.py`）原生构建全部 FLAG：K = k/2（½k 约定）、含氢/不含氢分列、
    参数类型按取值去重、LB 组合 A/B、前向排除表与 0 哨兵、1-4 载体规则不变；原子 0 落在
    第 3/4 位时整体反序（φ 不变）。
- **验证**：
  - 12 个固件：原生导入的 IR 与改动前 ParmEd 导入的 IR 快照逐项一致（浮点 1e-12，坐标逐位相同）。
  - 12 个固件：每个 FLAG 的原始数组与 ParmEd `parm_data` 完全相同（ParmEd 存在时作为对照测试，
    `importorskip`）。
  - IR → 原生 prmtop → 原生导入，逐项还原（E16.8 精度内），含顺序、方向、1-4 载体与排除表。
  - 篡改测试（非均匀 SCEE、缺失/重复 1-4 载体、非 LB 的 LJ）改为直接编辑原始数组，全部被拒。
  - 子进程检查：导入 importer 与 AMBER exporter 后 `parmed` 不在 `sys.modules`。
  - 全套 pytest（非 external 355、external 24）通过。
  - 正式 Gate：`jobs/gate_mvp0_20261005_3/`（G0 passed，2596 行 0 失败）、
    `jobs/gate_mvp1_20261005_4/`（G1 passed）。之前一次运行（`gate_mvp0_20261005_2`、
    `gate_mvp1_20261005_3`）引擎各腿全部通过，但 `test_importer_does_not_load_parmed` 的子进程
    没有继承 pytest.ini 的 pythonpath 而失败；已改为显式传 PYTHONPATH，失败运行保留。
- ParmEd 从 `runtime._PACKAGES` 与 `pyproject.toml` 可选依赖中移除；只在对照测试里可选使用。
- 不改变：单位、½k 约定、1-4 语义、IR 结构；provenance 的 `parmed_version` 改为 `prmtop_reader`。
