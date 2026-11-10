# 2026-10-28 审阅：10-24 代码（感知层 / fflib / localhess / CP2K / IR 周期盒）

只读审阅的证据目录；没有改源码。kasuga 在 node-A 上运行。

## 1. 单测 `pytest -m "not external"`

```bash
cd /home/ruigengji/CRUCIBLE
OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 PYTHONPATH=src \
  /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest -m "not external" \
  -p no:cacheprovider -q --junitxml=experiments/20261028_review_pytest/junit.xml
```

结果：1 failed, 573 passed, 59 deselected（125 s）。日志 `pytest_not_external.log`。
唯一失败：`tests/unit/test_torsion_mace.py::test_default_model_path_is_polar_m`。这个测试把
`/home/user/MLP/...` 写死在断言里，而代码用的是 `~` 展开。

## 2. 周期盒导出检查 `box_check/`

`box_check/check_box.py`：乙醇加 3 nm 正交盒，比较 NumPy 参考能量和 OpenMM 能量（分连续坐标、
被 wrap 打散的坐标两种），并试读 AMBER 导出文件。输出在 `box_check/check_box.log`：

- NumPy：连续坐标 −8.4486，wrap 后 −8.4486（kJ/mol）；
- OpenMM：连续坐标 −8.4486，wrap 后 8.84e6，各 Force 的 usesPeriodicBoundaryConditions 全是 False；
- sander 读 `box_check/amber/molecule.prmtop`（ntb=1）报 `ERROR: Flag "SOLVENT_POINTERS" not found in PARM file`，
  见 `box_check/amber/sp.out`。

## 3. 修复后（决策 0025）

- `pytest -m "not external"`：579 passed，0 failed（104 s），日志 `pytest_not_external_after_fix.log`；
- 相关 external 测试：`tests/unit/test_ir_box.py tests/unit/test_torsion_mace.py tests/unit/test_cli.py` 共 34 passed
  （设了 `FFC_MACE_MODEL=/home/ruigengji/MLP/mace/MACE-POLAR-1-M.model`，PATH 里加了 openmm_dev 的 bin，并设了 AMBERHOME）；
- 用同一个脚本复查（`check_box_after_fix.log`）：wrap 坐标下 NumPy 和 OpenMM 都是 8837156.5497；ParmEd 能读入；
  sander 在 `ntb=1` 下读 `box_check/after_fix/amber/molecule.prmtop` 正常运行，见 `after_fix/amber/sp.out`。

## 4. G0–G3 重跑（`jobs/gate_mvp{0,1,2,3}_20261028/`）

- 各 Gate 的数值判据与 10-18 第 2 轮一致（G0 report.md 逐字相同，G2 metrics.csv 逐字相同）；
- 四个 Gate 都判为未通过，唯一原因是 external 测试 `test_cli_derive_bonded_ethanol_gfnff` 失败：
  测试的子进程没有 HOME，`~` 被解析成 `/home/kasuga`；
- `gfnff_cli_home_check.log`：同一条 CLI 命令，不传 HOME 时 exit=1，传 `HOME=/home/ruigengji` 时 exit=0。

## 5. 修好测试后的第 2 轮（`jobs/gate_mvp{0,1,2,3}_20261028_2/`）

- `tests/unit/test_assign.py` 的 CLI 子进程现在会继承 HOME 和 FFC_XTB；在 kasuga 账号下分别只设 FFC_XTB、只设 HOME，都能通过（`test_assign_fix_check.log`）；
- G0–G3 全部 passed；四个 Gate 的 report.md 与 10-18 第 2 轮逐字相同；pytest 579 + 60，0 失败。
