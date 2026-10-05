# GROMACS 2026.3 补充验收（2026-10-04，Asia/Tokyo）

结论：用户提供的 GROMACS **已安装且可在 node-A 运行**。首轮审查根据 PATH 未找到 gmx 推断“本机未安装”，该判断不正确。加载用户的现有模块后，原预留测试失败；独立数值审查完成 12 个分子 × 4 个构型，按现有文件比较容差 **40 个通过、8 个超差**。正式验收仍不通过。

本次只增加审查脚本、运行产物和报告，未修改 compiler 的源码、测试、modulefiles、系统配置、原始固件或环境中的包。此前发现的两个 P1 缺陷仍未修复。

## 安装与运行环境

- 可执行文件：`/home/user/gmx26.3/bin/gmx`，实测版本 **2026.3**。
- 当前 CPU：Intel Xeon Gold 6138；按用户 modulefile 选择普通版 `AVX2-fallback`，不使用 `gmx26.3_512`。
- 构建：**mixed precision、thread_mpi、AVX2_256**；审查时指定 CPU nonbonded/bonded、1 MPI rank、1 OpenMP thread。使用 CPU 是单点正确性检查的实现选择，不代表 GPU 检查或性能验证。
- modulefiles 目录：`/home/user/modulefiles`；加载 `gromacs2026.3`，自动加载 cuda13.3、openmpi5.0.10、HDF5。
- 直接运行而未加载模块时，动态加载器报 `libhdf5.so.310` 缺失；现有 HDF5 模块提供所需路径，随后 gmx --version 成功。没有创建或修改任何兼容链接。
- 版本原始输出：`version.txt`；完整矩阵运行 manifest：`../gromacs_engine_audit_20261004_3/manifest.json`。

可在独立 shell 中复现，不需修改用户启动文件：

```bash
source /etc/profile.d/env-modules.sh
module use /home/user/modulefiles
module load gromacs2026.3
PATH=/home/user/miniforge3/envs/openmm_dev/bin:$PATH \
PYTHONDONTWRITEBYTECODE=1 \
/home/user/miniforge3/envs/openmm_dev/bin/python \
jobs/gromacs_audit_20261004/run_real_gromacs.py
```

## 原预留测试：真实失败

`tests/unit/test_export_gromacs.py::test_real_gromacs_rerun` 加载 gmx 后不再 skip，结果 **1 failed**，在 grompp 阶段退出码 1。独立重复同一 grompp 命令保存了完整 stderr：

1. `With Verlet lists only full pbc or pbc=xy with walls is supported`：测试写的 `cutoff-scheme=Verlet, pbc=no` 不受支持。
2. 非周期维度没有盒子体积，却启用了默认正值的 verlet-buffer-tolerance，自动计算 rlist 失败。

这不是“不存在 gmx”的问题。除此之外，该测试的 mdrun 命令未指定生成的 `run.tpr`；即使前处理修好，还需传 `-s run.tpr`，不能默认读不存在的 topol.tpr。最终只断言日志字符串存在，也仍不能证明能量一致。

证据：`existing_test.log`、`existing_test.xml`、`existing_grompp.stdout`、`existing_grompp.stderr`、`existing_test_processed.mdp`。

## 原 Gate 的真实假阳性

加载 gmx 后还原样运行了 `scripts/gate_mvp0.py`，输出到 `../gate_mvp0_20261004_2/`，没有覆盖原始 Gate。其 `acceptance.json` 同时包含：

```json
{
  "passed": true,
  "pytest_summary": {
    "external": {"passed": 15, "failed": 0, "skipped": 0, "returncode": 1}
  },
  "skipped_legs": {
    "gmx": "gmx leg not implemented in MVP-0; run test_real_gromacs_rerun"
  }
}
```

此处为字段节选。非 external 的 143 项通过，external 的 15 项通过且进程失败；结合单项原测试的失败日志，确认真实环境中的测试失败被摘要解析漏报，并且 Gate 忽略了 returncode=1。Gate 的外层进程仍退出 0。这将此前 mock 复现的 P1 提升为实际复跑证据，不能把该文件的绿色状态解释为通过验收。

证据：`gate_with_gromacs.log`、`../gate_mvp0_20261004_2/acceptance.json`、`../gate_mvp0_20261004_2/manifest.json`。

## 独立数值矩阵

脚本：`run_real_gromacs.py`。完整产物：`../gromacs_engine_audit_20261004_3/`；摘要：`matrix_summary.json`；外层进程日志：`matrix.log`。

方法：使用当前导出器生成的 top/itp，通过真实 grompp + mdrun -rerun 求能量，再从 edr 经 gmx energy -dp 提取各分项及 Potential，逐行对照 numpy-IR。12 个固件、每个原始构型加 3 个 N(0, 0.005 nm) 扰动，使用原 Gate 的分子顺序和共享 RNG 种子 20261004。每个分子的 4 帧批量 rerun，检查帧数和时间索引，未忽略失败分子。

为满足现代 GROMACS 的 Verlet 要求，使用 10 nm 周期盒、4 nm cutoff，关闭 Coulomb/LJ potential modifiers、PME 和 dispersion correction。实测所有分子内距离均小于 cutoff，坐标范围紧凑且远小于盒子，使截断包含所有分子内非键对，周期镜像不进入 cutoff。该单点设置与此次 gas-phase NoCutoff 参照比较具有相同相互作用集合，不是凝聚相模拟。

保留 compiler 原始 0.001 nm 精度的 `conf.gro`，另写审查用 `audit.gro`，坐标精度 1e-10 nm。对每个分子全部构型施加同一个刚性平移，使坐标在正值小范围内；检查平移前后 numpy 各分项变化不超过 1e-9 kJ/mol，并同时保存原始及平移后的坐标数组。该平移不改变分子构型，用于避免混合精度周期包裹在小分子负坐标上引入额外误差。

判据采用现有文件比较公式 `abs(delta E) <= 1e-4 + 1e-7*abs(E_ref)` kJ/mol，应用于 gmx 数值审查；计划 §2.4 未单列真实 gmx 的精度政策。没有新增 mixed-precision 豁免、放宽容差或跳过乙腈近线性角。总能量及所有分项均为 40/48 通过；所有 48 个构型完成了计算和解析，没有执行错误。

| 超差构型（config 从 0 开始） | 总能量绝对差 kJ/mol | 总能量容差 kJ/mol | 超差分项 |
|---|---:|---:|---|
| ethanol config 1 | 0.000163185 | 0.000111161 | bond、total |
| benzene config 3 | 0.000190849 | 0.000114671 | bond、total |
| cyclobutane config 3 | 0.000122304 | 0.000121699 | total |
| norbornane config 2 | 0.000134314 | 0.000119393 | bond、total |
| norbornane config 3 | 0.000153083 | 0.000118879 | bond、total |
| acetate config 1 | 0.000128670 | 0.000123860 | total |
| methyl_benzoate config 3 | 0.000218393 | 0.000127296 | bond、total |
| acetonitrile config 0 | 0.001123166 | 0.000102878 | angle、total |

最大分项差：bond 0.000168972；angle 0.001122927；proper 0.000018934；improper 0.000011732 kJ/mol。各非键分项均在当前判据内通过。含 improper 的 benzene、NMA 并未出现决策 0003 声称的函数字义分歧。

这些小差值与混合精度/近线性角的数值问题相容，但当前只证明存在超差，**未证明全部超差的根因**，不能据此将其自动归类为可接受误差。后续应在双精度 GROMACS 上对同一输入补验，或在有定量证据并获得接受的情况下明确 mixed-precision 验收政策；本次没有编译或安装新 GROMACS。

## 保留的诊断历史

- `../gromacs_engine_audit_20261004/`：首个乙醇诊断未施加平移，保留负坐标，原始高精度输入仍出现 total 差 0.000387274 kJ/mol。
- `../gromacs_engine_audit_20261004_2/`：施加能量不变的刚性平移后，ethanol、benzene、NMA 的原始构型三者均通过。
- `../gromacs_engine_audit_20261004_3/`：上述完整矩阵，40/48 通过。
- 对这几次脚本，外层 shell 退出码以 `pilot.log`、`translated_pilot.log`、`matrix.log` 为准。当前 RunContext 只依据 with 是否抛异常写 exit_status；审查脚本最初在正常退出 with 时返回 1，因此失败诊断目录的 manifest/run.log 显示 exit=0，而外层日志正确显示 exit=1，acceptance.passed=false。该历史不被覆盖。

## 对验收结论的影响与下一步

将首轮报告的“本机缺 gmx”改为“首轮 PATH 未发现；加载用户现有模块后已可运行”。G0.2 当前不能通过：生产 compare 在找到 gmx 时仍无条件记录 SKIPPED，原预留测试失败，独立数值矩阵也未全部达到采用的文件容差。

先修此前两个 P1（Gate 对 pytest 失败的假阳性、源 prmtop 1–4 语义改变），再将真实 GROMACS 适配器、可运行的 MDP、实际能量解析比较接入 Gate，补记录精度政策。保留本次 8 个超差用例作为回归材料；不要把审查脚本的独立检查当作生产 compare 已实现。
