# forcefield_compiler 验收审查（2026-10-04，Asia/Tokyo）

**GROMACS 状态更正（同日补验）：** 用户指出 `/home/user/gmx26.3` 后，已加载 `/home/user/modulefiles/gromacs2026.3` 并运行真实 2026.3 引擎。首轮 PATH 未找到 gmx，不能据此推断未安装。原预留测试加载后失败；独立 12×4 数值矩阵完成，按现有文件比较容差 40 个通过、8 个超差。详见 [GROMACS 补充验收报告](../gromacs_audit_20261004/report.md)。此前两项 P1 仍未修复。

结论：**MVP-0 正常固件路径已通过复测，但本次正式验收不通过，需修复两个 P1 缺陷。** 当前交付范围为 IR、AMBER 导入、OpenMM/AMBER/GROMACS 导出和等价检查。MVP-1 Hessian、MVP-2 torsion provider、MVP-3 charge provider 等仍是计划或参考代码，不属于已实现成果。正常转换测试通过不代表参数准确度或凝聚相科学验证。

本次只增加审查材料、运行日志和诊断产物，未修改实现、既有决策、原始固件或原始 G0 报告；未安装或更改环境。git 未初始化，按决策 0001 不处理版本控制。

## 证据与环境

- 本次主机 `node-A`，使用已有 `/home/user/miniforge3/envs/openmm_dev/bin/python`，Python 3.12.13、numpy 2.4.3、OpenMM 8.5.2、ParmEd 4.3.1、pytest 9.1.1。
- 完整 pytest：**158 passed, 1 skipped，26.10 秒，退出码 0**。唯一 skip 是缺少 gmx。日志：`pytest.log`；结构化结果：`pytest.xml`。
- 复测包括 12 个 GAFF2 固件各 4 个构型的跨引擎测试、故意篡改 bond 的检测、多路径/缺失 torsion 的 1–4 测试、带电分子往返、多项 torsion 及 improper 顺序测试。
- 原始 `jobs/gate_mvp0_20261004/`：12×4 构型，所有已执行腿通过；gmx 和部分乙腈 sander-vs-numpy angle/dihed 行跳过。后者是决策 0004 的豁免，不能算作已经验证的行。
- 原始报告来自 `login-nodegw2.cluster.local`，numpy/scipy 与本次环境不同，故未把原始计时或版本当成本机复跑数据。
- 核对原始 Gate 的 12 个 prmtop SHA256 和固件 manifest 的全部 60 个文件 SHA256，全部一致。
- 原始 Gate 未保存其 pytest 子进程的完整 stdout/stderr，也未记录源码摘要；本次提供独立日志和源码摘要。原始 Gate 的随机数生成序列与参数化 pytest 的逐分子重置序列不同，本次没有声称逐点复跑原 Gate。

## 发现与修复验收条件

### 1. [P1] Gate 忽略 pytest 退出码，测试失败仍可能通过

位置：`scripts/gate_mvp0.py:38-55`、`:129-132`。

`_pytest_summary()` 假设摘要顺序是 passed 后 failed；实际失败摘要可以是 `1 failed, 142 passed ...`，正则只抓到 passed，failed 被置 0。收集报错或空输出也把 failed 置 0。最终判定只检查 failed==0，完全没有使用已经保存的 returncode，也未约束 skip。

已用 mock 隔离 pytest 子进程返回值复现：returncode=1 的测试失败和 returncode=2 的收集报错均满足现有 Gate 的通过条件。任意非 gmx skip 也会被接受。证据：`adversarial_results.json`。

同日加载用户 GROMACS 后的真实复跑也证实：`jobs/gate_mvp0_20261004_2/acceptance.json` 的 external pytest returncode=1，却记录 failed=0、passed=true；外层 Gate 退出 0。详见 GROMACS 补充验收报告。

验收条件：任何非零 pytest returncode、收集错误、无测试/空或不可解析结果都让 Gate 非零退出、acceptance.passed=false；skip 只允许确定的 gmx 缺失用例。使用 pytest JUnit XML 或等价结构化结果统计，并保留两轮完整 stdout/stderr。补回归用例覆盖失败、收集错误、无测试、异常 skip 和正常通过。

### 2. [P1] 导入器改变原始 1–4 Hamiltonian，CLI 等价检查仍通过

位置：`src/ffcompiler/importers/amber_prmtop.py:214-228`；`src/ffcompiler/cli.py:44-49`。

导入器用键图重建 pairs14，但对源文件的检查仅比较 `prmtop_pairs | from_dihedrals` 的并集。源 prmtop 的排除表已经包含 1–4 对，因此这个并集不能证明每对确实有一个有效的、ignore_end=false 的 proper 载体，也不能检查载体重复计数。

复现：复制乙醇 prmtop，仅把所有 proper 条目的 ignore_end 设为 true，禁用源文件的 1–4 作用。`from_amber()` 仍成功，IR 恢复了 12 个 pairs14。OpenMM 读取源文件的总能量为 **25.7234629402 kJ/mol**，IR 为 **−8.4485739838 kJ/mol**，差 **−34.1720369239 kJ/mol**。

CLI 调用 compare 时没有把原始 prmtop/inpcrd 作为参照传入，各后端只互相比转换后的 IR，因此 `cli_style_compare_passed=true`。显式提供源参照后，sander 的 vdw14 和 eel14 两行失败；源值均为 0，导出值分别为 0.1708 和 −8.3378 kcal/mol。证据：`ignored14.prmtop`、`adversarial_results.json`、`source_reference_results.json`、`ignored14_source_check/`。

验收条件：若源文件的有效 1–4 载体、计数或缩放不符合 IR 支持的规则，导入应明确拒绝，不应重建后静默接受。检查源 exclusions 和活跃 1–4 载体的语义，不能只比较并集。CLI 的 check-equivalence 必须包含源 Hamiltonian 的独立参照。上述文件应拒绝导入或让 CLI 非零退出；补部分缺失载体、重复载体、异质缩放及正常环体系回归。

### 3. [P2] IR 校验允许非有限数和无效原子参数

位置：`src/ffcompiler/ir/validation.py:112-139`、`:163-node-A`。

独立修改一个原子、分配独立 atom_type 后，partial_charge=NaN、sigma=NaN、epsilon=−0.5、mass=−1 均被 validate 接受。NaN 的电荷守恒比较不触发错误；atom_type 一致性只比较同类型之间是否一致，不能替代数值范围检查。Bond 的正值检查也会允许正无穷，torsion phase/force_constant 没有有限性检查。

验收条件：对所有数值物理字段检查有限性和合法范围；mass>0、sigma/epsilon>=0（保留 GAFF 零 LJ 参数）、电荷有限；按锁定模型确认 torsion 系数范围。抛带参数位置的 IRValidationError，并用反例回归验证。

### 4. [P2] gmx 可用时仍跳过，预留测试不比较能量

位置：`src/ffcompiler/validation/engines.py:364-368`；`tests/unit/test_export_gromacs.py:171-196`。

compare 在有 gmx 时依旧记录 SKIPPED。预留测试只找日志中的 `Potential Energy` 字符串，没有读取数值或比较 numpy-IR；仅 ethanol 一个分子，且不含 improper。不满足 G0.2 的“如果有 GROMACS 可用，则必须通过”。首轮未在 PATH 找到 gmx；同日加载用户模块后已证实安装可运行，并完成补验，详见顶部更正。

验收条件：gmx 可用时实际求能量并与同一构型的独立参照比较；缺失时才按 D1 跳过。解析失败或超容差必须失败，含 benzene/NMA improper 和带电/环固件。正式比较应使用保留原始精度的坐标；当前 conf.gro 只写 0.001 nm 精度，若直接作为 rerun 坐标，必须处理量化引入的构型差异。

### 5. [P2] 决策 0003 把 GROMACS funct 4 与 funct 2 混淆

位置：`docs/decisions/0003-gromacs-improper-funct4.md:4-9`。

GROMACS 官方拓扑表明确规定：funct 2 是 harmonic improper，funct 4 是 periodic improper，参数为 phase、k、multiplicity。当前导出器写 funct 4 的周期形式是正确的；文档所称 OpenMM/真 gmx 在该函数编号上的分歧没有成立依据。不要据此把导出器改成谐振近似。

依据：[GROMACS 官方拓扑格式表](https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html)。验收条件：以追加勘误或可追溯修订纠正文档，区分已经证实的函数形式与仍未执行的真实 gmx 测试。

## Gate 判据映射

| 判据 | 审查结论 |
|---|---|
| G0.1 | 正常固件已执行腿通过；部分乙腈行有豁免，Gate 聚合器有假阳性缺陷 |
| G0.2 | 首轮 PATH 缺 gmx；补验已找到并运行 2026.3。原测试失败，独立矩阵 8/48 超差，“可用时必须验证”生产分支未完成 |
| G0.3 | 本次独立 pytest 158 passed / 1 gmx skip；原提交缺完整 pytest 日志 |
| G0.4 | 现有 Review Focus 1–5 对应正常/既定反例通过；新发现 1–4 语义和有限性漏洞 |
| G0.5 | 原计划 bond k×1.001 篡改检测测试存在且本次通过 |
| G0.6 | 决策 0001–0004 存在；0003 有事实错误，不能据此认证全部偏离合理。当前文件不能证明 0002–0004 已获用户确认 |

建议先修复 1、2，再处理 3–5 并重新提交 G0 材料。MVP-1 可继续工程工作，但不得把现有 acceptance.passed=true 当作可靠的正式验收签字。

## 复现命令

已有完整套件复跑命令保存在 `pytest.log` 的 START 行。诊断探针：

```bash
PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1 /home/user/miniforge3/envs/openmm_dev/bin/python jobs/acceptance_audit_20261004/reproduce.py
```

`reproduce.py` 打印当前缺陷，输出是诊断证据，不应解释为成功验收。
