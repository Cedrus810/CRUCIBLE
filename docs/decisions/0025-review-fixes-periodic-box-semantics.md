# 决策 0025：10-24 代码审阅修复——周期盒语义、AMBER 周期 prmtop、CP2K 几何、Sobtop 匹配语义

2026-10-28。审阅 10-24 那一批代码（决策 0022/0023）时查出问题，用户指示"该修全修"。
证据目录：`experiments/20261028_review_pytest/`。

## 1. 周期盒只是模拟盒元数据，不进入参考能量（修订 0023「IR 周期盒」第 3 条）

**问题**（实测，`box_check/check_box.log`）：0023 让 `ir/energy.py` 在有盒子时先 unwrap、
再对非键取最小镜像。但三个导出器、以及跨引擎比对实际计算的都是非周期哈密顿量：OpenMM
导出用 NoCutoff，各 Force 的 `usesPeriodicBoundaryConditions` 全为 False；sander 比对用
`ntb=0`；GROMACS 比对把盒子放大到所有镜像都在截断之外。乙醇放进 3 nm 盒、坐标被 wrap
打散时，NumPy 得 −8.45 kJ/mol，OpenMM 得 8.8×10⁶ kJ/mol。0023 里「OpenMM 侧在无 cutoff 下
与最小镜像等价」这句不成立：NoCutoff 根本不使用周期边界。设计稿 §19 要求参考和三引擎
是同一个哈密顿量；而「最小镜像、无截断的普通库仑」没有任何引擎会算，所以这个语义没有
可对照的对象。

**裁定**：`box_nm` 是写进输出文件的模拟盒（gro 盒行、OpenMM 默认盒矢量、prmtop
IFBOX/BOX_DIMENSIONS 与 rst7 盒行），不改变 IR 的哈密顿量。`evaluate` 去掉 `box_nm` 参数、
unwrap 和最小镜像，求的是非周期能量，和各引擎腿是同一个量。坐标必须是完整分子，被 wrap
打散的坐标在每一端都得不到正确能量（现在两端结果相同）。修复后实测：wrap 坐标下 NumPy
与 OpenMM 都是 8837156.5497，连续坐标下两者差 1.3e-7。

周期体系真正的哈密顿量（PME / 截断加修正）取决于模拟设置，不属于本项目的 IR（设计稿
§2 的边界）。如果以后要做周期跨引擎腿，比较对象必须是各引擎实际使用的同一种静电方法，
那是新立项的事。

## 2. AMBER 周期 prmtop 补齐 sander 必读段

**问题**（实测）：只写了 `IFBOX=1` 和 `BOX_DIMENSIONS`，sander 读入时报
`ERROR: Flag "SOLVENT_POINTERS" not found in PARM file`。

**修复**：IFBOX=1 时同时写 `SOLVENT_POINTERS`（IPTRES=1、NSPM=分子数、NSPSOL=NSPM+1，
表示没有溶剂）和 `ATOMS_PER_MOLECULE`（按键图连通分量计数；分子的原子编号不连续时报错，
因为这个格式表达不了）。rst7 盒行改为 6 值 `a b c 90 90 90`。修复后 sander 在 `ntb=1`
下正常运行，ParmEd 也能读入；外部测试 `test_sander_reads_periodic_prmtop` 验证 ntb=0 单点
能量与无盒导出相同。

## 3. GROMACS 三斜盒检查下三角

gro 盒行只能表达 v1(y)=v1(z)=v2(z)=0 的盒子。现在与 OpenMM 导出做同一检查，不满足就报错，
不再写出非法盒行。

## 4. CLI 的 cp2k-out 不再用输入几何覆盖 Hessian 几何

**问题**：`derive-bonded` / `global-fit` 读 CP2K 时传了 `positions_nm=输入几何`，于是
`_check_topology` 拿输入几何和它自己比较，键长比恒为 1；文件里的 CP2K 坐标被丢掉，在别的
几何上算出的 Hessian 也会被静默接受。

**修复**：几何取文件里的 QUICKSTEP 坐标表（与 orca-hess、fchk 同样处理），交给既有的
拓扑检查与输入几何比对。文件没有坐标表时，读取器抛 `MissingCoordinatesError`，CLI 这时
才退回输入几何，并在 stderr 写明「Hessian 几何与输入几何的一致性未核对」（与 `--ir-in`
缺三件套时的提示同一做法）。

## 5. bondcrit 窗口是替换，不是并集（修订 0022 的实现，0022 的文字本来就对）

**问题**：`perceive_bonds` 实现成「半径判据 ∪ 窗口」，测试也按并集写。0022 和逆向报告
§4.1 写的是「在窗口内**才**算成键」。

**证据**：`/tmp/sobtop_re/decompiled.c` 的 `sobtop_IP_setbondcrit_`：先为每个元素对建
下限矩阵（初值 0）和上限矩阵（初值 (covr_i+covr_j)×bondcrit）；bondcrit.dat 的每一行把
该元素对的下限、上限两个矩阵元**整体改写**（÷0.529177 换成 Bohr）。所以列出的元素对
只看窗口。

**修复**：列出的元素对改为只按窗口判定；测试补上「满足半径判据但落在窗口外 → 不成键」
的情形。

## 6. 负 pn 只取绝对值，phase 不变（修订 0022 第 2 条）

0022 写的是「负 pn 按 cos 偶函数等价改写：pn→|pn|、phase→−phase」。AMBER 格式里负
周期性只是「这一组还有下一项」的标记，phase 原样使用；逆向报告 §7.2 记录的 Sobtop
行为也是「直接取 |pn|」，phase→−phase 是当时补出来的推断，不是观测结果。查过出厂
`bonded_param.dat`，没有 pn<0 且 phase∉{0,180} 的条目，所以出厂数据的结果不变。

## 7. improper 库匹配按 Sobtop findffimp（修订 0022 的实现）

原实现只试正序和反序，反序会把中心原子挪出第三位；外围原子的顺序也会影响能否命中
（例如中心 c、外围按编号 (c3, o, o) 时匹配不到 X-o-c-o）。按 `sobtop_IP_findffimp_`
的比较结构改为：第三列等于中心类型；第四列必须等于某个外围原子的类型（这一列不和 X
比较）；前两列匹配其余两个外围原子（可以是 X）；外围顺序任意；按文件顺序取第一条。

## 8. 其他

- `test_default_model_path_is_polar_m` 原先把 `/home/user/...` 写死在断言里（发布时 filter-branch
  改写所致），换任何用户运行都会失败。改为断言默认路径相对 HOME，并用 monkeypatch 检验回退逻辑。
- `hessian/cp2k.py` 删掉两处死代码（无效的 `if …: pass`）。
- `test_assign.py::test_cli_derive_bonded_ethanol_gfnff` 的 CLI 子进程环境原先只有 PYTHONPATH 和 PATH，`~` 会按当前账户的家目录解析（非原作者账号下 xtb 找不到，10-28 第 1 轮 Gate 因此判未通过）。现在把 HOME 和 FFC_XTB 传下去；PATH 仍保持精简。在非原作者账号下分别只设 FFC_XTB、只设 HOME，两种都能通过（`experiments/20261028_review_pytest/test_assign_fix_check.log`）。

## 不做的事

- 默认路径相对 HOME 是用户为发布做的选择（10-25 提交），不改。以非原作者账号运行时，
  用 openmm_dev 的 python 绝对路径，并设置 `FFC_XTB` / `FFC_GMX` / `FFC_MACE_MODEL`。
- localhess 的数学（键、角的 B 行和柔量取倒数）审过，没有问题，不动。
