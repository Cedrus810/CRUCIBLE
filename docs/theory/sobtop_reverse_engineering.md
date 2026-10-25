# Sobtop 2026.1.16 逆向工程报告

**目标二进制**：`sobtop`（Linux x86-64 ELF，17.8 MB）与 `sobtop.exe`（Windows PE64，14.4 MB）
**作者/发行**：Tian Lu (Sobereva)，http://sobereva.com/soft/sobtop，2026-Jan-16 发行
**逆向方法**：静态（Ghidra 11 headless 反编译 86 个自有符号 + radare2/objdump 交叉验证）+ 动态（stdin 驱动菜单的黑盒实验，乙醇/乙烯/苯/SiO₂ 体系，输出逐字段比对）
**产物位置**：`/tmp/sobtop_re/decompiled.c`（85 函数）、`/tmp/sobtop_re/outitprtp.c`（核心写出器）、`/tmp/sobtop_re/sandbox/`（动态实验沙箱）
**用途声明**：为实现互操作与算法核对的行为级还原；参数数据本身来自公开的 AmberTools (gaff.dat 1.81 / parm19)、UFF 论文 (Rappe 1992) 与文献共价半径 (Cordero 2008)。

---

## 1. 二进制技术底细

| 项 | 内容 |
|---|---|
| 语言/编译器 | Intel Fortran（ifort，`for_write_seq_fmt` 等 runtime 符号） |
| 数学库 | MKL（DGEEV、DGESVD、LAPACK/BLAS）、OpenMP 并行（`OMP_STACKSIZE` 必须 ≥64M，否则 OMP Error #34） |
| 符号 | **未 strip**：16,799 个符号；自有代码 = `MAIN__`（全部交互逻辑）+ 85 个 `sobtop_IP_*` 前缀函数 |
| 最大函数 | `sobtop_IP_outitprtp_`（51 KB，.itp/.rtp 写出器） |
| 数据文件 | `sobtop.ini`（运行参数）、`LJ_param.dat`（LJ 库）、`bonded_param.dat`（成键库）、`bondcrit.dat`（成键判据覆盖）、`assign_AT.dat`（自定义类型规则）、`ATOMTYPE_GFF.DEF`/`ATOMTYPE_AMBER.DEF`+`atomtype` 二进制（AmberTools 的 antechamber atomtype，shell out 调用） |

菜单树（从字符串还原，与动态实验一致）：

```
主菜单: -2 k 缩放因子 | -1 选 Hessian 方法 | 0 退出 | 1 生成 GROMACS 拓扑
        2 .gro | 3 .rtp(pdb2gmx) | 4 .gjf | 5 单项推导 | 6 全部刚性参数
        7 电荷(EEM/Gasteiger←Multiwfn; MMFF94←OpenBabel; 手动; .chg 导入导出)
        10 其他(B 矩阵/RIC Hessian 导出等)
类型子菜单: UFF 全部分配 / GAFF→UFF / 仅 GAFF / 仅 AMBER / assign_AT 规则 / 手动 / 重置
装配子菜单(7 模式): 1 无参数 | 2 全部方法推导 | 3 预置+方法补缺 | 4 =3+乱猜
        | 5 刚性区方法+其余预置+方法补缺 | 6 =5+乱猜 | 7 键角全方法+扭转预置+方法补缺
```

---

## 2. 全局数据结构（`defvar` 模块，从符号与交叉引用还原）

| 变量 | 含义 |
|---|---|
| `a_` | 原子记录数组，每原子 0x50=80 字节：xyz 在 +8/+0x10/+0x18（**Bohr**），类型字符串在 +0x38（10 字符），名称在 +0x42 |
| `hess_` / `hessric_` | Cartesian Hessian (3N×3N，E_h/Bohr²) / RIC Hessian (m×m 全矩阵) |
| `bmat_` | Wilson B 矩阵（m×3N，Bohr 构造） |
| `bond_/angle_/dih_` | 拓扑项数组 {i,j,k,l,...} |
| `ncenter_/nbond_/nangle_/ndih_/nric_` | 计数；`nric = nbond+nangle+ndih`（**不含 improper**） |
| `covr_` | 共价半径表（Å，见 §4） |
| `connmat_` | 连接矩阵（int16 稀疏存储） |
| `cellv1/2/3_`, `ifpbc_` | 晶胞向量与周期开关 |
| `k_method_/bkfac_/akfac_/dkfac_/ikfac_` | 方法选择与四类 k 缩放因子 |
| `rigidatm_/nrigidatm_` | 刚性区原子集 |
| `nffbond_/nffang_/nffdih_/nffimp_` | bonded_param.dat 已加载条目数 |

---

## 3. 输入层

- **mol2**：元素从 **TYPE 列**解析（`elename2idx`）。因此 antechamber 输出的 GAFF 类型（c3 等）**直接拒载**（"Unable to recognize element"）；正确输入是 TYPE 列写元素符号（CRUCIBLE 的 `write_ir_mol2` 正是这么做的）。BOND 块的键级只透传给外部 atomtype 程序，sobtop 自身的扭转处理**不使用键级**。
- **pdb/pqr**：读 CRYST1 → `abc2cellv`；半径判据成键（§4）；>1000 原子不屏显坐标（`ioutatminfo` 可改写文件）。
- **.chg**（Multiwfn）：坐标+电荷；`ichggeom=1` 时**用 .chg 的几何替换当前几何**。
- **Hessian**：Gaussian `.fch/.out`、ORCA `.hess`、CP2K `.out`、xtb `hessian`。`loadhess` **零单位换算** → 输入必须已是 E_h/Bohr²。
- 坐标进入内部即 ÷0.529177210903 转 Bohr；输出时再转回。

---

## 4. 感知层（成键 / 角 / 二面角 / 周期）

### 4.1 成键判据（`setbondcrit` + `genconnmat` + `genneighlist`）

```
bonded(i,j) ⟺ d(i,j)[最近镜像] < (covr(Z_i) + covr(Z_j)) × bondcrit
默认 bondcrit = 1.15（sobtop.ini 可改）
bondcrit.dat 可按元素对覆盖：行格式 "元素1 元素2 下限 上限"（Å 距离窗口，
在窗口内才算成键——注意是绝对距离区间，不是因子），加载时展开进 150×150 判据矩阵
```

**共价半径表**（从 `.rodata` @0x160c2e0 直接抠出，= Cordero et al. 2008，Å）：

```
H .31 He .28 Li 1.28 Be .96 B .84 C .76 N .71 O .66 F .57 Ne .58
Na 1.66 Mg 1.41 Al 1.21 Si 1.11 P 1.07 S 1.05 Cl 1.02 Ar 1.06
K 2.03 Ca 1.76 Sc 1.70 Ti 1.60 V 1.53 Cr 1.39 Mn 1.39 Fe 1.32 Co 1.26
Ni 1.24 Cu 1.32 Zn 1.22 Ga 1.22 Ge 1.20 As 1.19 Se 1.20 Br 1.20 Kr 1.16
Rb 2.20 Sr 1.95 Y 1.90 Zr 1.75 Nb 1.64 Mo 1.54 Tc 1.47 Ru 1.46 Rh 1.42
Pd 1.39 Ag 1.45 Cd 1.44 In 1.42 Sn 1.39 Sb 1.39 Te 1.38 I 1.39 Xe 1.40
Cs 2.44 Ba 2.15 La 2.07 Ce 2.04 Pr 2.03 Nd 2.01 Pm 1.99 Sm 1.98 Eu 1.98
Gd 1.96 Tb 1.94 Dy 1.92 Ho 1.92 Er 1.89 Tm 1.90 Yb 1.87 Lu 1.87 Hf 1.75
Ta 1.70 W 1.62 Re 1.51 Os 1.44 Ir 1.41 Pt 1.36 Au 1.36 Hg 1.32 Tl 1.45
Pb 1.46 Bi 1.48 Po 1.40 At 1.50 Rn 1.50 Fr 2.60 Ra 2.21 Ac 2.15 Th 2.06
Pa 2.00 U 1.96（表共 110 项，到 109 号）
```

### 4.2 角与二面角生成（`genall_ang_dih`）

- **角**：每个中心原子的邻居两两组合。含氢、不去重、**不剔除线性角**（乙腈 C-C≡N 的线性角保留，只是 Hessian 方法跳过）。
- **二面角**：所有 i-j-k-l 路径（甲基-甲基 9 条全出）。剔除条件（精确代码）：

```
跳过 ⟺ (180° − θ_ijk < 1.0°) ∨ (180° − θ_jkl < 1.0°)
```

  存储规范化为 i<j 方向。
- 邻居表 + 连接矩阵一次生成；角数/二面角数 >10000 时跳过输出排序。
- `iskipgendih=1`（ini）可整体跳过二面角生成（大体系提速）。

### 4.3 周期处理

- `getbondvec(v, A, B)`：v = B 相对 A 的**最近镜像**位置差，全程 PBC 安全（键向量、距离、角、二面全部走它）。
- 三斜晶胞支持（CRYST1 α/β/γ → 向量）；`.gro` 写 GROMACS 盒向量（FeCl2 实测六方盒正确）。
- 实证：SiO₂ 3×3×3（324 原子）→ 432 键 = 108×4，全部跨边界键正确。

---

## 5. Hessian 参数推导层

### 5.1 方法分派（`k_method`）

| # | 名称 | 键 k | 角 k |
|---|---|---|---|
| 1 | Seminario | `Σλ_t·|u·v_t|` | `1/(1/k_A+1/k_C)`（无 f 因子） |
| 2 | mSeminario | **同 1**（键无区别） | 同 1 **但含 f 因子** |
| 3 | m2Seminario | 二次型 `uᵀ(−H_AB)u` | 二次型 + **同样含 f 因子** |
| 4 | DRIH | `hessric(idx,idx)` | 同左 |
| 5 | **隐藏方法**（菜单无入口） | `genlochess` 局部 6×6 块 + B 行 + 双 SVD 伪逆投影 | 同左（3/4 原子版） |

核心投影（`prjhess_frc`）：

```
H3 = hess(3i−2:3i, 3j−2:3j)         # 3×3 Cartesian 块
method 1/2: A = −H3 → DGEEV("N","V") 非对称本征分解（虚部丢弃）
            k(u) = Σ_t λ_t |u·v_t|    # v_t 逐个再归一化
method 3:   k(u) = uᵀ(−H3)u
method 4/5: 返回值不写入（未定义）→ 见 §8 bug 清单
```

**角公式**（`calcangleparm`，u_PA/u_PC 来自 `getpapc`）：

```
kA = |r_AB|² · k_proj(u_PA)；kC = |r_CB|² · k_proj(u_PC)
method 1: k = 1/(1/kA + 1/kC)
method 2/3: f_A = (Σ_{以 B 为中心、含 A 端的邻角} cos²(u_PA, u_PA′) − 1)/(n−1)   # 扣除本角自身
            （n=1 时 f=0）；f_C 同理
            k = 1/((1+f_A)/kA + (1+f_C)/kC)
```

> **对既有理解的两处修正**：(a) f 因子对 method 2 和 3 **同时生效**，不是 m2Seminario 特有；
> (b) f 的精确形式是"其余共端邻角 cos²(u_PA) 的算术平均"。

`getpapc(a,b,c)`：θ≤0.1° 或 ≥179.9° 时给 B 加位移 **(0.352, −0.192, 0.737) Bohr**（模≈0.839）防退化；u_PA = normalize(t×v₁)、u_PC = normalize(v₂×t)，t = 角平面法向。`crossprod` 平行兜底阈值 |cos|>0.99999（平行时替换 b 为 (1,−1,1)）。

### 5.2 DRIH（`formb`+`formhessric`+`diaggemat`）

```
RIC 集 = 全部 bond ∪ 全部 angle ∪ 全部 proper dihedral（不含 improper）
B 用 Bohr 构造（实证：Å/nm 构造时环系错 1–8 倍）
pinv(Bᵀ) = DGESVD 两步（lwork 查询+实算），奇异值截断 |σ| > 1e-10【绝对阈值】
H_q = A·H_x·Aᵀ，全 m×m 矩阵 —— 非对角块 = 键-键/键-角交叉项（→ §7 funct 3）
参数取对角元；键 hessric(b,b) 单位 E_h/Bohr²，角/二面 E_h/rad²
```

> **阈值差异**：sobtop 是 B 奇异值绝对阈值 1e-10；等价于 G=BBᵀ 本征值 λ<1e-20 置零。
> 数值上等价于相对阈值 σ<1e-4·σmax 只在谱双峰时成立——病态几何（近线性、近冗余 RIC）
> 下不同阈值策略会分叉，是复现的边界条件。

### 5.3 二面角与 improper

- **二面角**（`calcdihparm`）：method 1/2/3 走两端投影调和平均
  `k = 1/(1/k_A + 1/k_C)`，`k_A = |e₁×(−e₂)|²·r₁₂²·k_proj(m₁)`；method 4 取 `hessric(d,d)`。
  平衡值 = **当前几何的带符号二面角**（`xyz2dih_sign`，±(0,180]°）。
- **improper**（`calcimpparm`）：`u = normalize((r_j−r_i)×(r_l−r_k))`，
  `k = k_proj(i,j,u) + k_proj(i,k,u) + k_proj(i,l,u)`（三个 1-2 块求和）。
  **距离型参数化**：d0 = i 到 jkl 平面垂距（`pointprjple` 投影 + `potpledis` 求距），
  k 单位 E_h/Bohr²（×9375.8307 → kJ/mol/Å²）。**不分派 k_method** → method 4/5 下未定义（§8）。

### 5.4 单位换算（代码字面量，精确复现用）

| 换算 | 常量 |
|---|---|
| Bohr↔Å | 0.529177210903 |
| E_h/Bohr² → kJ/mol/Å² | ×9375.8307 |
| E_h/Bohr² → kcal/mol/Å² | ×2240.8771 |
| E_h/rad² → kJ/mol/rad² | ×2625.5 |
| E_h → kcal/mol | ×627.5095 |
| kJ/mol/Å² → kJ/mol/nm² | ×100 |
| bonded_param.dat 的 AMBER 约定 | K = k/2（E=KΔq²，GROMACS ½k 的两倍） |

`prjhess_frc` 输出 ×`bkfac/akfac/dkfac/ikfac`（键/角/二面/improper 缩放，默认 1）。

---

## 6. 原子类型与 LJ 层

- **UFF typing = `"UF_" + 元素符号`**（`ind2name` 表拼接）。纯元素级——真 UFF 的 C_3/C_R/N_2 杂化子类型**完全没有**；等于只借 UFF 的 LJ 表（103 类型）。
- **GAFF/AMBER typing = shell out** AmberTools `atomtype`：`atomtype -i in.mol2 -f mol2 -o MOL.ac -d ATOMTYPE_GFF.DEF`（AMBER 用 ATOMTYPE_AMBER.DEF），读回 MOL.ac。"GAFF 优先然后 UFF" 模式逐原子兜底（FeCl2 lazy 实测：cl 用 GAFF、Fe 用 UF_Fe 混排）。
- **LJ 转换**（`lj_convert`）：库存 (R_half[Å], ε[kcal/mol]) → GROMACS
  `σ[nm] = 2·R_half·2^(−1/6)·0.1`，`ε[kJ/mol] = ε×4.184`。
  实证：UF_C (1.9255, 0.105) → σ=3.430851E-01、ε=4.393200E-01，分毫不差。
- 类型在 LJ_param.dat 缺失 → **报错并要求手工编辑 LJ_param.dat**（不猜测）。
- `[ atomtypes ]` charge 列恒 0；质量表 CODATA 系（C=12.010736, O=15.999405, H=1.007941）。

---

## 7. 参数库匹配与拓扑装配层

### 7.1 bonded_param.dat / LJ_param.dat 解析（`loadbonded_param`）

字段：`$BOND / $ANGLE / $DIHEDRAL / $HARM_DIH / $IMPROPER`，`$end` 结束，`;` 注释。
条目格式（AMBER 风格）：

```
$BOND:     type1-type2                k[kcal/mol/Å²]  r0[Å]
$ANGLE:    type1-type2-type3          k[kcal/mol/rad²] θ0[°]
$DIHEDRAL: type1-type2-type3-type4    IDIVF  barrier[kcal]  phase[°]  pn（可负）
$HARM_DIH: type1-...-type4            K[kcal，AMBER K=k/2 约定]  d0[°]   （出厂为空）
$IMPROPER: type1-type2-type3-type4    barrier[kcal]  phase[°]  pn
```

### 7.2 匹配优先级（`findff*`，实证 hc-c3-c3-oh 压过 X-c3-c3-X）

1. 先试**全具体类型**（正向 type_i-type_j-type_k-type_l；正反两个方向）；
2. 失败再试 **X 通配**条目（X 在任意位）；
3. $DIHEDRAL 的 `barrier/IDIVF` 后换算 kJ（×4.184）；同一类型多条 → **展开为多条 funct 9 行**（不同 pn）；
4. **负 pn 处理**：直接取 |pn|；phase=0 时数学严格等价（cos 偶函数）；非零 phase 的等价改写为 phase→−phase。

### 7.3 improper 检测（`outitprtp` 内联，数据驱动，无硬编码）

```
对每个原子 i：nneigh(i)==3 且 type(i) ∈ {已加载 $IMPROPER 条目的第三列(中心)类型集合}
→ 生成 improper (j₁, j₂, i, j₃)，j 为邻居升序（GROMACS funct 4，第三原子为中心）
```

中心类型集合（当前库）= GAFF {ca, c, c2, n, na, n2} ∪ AMBER {CA, C, N, N*, CM, CJ, CS, CB, CC, C4, NA, N2, CW, CV, CR, CQ, CP, CN, CK, C5}。
UFF 类型不触发**只因库里没有 UF 条目**——往库里加一条带 UF 中心的 improper 即触发。
参数缺失时兜底 "guess (same as GAFF X-X-ca-ha)"：1.1 kcal → kd=4.6024 kJ、phase 180、pn 2。

### 7.4 缺参"乱猜"常量（`add*miss`，实证）

| 项 | 猜测值 |
|---|---|
| 键 | r0 = 当前距离，k = **2.5×10⁵ kJ/mol/nm²**（写死） |
| 角 | a0 = 当前角，k = **500 kJ/mol/rad²**（写死） |
| 二面 | funct 9，phase 0，kd **0**（自由旋转！），pn 1 |
| improper | X-X-ca-ha 等价（1.1 kcal, 180°, pn 2） |

### 7.5 装配决策（funct 2 vs 9、pairs、交叉项）

- **funct 2（谐振）** ⟸ 参数来自 Hessian 方法（d0=当前几何，k=推导值；"全方法"模式把**可旋转键也锁死**）。
- **funct 9（周期）** ⟸ 预置库或缺失（缺失→零势垒）。
- **pairs**：每条 funct 9 二面角的 (i,l) 端点；**剔除**互为 1-2（connmat≠0）或 1-3（同角）的对（环闭合正确豁免）；无序去重；全谐振模式下 pairs 为空。
- **[ angles ] funct 3 键-键交叉项**：仅当 `k_method==4 ∧ 装配模式==2(全方法) ∧ 存在 DRIH 推导的角` 时写出，数据 = `hessric` 的键-键非对角块（实证乙醇 DRIH 全方法模式输出 13 条；mSeminario 模式无）。
- `[ defaults ] 1 2 yes 0.5 0.8333` 写死；nrexcl=3；cgnr=原子号；单残基 MOL。
- FFderiv.txt：按类型类统计 k/r0 的 avg/min/max（插回 bonded_param.dat 用平均值；itp 内写逐项值不平均）。
- `5 单项推导` 的插入提示（"Information inserted into $HARM_DIH"）只**打印**行（`[类型名] K[AMBER] d0`），不改文件。

### 7.6 输出面

- .top/.itp/.gro（GROMACS 专属，无 AMBER prmtop / OpenMM 输出）。
- .rtp = **GROMACS pdb2gmx 残基模板**（[bondedtypes] `1 1 9 4 1 3 1 0`），要求手工搬运 atomtypes、质量到 atomtypes.atp、手工填电荷；且警告 pdb2gmx 后 functype 重复需手动删。
- 排序：键按 (i,j) 升序；角度/二面 >10000 跳过排序。

---

## 8. 隐藏功能与 bug 清单

1. **[BUG] DRIH/方法5 模式下 improper 未定义**：`calcimpparm` 不分派 k_method，`prjhess_frc` 在 method 4/5 不写返回值 → k = 未初始化栈值 × ikfac。交互式装配路径无条件调用它。
2. **[隐藏] k_method=5**：`genlochess`（抽 3n×3n 局部 Cartesian Hessian 块）+ B 行 + 双伪逆投影的局部方法；键/角/二面三处实现完整，**菜单无入口**（作者半成品）。
3. **[物理缺陷] "全方法"模式锁死构象**：可旋转键也变 funct 2 谐振项钉在当前几何；全谐振模式下 1-4 pairs 全部消失。
4. **[物理缺陷] 缺失扭转势垒置零**（自由旋转）；缺键/角用写死常数 2.5e5/500 乱猜。
5. **[互操作缺陷] mol2 的 GAFF 类型拒载**（元素从 TYPE 列解析）。
6. **[退化] UFF typing 无杂化**，与 UFF 论文的类型学不符。
7. **[移植性] funct 3 交叉项只有 GROMACS 能表达**；无跨引擎验证机制。
8. **[边界] 伪逆绝对阈值 1e-10** 在病态几何下与相对阈值策略分叉。

---

## 9. 附录 A：动态实验证据索引（沙箱 /tmp/sobtop_re/sandbox/）

| 实验 | 文件 | 验证点 |
|---|---|---|
| 乙醇裸拓扑 | et.mol2 → 0(MOL.itp) | 8键/13角/12二面计数；[defaults]；LJ 转换；pairs=12 |
| 乙烯+GAFF | ethylene.mol2 → ethylene.itp | improper 触发(2条)、X-c2-c2-X 匹配、guess 兜底 |
| 苯+GAFF | benzene.mol2 → benzene.itp | 环内 improper 6 条、X-ca-ca-X 14.5/4×4.184=15.167 |
| 乙醇+DRIH 全方法 | etopt.mol2 → etopt.itp | 全 funct 2、funct 3 交叉项 13 条、无 pairs |
| 乙醇+DRIH 混合模式7 | 同上（GAFF 类型） | 预置多项 funct 9 展开（ho-oh-c3-c3 两条 pn 3/1） |
| SiO₂ 周期 | SiO2_3x3x3.pdb | CRYST1、432 镜像键 |
| 单项推导 | 菜单 5，4,1,2,7 | 0.015375 a.u.×2625.5=40.368 kJ/mol/rad²；$HARM_DIH 只打印 |

## 附录 B：与 CRUCIBLE 的算法对照（摘要）

| 算法 | sobtop | CRUCIBLE | 判定 |
|---|---|---|---|
| Hessian 四方法 | ✓（含隐藏方法5） | ✓ + 全局拟合 | 平 |
| f 因子 | method 2/3 均生效 | 已实现（表述需修正） | 平 |
| 伪逆阈值 | |σ|>1e-10 绝对 | 1e-8·λmax 相对 | 等价域内平 |
| 对称化 | 无 | 图自同构 | CRUCIBLE 胜 |
| 柔性扭转 | 预置/置零 | MACE-POLAR 拟合 | CRUCIBLE 胜 |
| improper 检测 | 数据驱动类型集 | prmtop 外包 | 可吸收 sobtop 规则 |
| 电荷 | 外部 EEM/Gasteiger/MMFF94 | AM1-BCC/NAGL | CRUCIBLE 胜 |
| 周期/金属/UFF-LJ | ✓ | ✗ | sobtop 独占（可吸收） |
| 导出 | 仅 GROMACS | 三引擎+验证 | CRUCIBLE 胜 |
| DRIH improper | 未定义行为 | — | sobtop bug |

*报告生成：2026-10-24，逆向与分析全程可复现（沙箱脚本与反编译产物在 /tmp/sobtop_re/）。*
