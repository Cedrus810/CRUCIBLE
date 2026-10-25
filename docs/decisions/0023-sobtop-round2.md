# 决策 0023：Sobtop 吸收第二批——localhess 方法、CP2K provider、IR 周期盒

## 范围

决策 0022 吸收了感知层与数据文件解析。本批按用户指示完成三件（并明确
一件不做）：

1. **`bonded/localhess.py`（experimental）**：复刻 Sobtop 隐藏的 k_method=5；
2. **`hessian/cp2k.py`**：CP2K `.out` VIB Hessian 读取器（`--hessian-provider cp2k-out`）；
3. **IR 周期盒**：`ForceFieldIR.box_nm` 字段 + 三导出器 + NumPy 能量镜像；
4. **mol2 结构入口：明确不做**（用户决定）。管线里 `charge/loaded` 已消费
   mol2 的电荷语义，再开一个 mol2 结构入口会让"两种 mol2"在管线中共存；
   结构输入仍走 prmtop/inpcrd/sdf 或管线 sidecar。

## localhess（Sobtop k_method=5 复刻）的歧义与裁决

逆向（报告 §5.1、§8.2）：`calcbondparm/calcangleparm/calcdihparm` 的
`k_method==5` 分支调用 `genlochess`（抽项所涉原子的 3n×3n Cartesian 局部
块），取该项 Wilson B 行 b，SVD 伪逆（Sobtop 绝对阈值 |σ|>1e-10），按
`k = bᵀ·M·b` 二次型收尾（汇编 0x43601e–0x436086 双循环确证）。**该方法在
菜单无入口、sobtop.ini 写 5 也被拒绝（实测回退 mSeminario）→ 无 oracle。**

两处无法从汇编裁决的歧义：(a) M 是 pinv(H_loc) 还是 pinv(pinv(H_loc))；
(b) 收尾是否漏取倒数——`bᵀ·pinv(H_loc)·b` 量纲是柔量 1/k，直接当 k 用
量纲不自洽（这也可能是作者未放出该方法的原因之一）。

**裁决**：采用量纲自洽的局部柔量（compliance，Decius 1963 系）解释：

    k_t = 1 / (bᵀ · pinv(H_loc, σ > 1e-8·σ_max) · b)

- 截断阈值用相对 1e-8（与 drih 决策 0012 同精神；Sobtop 的绝对 1e-10 是
  E_h/Bohr² 体系下的值，换算后不等价）；
- 刚体平移零空间由 SVD 截断处理；对 (i,j) 交换严格不变（无 mSeminario
  方向不对称）；线性角 → LinearAngleError（与 mseminario 一致）；
- 只做 bond/angle（柔性扭转走 MLP 管线，设计 §10.1）；非正 k → 
  NonPositiveForceConstantError（决策 0014 同款，不回退）；
- 乙醇 GFN2 烟测：键 2.7–4.2×10⁵ kJ/mol/nm²、角 490–1030 kJ/mol/rad²，
  与 mseminario 同量级偏刚（局部块忽略环境耦合，符合柔量法特性）；
  硬度排序（O-H > C-H > C-C）化学正确。

无 oracle ⇒ 标注 experimental（CLI choices 注释），不做 Gate 标准；
测试策略为内部一致性（B 行数值梯度、交换不变性、量级窗口、错误路径）。

## CP2K provider 的数据质量边界

格式（FeCl2 官方示例逐位核对）：`VIB| Hessian in cartesian coordinates`
后 5 列分页矩阵（页头列号 + 3N 行"行号 元素 值×5"），单位 E_h/Bohr²；
元素取自行序（每原子 3 行，按 r//3 折叠）；坐标取最后一次 QUICKSTEP
坐标表（Å），`positions_nm` 参数可覆盖（CLI 传管线的权威几何）。

**该文件打印的是前向差分 raw Jacobian**：与 CP2K 自己报告的低频/本征
向量不自洽（Rayleigh 商验证：打印矩阵 + 打印本征向量给不出打印频率，
-735 vs -2.2 cm⁻¹）；平移不变性破坏达对角元量级。Reader 读入时对称化
并把 raw 不对称度记入 meta；`check_hessian` 的平移检查会拒绝这类文件
（正确行为）。Sobtop 能用它是因为 Seminario 投影只取 3×3 局部块、对
全局低频噪声不敏感。**结论：解析职责 = 忠实打印值（已逐位验证）；数值
质量是文件属性，由既有质量检查把关。**夹具为裁剪版
`tests/fixtures/hessian/cp2k/FeCl2_5x5x1_freq.out`。

## IR 周期盒（box_nm）

- `ForceFieldIR.box_nm: 3×3 | None`（行 = 晶胞向量 nm），JSON 往返逐位；
- validation：形状 (3,3)、有限、det > 0；
- **NumPy 参考能量**（`ir/energy.py`）：box 非 None 时先沿键图 BFS unwrap
  （wrap 打散的坐标重整为连续分子坐标；对连续坐标恒等），非键/1-4 距离
  用最小镜像。wrap 不变性测试：`mod(pos, L)` 后能量逐位恢复；
- **GROMACS**：gro 盒行——正交 3 值、三斜 9 值（GROMACS 约定顺序）；
- **OpenMM**：`setDefaultPeriodicBoxVectors`（要求下三角形式，非下三角
  报错提示规范化）；
- **AMBER**：prmtop `IFBOX=1` + `BOX_DIMENSIONS`（beta=90、box x/y/z Å）
  + rst7 末行盒；**只支持正交盒**（prmtop 格式限制），三斜 → ValueError。
- 跨引擎验证（check-equivalence）暂不扩展周期腿（OpenMM 侧能量在无
  cutoff 下周期最镜像等价，后续按需加）。

## 测试

新增 33 个：localhess 7、CP2K 11、周期盒 8（含 wrap 不变性、三导出器、
validation、JSON 往返）。全量见 CHANGELOG。

## 关联

- 决策 0022（第一批吸收）、0011/0012（公式与单位）、0007（Hessian 质量
  检查）、0014（非正 k 处理）
- 报告：`docs/theory/sobtop_reverse_engineering.md` §5.1、§8
