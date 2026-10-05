# 决策 0012：DRIH 内部使用原子单位（Bohr、E_h）（取代 0010；用户已确认）

- 决策 0010 记录"DRIH 与 Sobtop 在环系上不一致，疑似未公开处理"，并称已排除
  单位换算——实际只排查了输出单位，没有排查 B 内部的长度单位。
- 机理：bond 行的 B 元素无量纲，angle/dihedral 行是 rad/长度。换长度单位等于
  对 B 的行做非均匀缩放。冗余关系只在同类坐标间时（无环分子）伪逆投影与缩放
  无关；冗余同时牵涉 bond 与 angle 时（环系）伪逆结果依赖单位。
- 实测（§2.8 容差，k/键长/键角，除 sf6_fchk）：

  | B 的长度单位 | 一致 case 数 |
  |---|---|
  | nm（原实现） | 5/12（仅无环） |
  | Å | 5/12 |
  | Bohr | **12/12**（最大相对 k 偏差 1.3e-5） |

- 处理：`bonded/drih.py` 把几何换成 Bohr、H 除以换算系数²后构造 B、G、A，
  bond k 再换回 kJ/mol/nm²。坐标集合、伪逆阈值 1e-8·λmax、A = G⁻B、取对角元
  均未改动。
- 性质：计划 §2.5 没有规定 B 的长度单位；这是**补全**未写明的定义，不是改定义
  凑数。
- **用户确认（2026-10-05）**：DRIH 内部约定为 Cartesian 坐标统一用 Bohr、Cartesian
  Hessian 统一用 E_h/a₀²（Gaussian/ORCA Hessian 本身就是原子单位体系）；H_q 的单位
  随坐标类型而定——stretch E_h/a₀²，angle/dihedral E_h/rad²，耦合项为混合单位。
  换到 IR 时 bond ×HARTREE/BOHR²（≈9.37583×10⁵ kJ mol⁻¹ nm⁻²），angle ×HARTREE
  （2625.50 kJ mol⁻¹ rad⁻²），**angle 不再除 a₀²**。不经过 Å。
- 实现按此字面执行（`bonded/drih.py`：入口换成 Bohr 与 E_h/a₀²，出口按坐标类型换回）。
  能量单位不影响 G 的伪逆，数值与只换长度单位的版本完全相同（oracle 仍 12/12）。
  理论文档 §4.1a 写入单位表与对应英文表述。
- 同时核对 ½k 约定：IR 统一按 E = ½kΔq² 储存；导出 GROMACS/OpenMM 用 k，AMBER 用
  K = k/2（`exporters/amber.py` 的 `_K_BOND = 2×4.184×100`）。代码一直正确，但理论文档
  §6 的 AMBER 换算漏了 ×2、SF6 核对值写成 700.28，已改正为 350.02 kcal/mol/Å²。
- 同批修正：drih 遇到输入里的线性键角改为抛 `LinearAngleError`（§2.4），经过
  线性键角的二面角按几何剔除（不再依赖 angle 列表，修复 library 策略下乙腈
  必然失败的问题）。
