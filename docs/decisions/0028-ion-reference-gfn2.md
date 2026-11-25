# 决策 0028：离子与两性离子的二面角参考势默认用 xtb GFN2（部分取代 0016 修订）

- 2026-11-25 用户同意。Claude 问"离子和两性离子的参考势改用 xtb GFN2、不再用 MACE，你同意吗？"
  （即 2026-11-20 计划 §7 待拍板 2），用户原话："可以"。
- **背景**：决策 0016 修订（2026-10-08）因 MACE-POLAR 有总电荷输入而让离子照常用 MACE 扫；但 MACE 对离子
  的适用性一直是在案的开放问题（A4 剩余项"甲铵正离子 MACE 回归"未做），覆盖集诊断把它记为缺项
  `charged_species(...; mace_applicability_open)`。
- **决定**：
  1. 判定：净电荷 ≠ 0，或任一原子带形式电荷（两性离子）→ 二面角扫描参考势用 xtb GFN2；其余仍用 MACE。
     实现 `torsion/xtb_ref.py`（`is_ion_or_zwitterion`、ASE 计算器 `XtbGFN2Calculator`：每次求能量/力调用
     `xtb --gfn 2 --grad --chrg Q --uhf 0`；openmm_dev 无 tblite/xtb-python，生产环境不装包）。
  2. CLI `derive-torsions --reference` 新增 `auto`（默认）与 `xtb-gfn2`；`auto` 按第 1 条选择并在 stderr 说明；
     `mace` / `xtb-gfn2` 显式强制。单命令 `run` 用 `auto`。中性分子行为与缓存键不变。
  3. 身份：模型名 `xtb-gfn2`，"模型文件"记 xtb 可执行文件的路径/大小/修改时间（同 MACE 模型文件记法，不用
     hash），进缓存键与 provenance；换电荷 refit 只读缓存时据此复原，不调用参考势。
  4. 覆盖集驱动：`refcheck` 对离子/两性离子记 `reference: xtb-gfn2`、按 GFN2 元素范围（H–Rn）检查，不再记
     `mace_applicability_open`；`scan_probe` 同一规则选参考势。
- **不变**：Gate G2/G3 直接构造 MACE（不经 CLI 默认），判据与缓存不变（acetate、methylammonium 在 G2 中仍是
  MACE 扫描）；中性分子一切照旧。
- **xtb 工作文件**放系统临时目录（`TMPDIR` 可覆盖），计算器回收时删除——是逐步瞬时文件，不是产物。
