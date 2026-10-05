# 决策 0002：improper 的原子顺序与连通性校验规则

- 日期：2026-10-04
- 背景：MVP-0 计划 §2.1 断言"AMBER 约定中心原子在第三位"，T2 校验项 2 据此要求 improper (i,j,k,l) 的第三个原子 k 与 i、j、l 都成键；T5 同时要求"原子顺序保持 prmtop 原样，不重排"。
- 实测（2026-10-04，tleap/AmberTools 生成、`tests/fixtures/gaff2/`）：
  - 苯的 improper 条目形如 `(H6, C0, C5, C1)`——中心原子 C0 在**第二位**；
  - NMA 的 improper 条目形如 `(C0, N3, C1, O2)`——中心原子 C1 在**第三位**。
  - 顺序随 GAFF improper 类型在 parm 文件中的定义顺序而定，并非统一"中心在第三位"。
- 二面角能量依赖四元组的**有序**次序：把中心原子移到第三位会改变 φ 与 E = k[1+cos(nφ−δ)]，导致导入后的 IR 与 sander / OpenMM 读同一 prmtop 的能量不一致，破坏 T5 用例 7 与 §19 的跨引擎等价。
- 决定：
  1. **导入/导出一律保持 prmtop（或 IR）中 improper 的存储顺序**，能量按该顺序求值（与 OpenMM `AmberPrmtopFile`、sander 对同一文件的求值一致）。计划的"不重排"照旧执行。
  2. T2 校验项 2 的 improper 连通性规则改为：**四元组中存在（且唯一）一个原子与其余三个原子都成键**（中心原子可以在任何位置）。其余校验项不变。
  3. T5 用例 3 相应改为：苯有 6 个 improper，每个四元组的连通中心是芳香 C（atom_type `ca`）。
- 影响：§2.1 中"AMBER 约定中心原子在第三位"一句在本实现中按上述实测修正；provenance、导出格式（GROMACS funct 4、ParmEd improper 标记）均按存储顺序处理，三后端能量一致性由 T9 的跨引擎测试兜底。
