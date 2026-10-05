# 决策 0013：mixed 精度 GROMACS 的验收容差；乙腈移出 Gate（取代 0006）

- 2026-10-05 用户裁定："错在精度，那不应该成为门的问题。"用户不使用近线性角
  体系，乙腈（C–C≡N）在常规 MD 中本身也难处理，不作为验收对象。
- **证据：失败来自精度，不是拓扑错误**（`jobs/gate_mvp0_20261004_3`）：
  - 同一份 `write_gromacs` 写出的 top/itp：openmm-top 腿（fp64 读取）288 行 0 失败，
    gmx 腿（mixed，fp32）432 行 11 失败。
  - 10 个失败行（bond/total）相对偏差 4.5e-7 – 2.7e-6，即 4–23 ε₃₂。
  - 拆分实验（乙醇 cfg1 bond）：只把参数压成 fp32 偏差 +2.1e-5；坐标、参数、运算
    全 fp32 偏差 −8.5e-5；GROMACS 实测 −1.7e-4。主因是 fp32 坐标/运算在 (r−r₀)
    相减中的放大（该构型 |r−r₀|/r 最小 5.8e-4），参数压缩只占一成左右。
  - 第 11 行是乙腈近线性角（θ = 179.98°）的 angle，2.4e-3 相对：acos 在 θ→180°
    病态，与 0004.2 中 sander 的近线性角问题同类。
- **决定**：
  1. gmx 腿按 `precision` 选容差（`validation/engines.py::gmx_energy_tolerance`）：
     - mixed：`1e-4 + 1e-5·|E|` kJ/mol（≈84 ε₃₂，对实测最大 2.7e-6 留 3.7 倍余量）；
       total 行的相对项按最大分量量级（同 0004.1）。
     - double 或精度未知：原经文件容差 `1e-4 + 1e-7·|E|` 不变。
     - `evaluation.json` 的 `tolerance_policy` 记录实际使用的政策。
  2. 乙腈移出 G0、G1 的评估矩阵（12 → 11 个固件），改为诊断用例
     （`test_acetonitrile_diagnostic`：各引擎须跑通，超差行只打印）。固件保留。
  3. 线性键角默认策略仍是报错（`LinearAngleError`）：真遇到线性角会明确失败，
     不会静默给出参数。`library` 策略代码与单元测试保留。
- **判别能力**：真实转换错误 ≥1e-3 相对（0004.1），对 mixed 容差仍有约 100 倍余量。
  新增 `test_compare_detects_tampered_gromacs_topology`：bond k×1.001 只写进 itp，
  gmx 腿的 bond 行必须失败；同几何未篡改时 gmx 行全部通过。实测通过。
- 不改变：其余各腿容差、OpenMM/sander 判据、气相单点设置；旧的失败运行保留不覆盖。
