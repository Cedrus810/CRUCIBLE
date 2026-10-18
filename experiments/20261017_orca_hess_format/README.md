# 2026-10-17 ORCA .hess 真实格式核对（路线图 D3）

目的：10-17 写的 `hessian/orca.py` 按猜测的格式实现（`$orca_hessian` 节、下三角、`$atoms` 每原子 4 个 token），
验收时怀疑与 ORCA 实际输出不符。用本机 ORCA 生成一个真实文件核对。

- 程序：`/home/user/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2/orca`（ORCA 6.1.1；
  注意 `~/.local/bin/orca` 是别的软件 orca-ide，不是量化程序）。node-A，单核，几秒。
- 输入：`water.inp`（HF/STO-3G Freq，未优化几何）。输出：`water.out`、`water.hess` 等。

## 实际格式（与原实现的假设全部不同）

| 项 | 原假设 | ORCA 6.1.1 实际 |
|---|---|---|
| 文件头 | `$orca_hessian` 节，内含原子数 + 数值 | `$orca_hessian_file`（无内容）；Hessian 在 `$hessian` 节 |
| `$hessian` 第一行 | 原子数 N | 维度 3N |
| 矩阵布局 | 下三角或连续满阵 | 按列分块的满阵：每块一行列号表头，再 3N 行"行号 + 该块各列值"（每块 5 列） |
| `$atoms` 每行 | 元素 x y z（4 个） | 元素 质量 x y z（5 个），坐标 Bohr，已平移到质心 |
| 其他节 | — | `$act_atom`、`$act_coord`、`$act_energy`、`$multiplicity`、`$vibrational_frequencies`、`$normal_modes`、`$actual_temperature`、`$frequency_scale_factor`、`$dipole_derivatives`、`$ir_spectrum`、`$end` |

原实现读这个文件时，第一步就会报"缺少 $orca_hessian 节——不是 ORCA .hess 文件"。

## 重写后的核对

用重写后的 `read_orca_hess` 读取 `water.hess`，再用 `$atoms` 里 ORCA 的质量，通过 `harmonic_frequencies` 算频率，与文件里 ORCA 给出的 `$vibrational_frequencies` 比较：

| 模式 | ORCA (cm⁻¹) | 我们 (cm⁻¹) | 相对差 |
|---|---|---|---|
| 6 | 2043.23419 | 2043.23431 | 5.737e-8 |
| 7 | 4488.35107 | 4488.35133 | 5.736e-8 |
| 8 | 4790.41182 | 4790.41210 | 5.736e-8 |

三个模式的相对偏差相同，都是 5.74e-8。这说明偏差来自两边的单位常数（Hartree、Bohr）略有不同，而不是布局或单位写错；布局或单位写错会给出 O(1) 的偏差。H 原子坐标 y = ±0.07572 nm，与输入一致。
`water.hess` 和 `water.inp` 已复制为测试固件 `tests/fixtures/hessian/orca/`（`test_real_orca_water_frequencies_match_orca`，容差为相对 1e-6）。

仍未覆盖：其他 ORCA 版本、开壳层、`$hessian` 每块不是 5 列的写法。读取器按表头列号读，不依赖每块的列数。
