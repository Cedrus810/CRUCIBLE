# 覆盖集 fallback 阶段 4 个虚频分子的诊断（2026-11-25）

来源：`jobs/coverage_20261124/` 中 acetone、methyl_acetate、styrene、glycine_zwitterion 的
fallback 阶段（fast preset：xtb GFN-FF `--ohess vtight` → mSeminario）被 `check_hessian`
以虚频 < −20 cm⁻¹ 拒绝。脚本 `diagnose.py`，结果 `results.json`，各分子重启运行在同名子目录。

## 检查

1. xtb.out 原文（acetone）："found 2 significant imaginary frequencies / writing imag mode distorted
   coords to xtbhess.xyz for further optimization"；acetone 鞍点对称性 c2v。
2. 虚频模式与各键刚性扭转位移的重叠 |cos|（g98.out 的笛卡尔位移）。
3. 区分性检查：从 xtb 写的 `xtbhess.xyz` 出发，用管线同一函数 `run_xtb_hessian(method="gfnff")`
   （--ohess vtight）重新优化，再过同一个 `check_hessian`。

## 结果

| 分子 | 原虚频 cm⁻¹ | 主扭转重叠 | 重启 − 鞍点 kJ/mol | 重启 check_hessian | 重启最低频率 cm⁻¹ | 对称性 |
|---|---|---|---|---|---|---|
| acetone | −56.1, −50.4 | C1–C2 0.71 / 0.66 | −0.66 | pass，0 虚频 | 19.7 | c2v → c2 |
| methyl_acetate | −62.8 | C2–C3 0.97 | −0.74 | pass，0 虚频 | 49.7 | cs |
| styrene | −41.2 | C1–C2 0.30（混合模式） | −0.18 | pass，0 虚频 | 49.2 | c1 |
| glycine_zwitterion | −221.9 | N0–C1 0.18（混合模式） | −9.42 | pass，0 虚频 | 50.0 | c1 |

## 结论（检查 3 可区分）

GFN-FF vtight 优化停在鞍点；沿 xtb 给出的虚频方向重启一次，四个分子都到达无虚频极小、能量更低、
通过同一质量门。问题不在 GFN-FF 势能面或 Hessian 本身，而在优化停在鞍点且管线没有重启
（`run_xtb_hessian` 不处理 xtbhess.xyz）。

## 未确认（Claude 推断）

- 停在鞍点的原因推测为起点（MMFF 几何）接近对称构象、对称破缺方向梯度为零；未验证起点对称性。
- styrene、glycine 的虚频不是单键纯扭转，模式性质未定。
