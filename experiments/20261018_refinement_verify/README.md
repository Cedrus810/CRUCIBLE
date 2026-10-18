# 2026-10-18 可选顺序重拟合（refinement_rounds）验收验证

接收方：代理，node-A，openmm_dev，`PYTHONPATH=src`。只读 G2 扫描缓存，不加载 MACE、不调 xtb、不重扫。

| 检查 | 命令 / 脚本 | 结果 | 日志 |
|---|---|---|---|
| 交接 §5 相关测试 | `pytest tests/unit/test_torsion_refinement.py test_torsion_refit.py test_pipeline_chain.py test_cli.py test_staleness.py test_torsion_scan.py -m "not external"` | 66 passed（54 s） | `related_tests.log` |
| 3 轮 vs B3 机制实验 | `compare_with_b3.py`（10-13 G2 缓存，cache_only，参考势调用即报错） | 6 个分子项布局相同，k、相位、逐键 best/min-shift 与拟合残差的差全部为 0（逐位相同） | `compare_with_b3.log` |
| 默认 0 轮 vs 改动前单遍 | `compare_rounds0_with_g2.py`（10-18 G2 缓存与存档 IR，同在 node-A） | 6/6 个分子 proper_torsions 逐位相同；`fit_protocol` = single_pass/0 | `compare_rounds0_with_g2.log` |
| 非 external 全量 | `pytest -m "not external"` | 524 passed（117 s；原 502 + 新 22） | `unit_full.log` |

解析耦合用例的期望值（3 轮后 0.9859375 / 1.00703125）已按递推 a ← 1.5 − b/2、b ← 1.5 − a/2 从单遍 1.45 手算核对。
之后在新源码上重跑 G0–G3（`jobs/gate_mvp*_20261018_2/`）：全部 passed，每个 Gate 的 pytest 524 + 59（含 external）。
