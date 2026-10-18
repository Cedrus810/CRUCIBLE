# 2026-10-17 全量 pytest（当前源码，去掉 Git 探测之后）

node-A，openmm_dev，`PYTHONPATH=src`，8 线程，不用 GPU。

| 部分 | 命令 | 结果 | 日志 |
|---|---|---|---|
| 非 external | `pytest -p no:cacheprovider -q -m "not external" --durations=5` | 502 passed，85 s | `unit.log` |
| external | `pytest -p no:cacheprovider -q -m external --durations=10`（setsid nohup 后台，PID 1357867，见 `pids.txt`，已结束） | 59 passed，546 s | `external.log` |

合计 561 个，0 失败、0 跳过。最慢的是 `test_cli.py::test_derive_torsions_mace_ethanol`（MACE CPU 扫描，415 s）。
314 条 warning 都是第三方库的（torch.load weights_only 等），和 10-13 Gate 的数量相同。
G0–G3 未重跑。
