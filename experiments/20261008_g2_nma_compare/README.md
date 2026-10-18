# G2.3 nma compare=False 诊断（2026-10-08）

- `run.py mlp`：用 `jobs/gate_mvp2_20261008_3/ir/nma.json` 重跑 compare。唯一不过的腿是 sander-ours dihed：sander 1.8638，numpy 1.863570 kcal/mol，偏差 2.3e-4，容差 2e-4（`run_mlp.log`）。
- `run.py gaff2`：原 GAFF2 IR 全部通过（`run.log`）。
- `snap_check.py`：假设 sander 把 |δ−π| ≤ 4e-3 的相位改成 π。照此重算得 1.863826，打印为 1.8638，与 sander 一致；不吸附的值打印为 1.8636，不一致（`snap_check.log`）。被吸附的拟合项是 0-1-3-4，n=2，δ=3.141059，k=13.09 kJ/mol。
- 阈值 4e-3 是我凭记忆写的 sander 源码中的值，【未核实于源码】。
