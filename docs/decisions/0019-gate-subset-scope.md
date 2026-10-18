# 决策 0019：Gate 子集运行只记 scope=subset，不能取得 passed=true

- 2026-10-11 用户同意（来源：`docs/superpowers/plans/2026-10-11-crucible-next-steps.md` §3 A1 的讨论）。
- 问题：G1 `--providers`、G2/G3 `--molecules` 只跑部分矩阵时，`acceptance.json` 仍可能是 `passed=true`
  （核实：`scripts/gate_mvp3.py:135-136` 按 `--molecules` 过滤预期集合，`:83` 不检查是否全集）。
  `manifest.debug_flags` 记了参数，但不能防止误读。
- 决定：
  1. 子集运行照常执行，不拒跑；`acceptance.json` 加 `scope`（`full` / `subset`）。
  2. `scope=subset` 时 `passed` 一律为 false，各判据结果照常写在 `criteria` 里供调试；report.md 写明"子集运行，不是正式通过"。
  3. 不新增数值阈值；现有容差与决策 0017 不变。
- 误拦代价：调试冒烟的退出码是 1、`passed=false`（`--skip-pytest` 已经是这样）；不影响任何正式运行。
