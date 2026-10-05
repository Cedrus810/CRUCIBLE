# 决策 0001：版本控制暂缓（用户指示"不管 git"）

> **已被决策 0005 取代**（2026-10-04）：用户要求恢复本地版本控制与分组提交。下文保留为历史记录。

- 日期：2026-10-04
- 决策人：用户（开工指令"准备开工 但是不管git"）
- 背景：MVP-0 计划 T0 要求 `git init --separate-git-dir=/home/user/.local/share/git-sep/ffcompiler.git`（NFS 上 git 写对象出权限错，git 目录放本地盘）。2026-10-04 实测 `/home/user` 在本机不存在，且用户随后明确指示开工但不管 git。
- 决定：
  1. 不执行 `git init`，不做任何提交；T0 的"任务完成即提交"要求同步豁免。
  2. `RunContext` 的 manifest 字段 `git_commit`、`git_dirty` 记 `null`（实现里保留探测逻辑，将来有仓库自动生效）。
  3. `scripts/run_logged.sh` 的 START 行 commit 显示 `unknown`（与上游 Amber-X 脚本的 fallback 行为一致）。
- 影响：Gate G0.6 涉及"偏离须有决策记录"——本文件即该偏离的记录。后续若要恢复版本控制，按原计划补 `--separate-git-dir` 即可，代码无需改动。
