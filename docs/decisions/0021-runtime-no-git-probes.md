# 决策 0021：运行记录禁止自动探测 Git

- 日期：2026-10-17。用户指出 runtime.py 的 RunContext 通过 `git -c safe.directory=…` 探测状态，违反本 NFS 工作区禁止 Git 的规则，要求在重跑 Gate 前移除。
- RunContext 删除 `_git_state` 及子进程探测；manifest 保留 `git_commit`、`git_dirty` 字段，值均为 null。不把未知状态记为干净提交。
- scripts/run_logged.sh 一并移除同类 Git 探测，START 行使用 `commit=null`。日志时间戳、被包装命令退出码与 END 行保持。
- 保留源码与输入清单，作为运行内容身份记录；本变更不改写历史运行产物。
- 回归：Python 探测路径由禁止创建子进程的拦截器保护；shell 包装器路径用 PATH 中的假 Git 留痕，避免测试失败时调用真实 Git。
- 2026-10-17，node-A：新增约束对旧实现为 2 failed；修复后 tests/unit/test_runtime.py 为 8 passed。未重跑 Gate。
