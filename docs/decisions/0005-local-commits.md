# 决策 0005：恢复本地版本控制与分组提交

- 2026-10-04 用户明确要求继续施工并创建多组本地 commit，覆盖决策 0001 的暂缓指示。
- 工作区仍为 `/home/user/forcefield_compiler`，不新建工作区或另开对话。
- 实测共享盘 Git loose object 关闭时 Permission denied。按原施工计划，Git 元数据使用本机盘 `/home/user/.local/share/git-sep/ffcompiler.git`，工作区 `.git` 是指向它的文件。
- 失败的空仓库初始化保留为 `.git.nfs-initialization/`，忽略提交；没有覆盖已有用户提交。
- 原环境没有提交身份，仓库局部使用 `本机代理`，不改用户全局 Git 配置；所有提交仅在本地。
- 混合挂载所有权导致 safe.directory 检查失败，命令和项目运行器对当前路径使用一次性 `-c safe.directory=...`，不改全局安全目录。
- baseline 02acbd9 保存现有 MVP-0 交付；后续按 Gate、源力场语义、IR 校验、GROMACS、文档、验证证据分组。
