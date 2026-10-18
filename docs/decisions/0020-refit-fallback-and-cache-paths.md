# 决策 0020：重拟合不升级原兜底键；--resume / --scan-cache 路径写错时报错

- 2026-10-13 用户选定（来源：A1/A2 最终审阅的推迟项，`.superpowers/sdd/2026-10-11-ffc-a1a2/progress.md` `Final: minor (deferred)`）。
- 决定：
  1. **重拟合只重做原来是 mlp 的键。** 原兜底键（`mlp_error`、`scan_not_converged` 等）原样沿用原行和原 GAFF2 项；
     缓存里后来补齐了匹配扫描也不升级为 mlp。换电荷只改变电荷，不改变哪些键被覆盖。
     实现：`MLPTorsionProvider(fit_bonds=...)`，`refit_torsions` 传入原 mlp 键（`src/ffcompiler/torsion/mlp.py`）。
  2. **G2 `--resume DIR`：`DIR/scans/` 不是目录就在解析参数时报错退出**（退出码 2，不建运行目录）；
     目录在但所选分子没有扫描文件时只打 warning（子集运行可能合法）。
  3. **G3 `--scan-cache` 必须显式给出且是目录**，去掉指向 `gate_mvp2_20261008_3/scans` 的默认值；
     正式 G3 指向同一轮 G2 运行目录下的 `scans/`。
- 误拦代价：2、3 只拦"给了路径但路径不存在"；要全新扫描就去掉 `--resume`，重新执行一次命令即可。
