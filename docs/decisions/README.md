# 决策记录索引 / Decision records

每条记录说明一次偏离计划或补全定义的原因与证据。被取代的记录保留原文作为历史。
Each record explains why the project deviated from, or completed, a plan definition.
Superseded records are kept verbatim as history.

| # | 主题 | 状态 |
|---|---|---|
| [0001](0001-git-layout.md) | 版本控制暂缓 | 已被 0005 取代 |
| [0002](0002-improper-order.md) | improper 的原子顺序与连通性校验 | 有效 |
| [0003](0003-gromacs-improper-funct4.md) | GROMACS funct 4 是周期 improper（勘误） | 有效 |
| [0004](0004-cross-engine-tolerances.md) | 跨引擎容差的数值地板；sander 近线性角系统差 | 有效 |
| [0005](0005-local-commits.md) | 本地版本控制与分组提交，git 元数据放本机盘 | 有效 |
| [0006](0006-gromacs-strict-precision.md) | GROMACS 继续用严格文件容差 | 已被 0013 取代 |
| [0007](0007-hessian-asymmetry-threshold.md) | Hessian 不对称度阈值取 §2.3 表（1e-3 报错 / 1e-6 警告） | 有效 |
| [0008](0008-write-all-code-before-review.md) | 先写完 MVP-1 代码，T15 审阅移到交付之后 | 有效（mSeminario 部分由 0011 结案） |
| [0009](0009-symseminario-naming.md) | m2Seminario 暂名 symseminario | 已被 0011 取代 |
| [0010](0010-drih-ring-mismatch.md) | DRIH 环系与 Sobtop 不一致 | 已被 0012 取代 |
| [0011](0011-seminario-formulas-and-m2seminario-naming.md) | mSeminario / m2Seminario 公式经 Sobtop 12/12 核对；方法名 m2seminario | 有效 |
| [0012](0012-drih-bohr-length-unit.md) | DRIH 内部使用原子单位（Bohr、E_h） | 有效（用户已确认） |
| [0013](0013-mixed-precision-gmx-tolerance.md) | mixed 精度 GROMACS 用 fp32 容差；乙腈移出 Gate | 有效 |
| [0014](0014-method-limitation-nonpositive-k.md) | 非正力常数经 Sobtop 复算确认 → 已知方法局限 | 有效 |
| [0015](0015-native-amber-io.md) | AMBER prmtop / rst7 原生读写，去掉 ParmEd | 有效 |
| [0016](0016-mvp2-torsion-decisions.md) | MVP-2 前置决定：代理接手施工；AMBER 周期二面角、GAFF2 兜底、暂用现有电荷；MACE-POLAR-1-M 直接经 ASE | 有效 |
| [0017](0017-openmm-reference-sander-report-only.md) | MLP 二面角 IR 以 OpenMM 为准；G2 的 sander 腿只报告（相位吸附为已知局限） | 有效 |
| [0018](0018-mvp3-charge-scope.md) | MVP-3 范围：RESP 往后放；载入 / AM1-BCC / NAGL；先电荷后二面角，换电荷用缓存重拟合 | 有效 |
| [0019](0019-gate-subset-scope.md) | Gate 子集运行记 scope=subset，passed 不为 true | 有效 |
| [0020](0020-refit-fallback-and-cache-paths.md) | 重拟合不升级原兜底键；G2 `--resume` / G3 `--scan-cache` 路径写错时报错，G3 无默认缓存 | 有效 |
| [0021](0021-runtime-no-git-probes.md) | RunContext 与 run_logged.sh 不探测 Git；提交与 dirty 状态记 null | 有效 |

新记录用下一个编号，文件名 `NNNN-short-topic.md`；取代旧记录时在旧记录标题下加一行说明。
New records take the next number; when superseding, add a one-line note under the old title.
