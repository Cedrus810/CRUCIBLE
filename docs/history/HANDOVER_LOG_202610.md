# HANDOVER 历史日志（2026-10-10 → 2026-10-17）

2026-10-17 文档整理时从 `HANDOVER.md` 原样移出（逐日"做了什么"、旧草稿说明与旧待办），只作历史；当前状态以 `HANDOVER.md` 为准。文中"今天""本文件"指当时的 HANDOVER。

## 2026-10-17 晚：接收方验收与修复（代理，node-A）

用户："我往里塞了一些新的代码 验收一下" → "交接给你 你负责修复"。
1. 按设计稿逐条验收，结论不通过：新测试首跑 10/43 失败（2 个源码 bug + 8 个测试错误），代码阅读另发现 7 处问题。
2. 全部修复（清单见施工计划 §0）：check-staleness 必崩、staleness 比错层级、ORCA 读取器按真实文件重写
   （`experiments/20261017_orca_hess_format/`：ORCA 6.1.1 水 HF/STO-3G，频率与 ORCA 相对差 5.7e-8）、
   check-equivalence --ir-in 不挂 tleap 腿、坐标 sidecar 升 ffcoords/0.2（带原子名）、MM 快照补 improper、
   learned_bcc base 必填、下界诊断 IR 改名、DFT 文件核对原子序数。
3. 验证（node-A）：新/改测试文件 49 passed；`-m "not external"` 502 passed（78 s）；external 的
   check-equivalence 与 derive-bonded gfnff 2 passed。**未跑**：MACE 的 `test_derive_torsions_mace_ethanol`、其余 external（与改动无关）。
4. 注意：ORCA 量化程序在 `/home/user/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2/orca`；`~/.local/bin/orca` 是 orca-ide，不是它。

## 2026-10-17 做了什么（代理，登录节点，只写未跑）

用户指令："没搭建的 你搭建 但是 然后别人负责交接"；过程中两次喊停登录节点上的
测试运行（"大哥没叫你 在登陆节点跑""不要乱动"），裁定分工：**登录节点只写
代码+静态验证；测试/Gate/一切计算由接收方在计算节点跑**（该规矩已入全局记忆）。

1. G1.8 收口：`docs/theory/hessian_methods.md` 四处修正（§4.3 补平方、§6
   rad⁻²、§5.3 措辞+凸性、§2 方向记号），状态改已审阅。
2. 新增 `hessian/orca.py`（D3）：.hess 读取器，**尚无真实 ORCA 样例核对**；
   `tests/unit/test_hessian_orca.py` 9/9 在喊停前已过（纯解析）。
3. C1 坐标合同：`ir/coords.py` sidecar（ffcoords/0.1）+ `cli.py` 重构
   （--ir-in/--coords-in、derive-bonded 出 ir/coords/report、convert/
   check-equivalence 吃最终 IR、参数组合 exit 2 语义保持）。
4. C3/C4：`bonded/staleness.py`（H_fixed 与 MM 扣除依赖快照+过期检查）；
   global_fit 记 `h_fixed_dependency`；mlp `_finish` 记 `mm_dependency`；
   CLI 新增 `global-fit`（下界命中→status=failed_lower_bound、不导出、
   exit 1）与 `check-staleness`。
5. MVP-5：`hessian/fallback.py`（gfnff→gfn2→dft 文件级，check_hessian
   质量门，Tier 注入；CLI 不暴露 auto）。
6. MVP-4：`charge/learned_bcc.py` 骨架（§13 形式+符号约定；工件缺失即
   报错不降级）+ 立项草案 `2026-10-17-mvp4-learned-bcc-proposal.md`
   （S1–S6 科学决策全部留白，批准前禁止训练）。
7. 新测试文件 5 个（coords/staleness/fallback/learned_bcc/pipeline_chain，
   共 34 用例）+ ORCA 9 例，全部只写未跑；静态验证：py_compile 全过、
   `--collect-only` 555 collected（512+43）、CLI 接线检查过。

**接收方任务（按序）**：
1. ~~全量 pytest~~ 2026-10-17 晚已做：非 external 502 passed（见上一节）；MACE 用例未跑。
2. ~~ORCA 真实样例核对~~ 2026-10-17 晚已做并重写读取器。
3. B1/B2 诊断与 C 的正式外部组合验收（10-15 计划 §4–§5，需报资源）。
4. MVP-4 立项讨论（读立项草案）；G1.8 收口后的文档状态确认。
5. 版本号与提交：`.commit_drafts/20261017_*.md` 6 份草稿（见下）。

## 2026-10-13 做了什么

用户给的顺序：先收 A1/A2 的两个运行细节，再做 A3（重跑 G0–G3、G1.8 审阅），之后 B（苯甲酸甲酯）/ C（串联管线）。
1. 决策 0020（用户选定）并实现：`refit_torsions` 只重拟合原 mlp 键（`MLPTorsionProvider(fit_bonds=...)`），
   原兜底键原样沿用；G2 `--resume DIR` 没有 `DIR/scans/` 时参数报错退出，0 文件时 warning；
   G3 `--scan-cache` 必填且须是目录（删掉指向 `gate_mvp2_20261008_3/scans` 的默认值）。
2. 新增 3 个用例（RED 已复现升级为 mlp）；相关 8 个测试文件 not external 80 passed（node-A，56 s）。
3. A4 逐条核实（只读）：结果表写在 next-steps §A4；没有会改变 Gate 结果或参数的条目，随 B/C 处理。
4. A3：用户同意后重跑 G0–G3，全部 passed（见“正式运行”）；按设计稿逐条验收，并完成 G1.8 审阅，写在 `docs/history/ACCEPTANCE_20261013_A3.md`。
   新发现的偏离：全局拟合没有把 H_fixed 的电荷/二面角来源记进 provenance（设计 §6.4.2），归到 C。
   注意：node-B 的 NFS 属性缓存很长，在 node-A 上看 node-B 写的 log 可能是旧的，要到 node-B 上直接读。

## 2026-10-11 做了什么

1. 整理后续工作顺序（`2026-10-11-crucible-next-steps.md`，A–F）；用户同意先做 A1/A2。
   核实：G0/G1 记录的源码 9/75 已变（含 `units.py`、`bonded/symmetry.py`），G2 5/91，G3 0/104。
2. 决策 0019：Gate 子集运行记 `scope=subset`，`passed` 不为 true。
3. A1/A2 施工（subagent-driven，T38–T42，每个任务独立审阅 + Opus 整体审阅 Ready）：
   - T38 G1/G2/G3 子集不能通过；T39 G2 `--resume` 复制旧 scans 到本次目录，旧目录只读；
   - T40 `assign_charges` 检查电荷形状/有限性，`provenance["charge"]["previous"]` 保留上一次来源；
   - T41 重拟合沿用原兜底原因（不再改记 `no_cached_scan`），不改写缓存；
   - T42 `--reference mm` 缓存键与 torsion provenance 记参考势身份（扫描时的电荷，字面值），MACE 键不变。
   - 验证：相关单元测试 77 passed（not external）；用现行代码重算 `gate_mvp2_20261008_3/scans` 12 根键的缓存键，12/12 与存档一致。
   - **正式 Gate 未重跑**（属 A3，等用户提交后另报计划）。

## 2026-10-10 做了什么

1. **文档整理**：三语 README、roadmap（加状态列、D1–D6 结论、过时事实补注）同步到实际状态；
   用户确认 `forcefield_compiler_design.md` 是现行设计稿，`.v0.md` 标"已作废"；MVP-1 交付清单移到
   `docs/history/DELIVERY_20261005.md`；删了根目录 sander 残留 `mdinfo`、`restrt`。
2. **G2 通过**：用户决定"OpenMM 说没问题就没问题"，G2.3 的 sander 腿改为只报告（决策 0017；
   sander 把接近 π 的拟合相位吸附成 π，AMBER 那一路留给以后测）。Gate 加 `judge_compare`，
   在 node-B 上 `--resume` 重跑 → `jobs/gate_mvp2_20261010_2/` passed。
3. **已知缺陷**：methyl_benzoate 1-2 的 MLP 二面角剖面比 GAFF2 差（3.36 vs 2.70 kJ/mol），用户定为
   已知缺陷、之后再想，写进三语 README 的已知限制。原因未查。
4. **MVP-3 电荷 provider**：讨论定范围（决策 0018：RESP 往后放；载入 / AM1-BCC / NAGL；先电荷后二面角，
   换电荷用缓存的扫描重拟合）→ 冒烟 `experiments/20261010_am1bcc_repro/`（AM1-BCC 12/12 与固件逐位相同）→
   计划 → Native 实现 T31–T37 → 整体审阅（2 个 Important 已修，9 个 Minor 推迟）→ 在 node-A 上跑
   `jobs/gate_mvp3_20261010_2/` passed。

## 待用户提交的草稿（`.commit_drafts/`，按顺序）

1. `20261008_*`：10-08 的改名、T24–T30、评审修复、G2 首跑（共 10 份，是否已提交未知）。
2. `20261010_1200_docs-tidy.md`：文档整理。
3. `20261010_1715_g2-openmm-reference.md`：决策 0017、Gate G2.3 判据、G2 通过、methyl_benzoate 缺陷。
4. `20261010_1800_mvp3-plan.md`：决策 0018、MVP-3 计划、AM1-BCC/NAGL 冒烟。
5. `20261010_1945_mvp3-charge-providers.md`：MVP-3 代码、测试、Gate G3、文档。
6. `20261011_2330_a1a2-gate-scope-provenance.md`：next-steps 整理、决策 0019、A1/A2（T38–T42）代码与测试、文档。
7. `20261013_1340_refit-fallback-cache-paths.md`：决策 0020 的代码、测试、文档。
8. `20261013_1530_a3-acceptance.md`：A4 核实、A3 验收报告、三语 README / HANDOVER / CHANGELOG 状态更新、启动脚本。
9. `20261017_1210_g18-doc-closure.md`：G1.8 四处修正与状态收口。
10. `20261017_1240_orca-hess-reader.md`：ORCA .hess 读取器与测试（D3）。
11. `20261017_1500_c1-coords-contract.md`：坐标 sidecar、CLI 重构、BondedReport 序列化。
12. `20261017_1540_c3-c4-staleness-globalfit-cli.md`：依赖快照/staleness/global-fit/check-staleness。
13. `20261017_1620_mvp5-fallback.md`：Hessian fallback 框架与测试。
14. `20261017_1650_mvp4-learned-bcc-skeleton.md`：learned BCC 骨架、立项草案、测试。
15. `20261017_1720_pipeline-tests-and-handover.md`：管线串联/组合测试、施工计划、台账、HANDOVER 本章。
16. `20261017_2000_unbuilt-modules-accepted-and-fixed.md`：**取代 10–15**；10-17 施工 + 当晚验收修复，一个提交。

**更新（2026-10-17 晚，用户确认修复后）：草稿 10–15 已作废，由
`20261017_2000_unbuilt-modules-accepted-and-fixed.md` 取代（一个提交，含 10-17 施工 + 验收修复 +
本文件 + CHANGELOG）；草稿 9（G1.8 文档）仍有效、与之不重叠。下面 10–15 两条旧说明只作历史。**

**2026-10-17 的草稿（9–15）在测试运行之前提交的话，提交信息里注明"测试未跑
（登录节点禁跑），验证待计算节点"——或者等接收方跑绿后再提交，由用户定。**

几份草稿改的是同一批文件（README、HANDOVER、roadmap）的工作区状态：先执行的那条会把后面的改动一起带上，后面的命令仍可执行。也可以把 2–5 合成一个提交。**本文件（HANDOVER）的这次重写不在任何草稿的说明里**，提交时一起 add 即可。

## 待办

**等用户定：**
1. MVP-2 计划 §7 第 2–4 条：扫描几何变形是否作兜底条件；同一键能量并到一个四元组的近似；`spin=1` 的含义。
2. methyl_benzoate 1-2 已知缺陷（§7 第 6 条）：迭代拟合，还是换剖面指标。候选原因：逐键拟合时其余键仍是 GAFF2；"各减最小值"的指标对单点敏感。都未核实。
3. MVP-3 推迟的 Minor：电荷来源覆盖、子集 Gate 通过、兜底原因改写、zip 非 strict、重拟合弱断言、
   `--reference mm` 缓存键不含电荷 —— 已在 A1/A2（T38–T42）修复。其余（`--scan-cache` 提示已随决策 0020 做完；重复
   `find_rotors`、mol2 报错）属 A4，未做。
3b. A1/A2 最终审阅留下的三条问题：2026-10-13 按决策 0020 已做。
4. MVP-4（learned BCC）立项讨论。

**人工 / 需要先提交：**

5. G1.8 理论文档审阅：2026-10-13 已审（验收报告 §2），实质内容通过；有 4 处措辞要修（§4.3 少平方、§6 rad² 笔误、§5.3 一句不通、§2 u_AB 方向记号），修完把文档开头的状态改为已审阅。
6. A3 重跑 G0–G3：2026-10-13 已做，全部 passed（见“正式运行”）。`git_commit` 为 null，提交后可以在实际机器上把哈希补记进验收报告。
7. 版本号（0.2.0）要不要改，等用户定。

**已知小问题（不阻塞）：**
- MVP-2：（`--resume` 写旧目录已由 T39 修复）其余见 `.superpowers/sdd/2026-10-08-ffc-mvp2-mlp-torsion/progress.md`。
- 非 external 的 pytest 在 node-A 高负载时要 72–131 s，超过 MVP-0 §1.5 的"1 分钟"约定；在 node-B 空闲时为 63 s（G2 实测）。
