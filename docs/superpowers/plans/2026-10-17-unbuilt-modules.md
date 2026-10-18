# 未搭建模块施工：ORCA 读取器、C1–C5、MVP-5 框架、MVP-4 骨架

日期：2026-10-17，login-node（登录节点，**只写代码 + 静态验证，未运行任何测试体**）。
依据：用户 2026-10-17 指示"没搭建的 你搭建 但是 然后别人负责交接"，随后
"你写交接文档 别人收到 继续跑"。规格来源：`2026-10-15-crucible-bc-detailed-next-steps.md`
（下称 10-15 计划）、设计稿 §8.1/§13–16/§23、路线图 D3。

**验证状态**：`py_compile` 全过；`pytest --collect-only` 555 collected（原
512 + 新增 43）；CLI 五个子命令 `--help` 与参数组合校验（coords-in 无
ir-in → exit 2）静态检查过。**全部测试体未运行**（登录节点禁跑，见
progress 台账）；运行归接收方。

## 0. 接收方验收与修复（2026-10-17 晚，node-A）

按设计稿逐条验收（第 3 章以下是施工方原文，**下列各处以本节为准**）。新测试首次运行
10/43 失败；修复后新/改测试文件 49 passed，全量 `-m "not external"` 502 passed
（78 s），相关 external 2 例（check-equivalence、derive-bonded gfnff）passed。

| # | 问题 | 处理 |
|---|---|---|
| 1 | `check-staleness` 必崩（`_validate_common_inputs` 读不存在的 `args.coords_in`） | 只对挂了公共输入的子命令校验 |
| 2 | `staleness` 拿整个 `h_fixed_dependency`（inputs/hessian/fd_step_nm）比快照 → global_fit 刚拟合完就 12 项过期 | 比 `dep["inputs"]` |
| 3 | ORCA 读取器格式全凭猜测，读不了真实文件 | 用 ORCA 6.1.1 真实输出重写（`experiments/20261017_orca_hess_format/`），频率与 ORCA 相对差 5.7e-8 |
| 4 | `check-equivalence --ir-in` 给齐三件套就挂 sander-tleap 腿，拿原 GAFF2 prmtop 比改过参数的 IR → 必然虚报失败 | `--ir-in` 模式不挂该腿（stderr 注明） |
| 5 | 坐标 sidecar 只记 (Z, 元素)，同元素原子换序查不出；三件套不全时交叉检查静默跳过 | sidecar 升 `ffcoords/0.2`，带原子名；只给 `--prmtop` 时 stderr 注明未做交叉检查 |
| 6 | MM 扣除快照漏了 improper | 补上；其他键 proper 有意不记（B2 的单遍近似），写进模块说明 |
| 7 | learned_bcc `base_provider=None` 取 IR 当前电荷 → 重复赋电荷叠加修正 | base 必填（S1 未定，不设默认） |
| 8 | global-fit 下界命中的诊断 IR 用标准名 `ir.json`；DFT 级读外部文件不核对分子 | 改名 `ir_failed_lower_bound.json`；核对原子序数序列 |

测试本身的错误（已改）：coords 测试辅助函数原子数不配；合成 Hessian 无二面角刚度却断言最低频 > 0；
下界用例的 fchk 文件不存在；错拓扑用例断言的报错文本不对。

**仍未做**：MACE 的 `test_derive_torsions_mace_ethanol`（external/mlp，未跑）；C2 串联测试里没有
torsion 一步；`derive-bonded --method` 仍无默认值（A3 #8，接口建议待用户定）；外部组合验收（10-15 §5 C5）。
（2026-10-18 更新：MACE 用例已在全量 pytest 中通过；C2 串联测试已补 torsion 一步，`test_local_route_chain_with_torsions`。）

## 1. ORCA .hess 读取器（路线图 D3）

> 2026-10-17 晚：本节的格式假设已证伪并重写，见 §0 第 3 条。

- 新增 `src/ffcompiler/hessian/orca.py`：`read_orca_hess(path)`。
- 格式假设：`$orca_hessian` 节（原子数 + Hessian 数值，Ha/Bohr²，下三角
  行主序或满阵按个数区分）、`$atoms` 节（元素符号或原子序数 + Bohr 坐标）、
  未知节跳过、D 幂宽容。梯度/能量不在 .hess 里（None）。
- **明确边界：尚未与 ORCA 实际产生的 .hess 逐值核对**（D3 原定等用户
  样例）。接收方拿到真实样例先跑 `tests/unit/test_hessian_orca.py` 的
  读取器本体核对；不符停报，不得改格式假设凑数。
- 测试 `tests/unit/test_hessian_orca.py`（9 个用例，2026-10-17 停止指令
  前已运行通过；纯解析无计算）。

## 2. C1：IR 与坐标传递合同

- 新增 `src/ffcompiler/ir/coords.py`：`Coordinates` sidecar
  （schema `ffcoords/0.1`，单位 nm，`atom_order` 逐原子 (Z, element)，
  `source` 来源；浮点与 IR JSON 同约定逐位往返）。
- CLI（`cli.py` 重构）：
  - `--coords-in` 只与 `--ir-in` 同用（argparse 报错 exit 2）；与 IR 显式
    核对原子序数/元素序列，不一致拒绝、不自动重排、不按同名配对；
  - `--ir-in` 模式下 prmtop/inpcrd/sdf 可选（三件套给齐则做拓扑交叉检查，
    沿用名字/元素/键集合规则）；无 `--ir-in` 时三件套必需（exit 2 语义保持）；
  - `derive-bonded` 补 `--ir-in/--coords-in/--ir-json/--coords-json`，
    输出 `ir.json + coords.json（hessian 几何）+ bonded_report.json`；
  - `derive-torsions` / `assign-charges` 输出 `coords.json`（记录实际
    使用坐标的来源 inpcrd:/coords:）；`convert` / `check-equivalence`
    可直接作用于最终 IR + 坐标（tleap 腿参数仅在三件套给齐时传）；
  - 中途失败不写 ir/coords/report（产物状态约定）。
- `BondedReport.to_dict()`（bonded/assign.py）：报告 JSON 序列化。
- 测试：`test_ir_coords.py`（7）+ `test_pipeline_chain.py` 的合同用例。

## 3. C3：H_fixed 依赖 provenance

- `bonded/staleness.py`：`h_fixed_snapshot(ir)`（电荷+电荷来源+proper/
  improper+LJ+混合规则+1-4 缩放+exclusions/pairs14，逐项数值，不用 hash）。
- `bonded/global_fit.py`：拟合把快照连同 Hessian 身份（provider/method/
  sha256/positions_sha256）与 `fd_step_nm` 写进
  `provenance["bonded"]["h_fixed_dependency"]`（inputs 取**拟合输入 IR**）。

## 4. C4：依赖过期规则 + global_fit CLI

- `staleness(ir)`：global_fit 的 IR 比对记录快照 vs 当前快照，报过期字段；
  依赖未记录（旧 IR）也算过期。MLP 二面角一侧：`mm_subtraction_snapshot`
  （bond/angle/电荷/LJ）记进 `provenance["torsion"]["mm_dependency"]`
  （`torsion/mlp.py::_finish`；refit 路径自动刷新）。
- CLI `global-fit`：标准导出路径；下界命中（n_at_bound>0，决策 0014）→
  report 记 `status=failed_lower_bound`、IR/coords 落盘供诊断、**不导出
  三后端**、exit 1；成功 → status=ok + 三后端。不静默截正、不换回 GAFF2。
- CLI `check-staleness --ir-in F`：exit 0 无过期 / 1 有过期（组合完成前的
  拒绝静默沿用）。
- 测试：`test_staleness.py`（7）+ `test_pipeline_chain.py` 的两个组合用例
  （真路径 ok + monkeypatch 契约用例；数学边界复用既有
  `test_global_fit.py::test_lower_bound_hit`）。

## 5. C2：局部路线串联（靠 C1 原语实现，不建调度器）

链条 = 既有命令 + C1 合同：基础 AMBER IR → `assign-charges` →
`derive-bonded --ir-in --coords-in` → `convert --ir-in --coords-in 三后端`
（Hessian 优化几何随 coords sidecar 传递，不回退 inpcrd）。
`test_pipeline_chain.py::test_local_route_chain_preserves_params_and_coords`
用合成谐振 fchk（H = Σ k_t b_tᵀb_t，§5.1，PSD 无虚频）驱动，不需 xtb。
默认局部方法仍是 mseminario（设计 §6.1；CLI 不改默认值）。

## 6. C5：组合验收（单元/组合级）

`test_pipeline_chain.py` 覆盖：完整参数保留（非默认电荷+Hessian 键角到
三后端不变）、坐标来源正确（hessian:gaussian-fchk）、换序 coords/错拓扑
sdf 拒绝、global-fit 后换电荷 → check-staleness 过期、下界命中产物状态。
正式外部组合验收（跑 Gate、sander/gmx 腿）按 10-15 计划 §5 C5 另报资源，
归接收方。

## 7. MVP-5：fallback / 不确定度调度框架（设计 §8.1/§16/§23）

- 新增 `src/ffcompiler/hessian/fallback.py`：`Tier`（可注入 runner）+
  `run_with_fallback`（顺序尝试，`check_hessian` 质量门；错误/拒绝/不可用
  分状态记录；接受即返回）+ `default_ladder`（gfnff → gfn2 → dft 文件级；
  DFT 无文件 → `TierUnavailableError`，本机无 QM 程序如实报告）。
- 第一版"可疑"判定 = check_hessian 报错条件 + runner 抛错；警告只记录
  不升级（阈值收紧属科学决策，不猜）。CLI 不暴露 auto（避免暗地连跑
  xtb 两级）；调度研究后置。
- 测试 `test_fallback.py`（7，假 runner，不调 xtb）。

## 8. MVP-4：learned BCC 骨架 + 立项（设计 §13–15）

- 骨架 `src/ffcompiler/charge/learned_bcc.py`：`BondCorrectionModel`
  协议（每键 Δq，符号约定 IR 键 (i,j) = i 端 +Δq / j 端 −Δq）、
  `LearnedBCCCharges(ChargeProvider)`（q = q⁽⁰⁾ + Σ ±Δq，工件缺失 →
  ChargeError 不降级）、`IdentityCorrections` 占位（Δq≡0，仅联调）。
- 立项草案 `docs/superpowers/plans/2026-10-17-mvp4-learned-bcc-proposal.md`：
  S1–S6 科学决策全部留白待批；批准前禁止训练/数据生成/正式参数。
- 测试 `test_learned_bcc.py`（6，守恒/符号/错误路径/provenance）。

## 9. G1.8 文档收口（顺带完成）

`docs/theory/hessian_methods.md` 四处修正（§4.3 对角元补平方、§6 首段
rad⁻²、§5.3 q⁰ 措辞 + 凸性说明、§2 方向记号与符号不变性），头部状态改
"已审阅（2026-10-13 A3 §2 实质通过；2026-10-17 收口）"。

## 10. 明确没做（不在本次范围）

- B1–B4 诊断（需要算力与 10-13 缓存的只读分析）；
- D 独立科学评价设计；MVP-4 的任何科学内容；
- RESP provider（决策 0018 往后放）；frcmod 路线（YAGNI）；
- 正式 Gate 运行与版本号改动（等用户）。
