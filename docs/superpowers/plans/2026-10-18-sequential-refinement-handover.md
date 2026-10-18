# 可选顺序重拟合：实现与计算交接

日期：2026-10-18，Asia/Tokyo。

**状态：实现、测试代码与交接已准备；只做文本/AST 静态核对，所有新增/修改测试均未运行。** 用户授权“自己干除了计算之外的事儿”，随后要求继续准备结案。本轮未运行 pytest（含 collect-only）、CLI、扫描、拟合、benchmark 或 Gate；未加载项目/模型来做检查，未安装环境、运行 Git 或修改历史科学产物。

## 0. 接收方验证（2026-10-18，代理，node-A）

按设计稿与 B3 规格核对实现（默认单遍；provenance 记模式/轮数/顺序；换电荷按记录重放；原兜底不升级；replaced_terms 不改写），未发现需要修改的代码问题。验证记录 `experiments/20261018_refinement_verify/`：

- §5 相关 6 个测试文件：66 passed；非 external 全量：524 passed。
- `refinement_rounds=3`（10-13 G2 缓存，cache_only，参考势调用即报错）与 B3 机制实验第 3 轮：项与逐键指标逐位相同。
- 默认 0 轮与改动前单遍（10-18 G2 存档 IR，同机）：6/6 分子逐位相同。
- G0–G3 在新源码上重跑（`jobs/gate_mvp*_20261018_2/`，08:50–09:15）：全部 passed，每个 Gate 的 pytest 524 + 59，数值与之前各轮一致。Gate 走默认单遍。下面第 5 节的"未运行"指施工时的状态。

## 1. 接口与范围

- Python：`MLPTorsionProvider(..., refinement_rounds=0)`。
- CLI：`derive-torsions --refinement-rounds {0,1,2,3}`；默认 0。
- 这里的数值是**原有单遍之后的额外轮数**。3 对应 B3 机制实验的“第 0 轮 + 3 轮”，不是把单遍算进三轮。
- 可选能力不等于推荐默认参数：单遍仍为默认，版本号、科学门槛和兜底规则不变。
- 顺序为中心键编号升序；每次仍只用本键扫描、扣除当前其余参数，立即替换本键项。有限轮之后不保证联合目标单调下降或等于联合最小二乘解。
- 更新只针对本次单遍成功的键；缓存、模型、Fourier 阶数和扫描几何不因启用额外轮数而改变。没有新增扫描质量阈值，也没有对 nma 特设自动兜底。

## 2. 初次拟合与换电荷共用实现

初次拟合：原有单遍 → 可选顺序更新 → 记录最终 IR 与报告。

换电荷：从 `replaced_terms` 恢复最初项 → 原有单遍 → 按 IR 记录的轮数顺序更新。复用同一个 provider 路径，参考势保持扫描时身份，严格只读扫描缓存；没有添加 assign-charges 的独立覆盖选项，避免两条路径无意使用不同协议。

`fit_protocol` 记录单遍/顺序方式、额外轮数与键顺序。旧 IR 没有这个字段时继续按旧单遍流程处理；未知协议、非法预算或键顺序与原 mlp 集合不一致时报错。没有 MLP 项时无需执行重拟合。

多轮过程中仅替换当前参数，**不覆盖最初 `replaced_terms`**；原兜底行和原参数在换电荷后仍保留，缓存补齐不自动升级（决策 0020）。

## 3. 报告与依赖

启用顺序更新时，torsion_report.json 和 IR provenance 包含：

- `refinement.extra_rounds` / `completed_rounds` / `order`；
- `objective=per_bond_scan`、`evaluation=in_sample`；
- `history`：初始单遍及每轮结束时的逐键拟合残差、最佳平移/原有最小值平移剖面 RMSE、滞后、几何变化；
- 每个分子的逐帧等权汇总残差。每条扫描各扣自身常数；跨分子总值应按帧数重新汇总，不能简单平均各分子 RMS；
- `fixed_proper_dependency`：不参与顺序更新的 proper，包括原兜底键。之后改动这些项，check-staleness 会报告过期。

单遍的既有 to_dict 报告形状保留；两种模式都会记录 fit_protocol。报告对象与已生成 IR 的嵌套来源记录分开复制，修改报告不会反向改写 IR。

没有成功键时不做空轮数求解，completed_rounds=0、汇总残差为 null。非有限剖面/拟合错误不返回部分成功的 refined IR；CLI 沿现有错误边界返回失败。建议所有运行使用新目录，防止把旧标准文件误当成本次结果。

## 4. 已准备的代码

| 文件 | 内容 |
|---|---|
| src/ffcompiler/torsion/mlp.py | 有界顺序更新、逐轮报告、协议记录/恢复、原始项保护、来源记录分离 |
| src/ffcompiler/torsion/mm.py | 完整 MM 求值函数；原扣除函数委托同一批量求值实现 |
| src/ffcompiler/bonded/staleness.py | 固定 proper 快照及非法协议/依赖更新识别 |
| src/ffcompiler/cli.py | derive-torsions 可选参数，MM/MACE 两条参考势路径同样传递 |
| tests/unit/test_torsion_refinement.py | 12 个测试函数，包含参数化非法预算用例；使用解析耦合曲线，不调用参考势 |
| tests/unit/test_pipeline_chain.py | 原完整串联测试扩展为额外轮数 0/3 两种模式 |
| tests/unit/test_cli.py | 非法 CLI 预算在加载输入/计算前拒绝 |

## 5. 接收方验证顺序（本轮均未执行）

以下命令是运行交接，不代表授权在任意主机执行。选择用户允许的计算节点/机器，先报资源与输出；登录节点不跑任何测试体。所有结果写新目录。使用既有 openmm_dev 与 PYTHONPATH=src，不安装包。

1. 先跑新增解析耦合回归与 CLI 参数回归，再跑 existing refit、完整链路和 staleness 文件。新解析用例包含可手算的递推：a ← 1.5 − b/2；b ← 1.5 − a/2。它用于确认更新使用当前其他键的项、顺序确定、历史保留和换电荷走同一协议。
2. 重点看：0 轮与旧数值行为一致；同电荷重复赋值幂等；3 轮 JSON 往返后仍沿用 3 轮；缺缓存明确失败；原兜底不升级；报告修改不污染 IR；固定 proper 修改后过期；非法协议不静默改回单遍。
3. 相关文件通过后，再决定全量单元/外部回归及新源码的正式 Gate。10-18 早间通过记录不自动覆盖本次实现。

```bash
# 仅供接收方在获准的计算资源上执行；此处未运行。
source scripts/env.sh
PYTHONPATH=src python -m pytest \
  tests/unit/test_torsion_refinement.py \
  tests/unit/test_torsion_refit.py \
  tests/unit/test_pipeline_chain.py \
  tests/unit/test_cli.py \
  tests/unit/test_staleness.py \
  tests/unit/test_torsion_scan.py \
  -m "not external" -q
```

接收方需要独立保存原始输出和身份记录；不能把本文件的“准备好”记成测试通过。新断言失败先保留反例与误差，不以放宽断言凑绿。

## 6. B4 与 C5 的运行准备

### B3 接口一致性复核

接收方可从 10-13 G2 缓存只读构造 `cache_only=True` 的 provider，比较 refinement_rounds=0/3 与已存 B3 实验结果。只比较相同源 IR、模型身份、扫描键和指标；不同 BLAS 的末位差异按既有规则记录，不能要求超出 CSV 精度的逐字一致。

该比较需运行数值求值，因此本轮未做。不要调用默认 derive-torsions 来假装只读缓存：CLI 没有 cache-only 选项，缓存不命中时会扫描。需要严格禁止新扫描的执行者应使用 provider 的 cache_only=True 并传入会报错的参考势工厂。

### B4 留出评价

执行前固定留出数据来源、分子/构象拆分、指标和预算。训练用的同一组扫描帧不能同时算独立测试。先比较默认单遍、顺序三轮和 GAFF2；保留能垒、极小值位置、能量差与力指标，逐键报告退化/失败。

nma 1-3 继续单列扫描质量问题；是否多起点/重扫/换分支处理单独决定，不从这批极差分布自动定通用阈值。本轮不替代用户的科学决策。

### C5 真实组合验收

执行者以新运行目录记录：Hessian bond/angle → 最终电荷 → 真实 MACE 二面角 → 最终 IR/坐标 → 三后端分项能量/力；并检查换电荷后的缓存重拟合、兜底集合、依赖刷新。先覆盖局部路线，再单独处理 global-fit 双向依赖。

0/3 两种模式分开标记；缓存命中与重扫分开记录。已定容差、子集 scope 规则与决策 0017 的 sander 报告范围不变；check-equivalence CLI 的原始退出码要保留，不能静默把 sander 超差改为通过。新的组合通过条件在运行前锁定。

global-fit 与 torsion 的交替策略、科学质量门槛、MVP-4 训练和 MVP-5 调度不在本次实现中自动决定。

## 7. 结案范围

本次完成的是无计算施工与交接。剩余工作是运行新回归、数值一致性验证、B4/C5 和原有科学决策。没有新增“通过”记录；旧 Gate、缓存、诊断与顺序实验的历史结果原样保留。
