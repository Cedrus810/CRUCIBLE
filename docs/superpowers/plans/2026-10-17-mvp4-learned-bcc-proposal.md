# MVP-4（learned BCC）研究立项草案

日期：2026-10-17。性质：**立项材料，未获用户批准；一切科学决策空白待补**。
用户 2026-10-17 指示"没搭建的你搭建"，本文件与 `src/ffcompiler/charge/learned_bcc.py`
的骨架是该指示下允许的准备动作（10-15 计划 §6："MVP-4 可以准备研究立项"）；
**训练、数据生成、任何正式参数产出都在批准之前禁止**。

## 1. 定位（设计稿 §13–15 已锁定）

- 不是"ML 直接预测 RESP"（那只是 surrogate/加速），而是
  **environment-aware、charge-conserving 的逐键电荷转移修正模型**
  （generalized / learned BCC）：
  `q_i = q_i^(0) + Σ_j Δq_ji`，`Δq_ij = −Δq_ji`，总电荷按构造守恒。
- 修正由连续感知完整局部化学环境的小模型给出，不再只由离散 bond type 决定。
- Phase 1（设计 §15）：新电荷 + **冻结 GAFF2 LJ**，先测
  ESP/偶极/HFE/溶剂转移/相互作用能/RBFE；Phase 2 才考虑 charge/LJ 联合。
  第一版不把 charge、LJ、torsion 一起训练。

## 2. 已有的基础设施（可直接复用）

- ChargeProvider 框架与 `apply.assign_charges`（均摊残差、provenance、
  有 MLP 二面角时用缓存扫描重拟合——换电荷的下游闭环已验证，G3）。
- 三套可比 baseline：载入电荷 / AM1-BCC（与固件逐位相同）/ NAGL。
- 符号约定已在骨架锁定：IR 键 (i, j) 书写顺序为 i 端 +Δq、j 端 −Δq；
  训练数据构造必须用同一约定。
- 交叉验证入口：D（独立科学准确度评价，10-15 计划 §6：不能用
  NAGL/AM1-BCC 相互一致替代物理准确度）。

## 3. 立项后必须先定的科学决策（全部空白，逐条要用户/评审拍板）

| # | 决策 | 备选与备注 |
|---|---|---|
| S1 | base charge q⁽⁰⁾ | EEQ / 半经验 / 廉价 QM 划分 / 载入值；影响可学性与可迁移性 |
| S2 | 模型结构 | 小图网络 / 等变网络；输入特征（元素、连接、芳香性、环、形式电荷） |
| S3 | 训练数据 | 设计 §14 的分层：气相 ESP/偶极、水相 HFE/相互作用、蛋白样片段、有机溶剂转移；规模与来源待定 |
| S4 | 损失权重 | 多目标 `w_ESP/w_μ/w_HFE/w_solv/w_int/w_reg`；先单目标消融还是直接多目标 |
| S5 | 留出与评价协议 | 分子拆分规则、参照方法（RESP2/ABCG2 对比）、D 的独立参照 |
| S6 | 不确定度输出 | 是否输出 per-bond 不确定度（与 MVP-5 的调度共享语义） |

## 4. 工程接口（已随骨架冻结，立项不改变）

- `charge/learned_bcc.py`：`BondCorrectionModel` 协议（`predict` 每键 Δq +
  `identity`）、`LearnedBCCCharges(ChargeProvider)`、占位
  `IdentityCorrections`（Δq≡0，仅联调用）。
- 工件缺失/未训练 → `ChargeError`，不静默降级；provenance 记
  `model`/`base`/公式与工件路径。
- base 电荷 provider 必须显式给出（2026-10-17 验收后改；原"默认取 IR 当前电荷"
  在重复赋电荷时会叠加修正）。S1 定下之前不设默认。
- CLI 接入在**有第一个真实工件之后**再做（`assign-charges --method
  learned_bcc --model ...`），避免暴露空壳选项。

## 5. 建议的立项后顺序（供讨论，不是决定）

1. S1+S2 的小规模原型（几百分子），单目标 ESP 起步；
2. S3/S4 数据与权重设计评审；
3. D 评价协议先行冻结，再开始训练；
4. Phase 1 冻结 LJ 的全量 benchmark（设计 §15 清单）；
5. 达标后开放 CLI 与 preset（设计 §17 的 `charge: new_charge`）。
