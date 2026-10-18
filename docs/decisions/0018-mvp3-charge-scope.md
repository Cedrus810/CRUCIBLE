# 决策 0018：MVP-3（电荷 provider 框架）的范围

- 2026-10-10 用户决定（对应 roadmap §1 MVP-3、§3 第 3 条、设计 §11/§23）：
  1. **RESP 往后放**，MVP-3 不做。用户原话："暂时不做 RESP，因为 RESP 有说法，你给他往后放"。
     设计 §23 MVP-3 列出的 RESP 推迟；接口给 RESP 留位置，不实现。
  2. **provider 三种**：载入电荷（保留现有或从文件读）、AM1-BCC（AmberTools antechamber/sqm）、
     NAGL（openff-nagl 0.6.1 + openff-gnn-am1bcc-1.0.0，已在 openmm_dev，不装包）。
  3. **电荷与二面角的耦合**：管线顺序固定为"先电荷、后二面角"。给已有 MLP 二面角的 IR 换电荷时，
     用缓存的扫描重新扣 MM、重新拟合（扫描几何与参考能量和电荷无关，不重扫）；provenance 记录
     二面角是在哪套电荷下拟合的。
