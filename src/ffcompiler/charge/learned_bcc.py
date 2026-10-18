"""learned BCC 电荷骨架（MVP-4，设计 §13–15；2026-10-17 只搭管线，科学部分未立项）。

数学形式（设计 §13，锁定）：

    q_i = q_i^(0) + Σ_j Δq_ji，且对每条 bond Δq_ij = −Δq_ji

即**基础电荷 + 逐键电荷转移修正**，总电荷按构造守恒（Σq = Q_formal 的
残差仍走 apply.equalize_residual 的均摊）。修正不由离散 bond type 决定，
而由感知完整局部化学环境的模型给出（generalized / learned BCC）。

**当前状态（2026-10-17）：只有管线骨架，没有任何科学内容。**

- :class:`BondCorrectionModel` 是模型接口（输入 IR/坐标/分子，输出每条键
  的 Δq，符号约定：IR 键 (i, j) 书写顺序为 i 端 +Δq、j 端 −Δq）；
- :class:`LearnedBCCCharges` 是 ChargeProvider 实现：注入 base 电荷来源与
  模型工件；工件缺失/未训练一律 :class:`ChargeError`，**不静默降级**成
  任何别的电荷。base 电荷 provider **必须显式给出**（S1 未定，不设默认）：
  2026-10-17 验收发现"取 IR 当前电荷"做 base 时连续赋两次修正会叠加；
- :class:`IdentityCorrections` 是 Δq≡0 的占位模型，仅用于管线联调与测试
  （输出=基础电荷），不代表任何已训练模型。

训练目标（设计 §14 的多目标损失）、数据集、网络结构、base charge 选择
都是**立项后**的研究决策，见
``docs/superpowers/plans/2026-10-17-mvp4-learned-bcc-proposal.md``；在那
之前本模块不得用于生成任何正式参数。
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol

import numpy as np

from ffcompiler.charge.base import ChargeError, ChargeProvider


class BondCorrectionModel(Protocol):
    """每键电荷转移修正模型接口（MVP-4 立项后实现）。"""

    def predict(self, ir, positions_nm: np.ndarray, mol) -> np.ndarray:
        """返回 (n_bonds,) 的 Δq（单位 e）；第 t 条对应 ir.bonds[t]。

        符号约定：对 ir.bonds[t] = (i, j)，+Δq 加在 i、−Δq 加在 j。
        模型训练时的约定必须与此一致。
        """

    def identity(self) -> dict:
        """进 provenance 的身份信息（结构、版本、训练摘要）。"""


class IdentityCorrections:
    """Δq≡0 的占位模型：输出等于基础电荷。仅用于管线联调/测试。"""

    def predict(self, ir, positions_nm: np.ndarray, mol) -> np.ndarray:
        return np.zeros(len(ir.bonds), dtype=np.float64)

    def identity(self) -> dict:
        return {"name": "identity-untrained", "trained": False}


class LearnedBCCCharges(ChargeProvider):
    """q = base + Σ_bond ±Δq 的电荷 provider（设计 §13）。

    base_provider 给 q⁽⁰⁾（必填；它若取 IR 当前电荷，如 ``LoadedCharges()``，
    对已含修正的 IR 再赋一次会叠加修正，调用方负责）；model 给逐键修正。
    model_path 指向训练工件；现在没有已训练工件，传入任何路径都会在
    compute 时检查存在性——文件不存在 → ChargeError（明确报"模型未训练/
    工件缺失"，不 fallback）。
    """

    method = "learned_bcc"

    def __init__(self, model: BondCorrectionModel, *, base_provider: ChargeProvider,
                 model_path: Path | None = None):
        if model_path is not None and not Path(model_path).exists():
            raise ChargeError(
                f"learned_bcc 模型工件不存在：{model_path}。MVP-4 尚未立项训练，"
                "当前只有骨架（charge/learned_bcc.py 的 IdentityCorrections 仅作联调）")
        self._model = model
        self._base = base_provider
        self._model_path = model_path

    def compute(self, ir, positions_nm: np.ndarray, mol) -> np.ndarray:
        q0 = np.asarray(self._base.compute(ir, positions_nm, mol), dtype=np.float64)
        n = len(ir.atoms)
        if q0.shape != (n,):
            raise ChargeError(f"learned_bcc: 基础电荷形状 {q0.shape} != ({n},)")
        delta = np.asarray(self._model.predict(ir, positions_nm, mol), dtype=np.float64)
        if delta.shape != (len(ir.bonds),):
            raise ChargeError(
                f"learned_bcc: 模型输出形状 {delta.shape} != ({len(ir.bonds)},)（每键一个 Δq）")
        if not np.all(np.isfinite(delta)):
            raise ChargeError("learned_bcc: 模型输出含非有限值（NaN/Inf）")
        q = q0.copy()
        for bond, d in zip(ir.bonds, delta, strict=True):
            q[bond.i] += float(d)
            q[bond.j] -= float(d)
        return q

    def source_info(self) -> dict:
        info = {
            "model": self._model.identity(),
            "base": {
                "method": getattr(self._base, "method", ""),
                "source": self._base.source_info(),
            },
            "formula": "q_i = q_i^(0) + sum_j dq_ji (bond-antisymmetric)",
        }
        if self._model_path is not None:
            info["model_path"] = str(self._model_path)
        return info
