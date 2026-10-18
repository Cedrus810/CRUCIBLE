"""Hessian provider 的三级 fallback / 不确定度调度框架（MVP-5，设计 §8.1/§16）。

设计锁定的升级链（§16）：

    GFN-FF →（可疑）→ GFN2-xTB →（可疑/特殊化学）→ DFT Hessian

实现约束（设计 §23 MVP-5 原话"不要一开始就搞复杂调度"）：

- 每一级是一个可注入的 runner（``Callable[..., CartesianHessian]``），
  框架本身只做**顺序尝试 + 质量门**，不做预算估计、不做模型选择；
- "可疑"的判定第一版固定为 :func:`~ffcompiler.hessian.base.check_hessian`
  的报错条件（不对称 / 平移破坏 / 非驻点 / 虚频），外加 runner 自身抛错。
  警告（如 −20…0 cm⁻¹ 的浅虚频）**只记录不升级**——阈值收紧属科学决策，
  留给立项后的调度研究，不在本框架里猜；
- DFT 级在本机没有 QM 程序（roadmap §0）：默认 runner 是"必须显式给
  fchk/.hess 文件"的读取器；没给文件 → 该级不可用，如实报告，不静默跳过。

模块级函数 :func:`default_ladder` 组装 xtb 两级 + 文件级 DFT；单测用假
runner 注入，不调用 xtb。CLI 不直接暴露 auto 模式（避免暗地里连跑两级
xtb 的资源消耗）；要自动升级的调用方自己组装 ladder。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Mapping

from ffcompiler.hessian.base import CartesianHessian, HessianQualityError, check_hessian

HessianRunner = Callable[..., CartesianHessian]


@dataclass(frozen=True)
class TierResult:
    tier: str
    status: str                       # "accepted" | "rejected" | "unavailable" | "error"
    detail: str = ""
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class FallbackReport:
    accepted_tier: str | None
    attempts: tuple[TierResult, ...] = ()

    def to_dict(self) -> dict:
        return {
            "accepted_tier": self.accepted_tier,
            "attempts": [
                {"tier": a.tier, "status": a.status, "detail": a.detail,
                 "warnings": list(a.warnings)}
                for a in self.attempts
            ],
        }


@dataclass(frozen=True)
class Tier:
    name: str
    runner: HessianRunner
    kwargs: Mapping = field(default_factory=dict)

    def __str__(self) -> str:  # 报告/报错信息用
        return self.name


def run_with_fallback(
    tiers: list[Tier],
    *,
    masses=None,
    allow_nonstationary: bool = False,
) -> tuple[CartesianHessian, FallbackReport]:
    """按顺序尝试各级 runner，用 check_hessian 做质量门。

    - runner 抛错（外部程序失败等）→ 记 error，升级下一级；
    - check_hessian 抛 HessianQualityError → 记 rejected，升级下一级；
    - 通过 → 记 accepted 并**立即返回**（不白跑更高等级）；
    - 某级 runner 声明不可用（抛 :class:`TierUnavailableError`）→ 记
      unavailable，升级下一级；
    - 全部失败 → FallbackError（带完整 attempts 报告）。
    """
    attempts: list[TierResult] = []
    for tier in tiers:
        try:
            h = tier.runner(**tier.kwargs)
        except TierUnavailableError as exc:
            attempts.append(TierResult(tier.name, "unavailable", str(exc)))
            continue
        except Exception as exc:  # noqa: BLE001 - 外部程序失败的种类不可枚举
            attempts.append(TierResult(tier.name, "error", f"{type(exc).__name__}: {exc}"))
            continue
        try:
            diagnostics = check_hessian(
                h, masses=masses, allow_nonstationary=allow_nonstationary
            )
        except HessianQualityError as exc:
            attempts.append(TierResult(tier.name, "rejected", str(exc)))
            continue
        attempts.append(
            TierResult(tier.name, "accepted", warnings=tuple(diagnostics.warnings))
        )
        return h, FallbackReport(accepted_tier=tier.name, attempts=tuple(attempts))
    raise FallbackError(
        "全部 Hessian 级别都未通过质量门：\n"
        + "\n".join(
            f"  {a.tier}: {a.status} {a.detail}".rstrip() for a in attempts
        ),
        FallbackReport(accepted_tier=None, attempts=tuple(attempts)),
    )


class FallbackError(RuntimeError):
    """所有级别失败。``.report`` 带逐级尝试记录。"""

    def __init__(self, message: str, report: FallbackReport):
        super().__init__(message)
        self.report = report


class TierUnavailableError(RuntimeError):
    """该级别在当前配置下不可用（如 DFT 级没有给文件）。"""


def default_ladder(
    atomic_numbers,
    positions_nm,
    *,
    charge: int,
    workdir: Path,
    dft_file: Path | None = None,
    xtb_binary: Path | None = None,
) -> list[Tier]:
    """默认三级：gfnff → gfn2 → dft（文件级）。

    xtb 两级共享 run_xtb_hessian 的接口；DFT 级必须显式给 fchk 或 .hess
    文件（按后缀选读取器），否则该级 unavailable。
    """
    from ffcompiler.hessian.xtb import XTB_BINARY, run_xtb_hessian

    binary = xtb_binary or XTB_BINARY
    tiers = [
        Tier("gfnff", run_xtb_hessian,
             dict(atomic_numbers=atomic_numbers, positions_nm=positions_nm,
                  method="gfnff", charge=charge, workdir=workdir / "gfnff",
                  binary=binary)),
        Tier("gfn2", run_xtb_hessian,
             dict(atomic_numbers=atomic_numbers, positions_nm=positions_nm,
                  method="gfn2", charge=charge, workdir=workdir / "gfn2",
                  binary=binary)),
    ]

    expected_z = [int(z) for z in atomic_numbers]

    def _dft_runner(path: Path | None = dft_file) -> CartesianHessian:
        if path is None:
            raise TierUnavailableError(
                "DFT 级需要 --hessian-file（fchk 或 ORCA .hess）；本机没有 QM 程序，"
                "不能自动生成 DFT Hessian（roadmap §0）")
        if path.suffix == ".hess":
            from ffcompiler.hessian.orca import read_orca_hess
            h = read_orca_hess(path)
        else:
            from ffcompiler.hessian.gaussian import read_fchk
            h = read_fchk(path)
        # 文件是外部给的：原子序数序列必须与本分子一致（xtb 两级由本分子直接算出）
        if [int(z) for z in h.atomic_numbers] != expected_z:
            raise ValueError(f"DFT Hessian 文件 {path} 的原子序数序列与本分子不一致")
        return h

    tiers.append(Tier("dft", _dft_runner))
    return tiers
