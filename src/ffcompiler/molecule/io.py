"""SDF 分子入口的读入与检查（设计 §3.1：只检查、不优化、不生成几何）。

入口合同（§3.1 + 2026-11-11 计划 §2）：
- 必须能被 RDKit 解析（removeHs=False）；
- 必须自带 3D 构象——缺 3D 直接报错，不做 ETKDG 兜底（构象生成是显式步骤，§10.4/§6.4.5）；
- formal charge 总和 = 声明净电荷（未声明时取 SDF 总和并回显给 antechamber -nc）；
- 无严重原子穿距：任意两原子距离 < :data:`CLASH_MIN_A` 报错；
- 记录 SDF sha256 供 provenance。
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np

CLASH_MIN_A = 0.5  # 任意两原子最小距离（Å）：§3.1"无严重原子穿距"的量化门


class MoleculeInputError(ValueError):
    """SDF 不满足入口几何合同（§3.1）。"""


@dataclass(frozen=True)
class MoleculeInput:
    """入口检查通过后的最小身份信息（provenance 与 -nc 之用）。"""

    path: Path
    name: str
    net_charge: int
    n_atoms: int
    sha256: str


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_sdf(path: Path, *, name: str | None = None,
             net_charge: int | None = None) -> MoleculeInput:
    """读入并检查 SDF；不修改、不写出任何坐标。"""
    from rdkit import Chem

    path = Path(path)
    if not path.is_file():
        raise MoleculeInputError(f"SDF 不存在：{path}")
    mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=True)
    if mol is None:
        raise MoleculeInputError(f"SDF 解析失败（RDKit 无法读取）：{path}")
    if mol.GetNumConformers() == 0:
        raise MoleculeInputError(
            f"{path}：没有构象。§3.1 入口不生成构象（不允许静默 ETKDG），"
            "请提供带 3D 坐标的 SDF；构象生成是显式步骤（设计 §10.4、§6.4.5）")
    conf = mol.GetConformer()
    pos = np.asarray(conf.GetPositions(), dtype=float)
    if not conf.Is3D() or float(np.abs(pos).sum()) == 0.0:
        raise MoleculeInputError(
            f"{path}：构象不是 3D（Is3D={conf.Is3D()} 或坐标全零）；入口拒绝（§3.1）")
    n = mol.GetNumAtoms()
    dist = np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(-1))
    iu = np.triu_indices(n, k=1)
    pair_d = dist[iu]
    k = int(np.argmin(pair_d))
    if float(pair_d[k]) < CLASH_MIN_A:
        i, j = int(iu[0][k]), int(iu[1][k])
        raise MoleculeInputError(
            f"{path}：原子 {i} 与 {j} 距离 {pair_d[k]:.3f} Å < {CLASH_MIN_A} Å"
            "（严重原子穿距，§3.1 拒绝）")
    formal = [a.GetFormalCharge() for a in mol.GetAtoms()]
    sdf_q = int(sum(formal))
    if net_charge is None:
        net_charge = sdf_q
    elif net_charge != sdf_q:
        raise MoleculeInputError(
            f"{path}：SDF formal charge 总和 {sdf_q} 与声明净电荷 {net_charge} 不一致（§3.1）")
    if name:
        mol_name = name
    elif mol.HasProp("_Name") and mol.GetProp("_Name"):
        mol_name = mol.GetProp("_Name")
    else:
        mol_name = path.stem
    return MoleculeInput(path=path, name=mol_name, net_charge=net_charge,
                         n_atoms=n, sha256=_sha256(path))
