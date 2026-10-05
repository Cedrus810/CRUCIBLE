"""
配体二面角 MLP 扫描（取自 ABFE_IBS）

来源：ABFE_IBS 仓库 commit 56c2a92 的 `abfe_core.py`
  - `OrbScanner`（L2772）：真空 MLP 势下对单个二面角做刚性 1D 扫描，-180..170°，步长 10°
  - `OrbMMHybridFactory.get_rotatable_torsions_rdkit`（L3013）：枚举配体全部可转键的四原子二面角
  - `_build_openmmml_kwargs` / `_build_mace_potential`（L52 / L82）：openmm-ml 构建辅助
代码原样搬过来，只把 device 默认值从仓库的全局设备探测换成 "cuda"。

说明：
  - 类名叫 Orb，但默认模型是 `mace-off24-medium`（走 openmm-ml 的 MLPotential），
    换成别的 openmm-ml 支持的模型名即可。
  - 扫描是刚性扫描（RDKit `SetDihedralDeg` 直接转，不做约束优化）。
  - 返回的是真空 MLP 能量 E(φ)（kJ/mol），没有扣 MM 项，也没有拟合。
  - ABFE_IBS 那边是把修正项写成 Fourier 级数 offset + Σ c_n cos(nφ) + s_n sin(nφ)，
    用 CustomTorsionForce 加回 System，不写回拓扑。
    写回拓扑（[ dihedrals ] 的 k/n/phase）需要你那边自己做。
"""

import os
from typing import Any, Dict, Optional

import numpy as np
import openmm
import openmm.app
from openmm import unit

try:
    from openmmml import MLPotential
    HAS_ORB = True
except ImportError:
    MLPotential = None
    HAS_ORB = False


def _build_openmmml_kwargs(
    device: Optional[str] = None,
    precision: Optional[str] = None,
    return_energy_type: Optional[str] = None,
    charge: Optional[int] = None,
    multiplicity: Optional[int] = None,
) -> Dict[str, Any]:
    """
    统一按 openmm-ml 官方接口组织 MLPotential.createSystem() 参数。
    是否显式传 precision 由调用侧决定；若不传则遵循 openmm-ml 的模型默认精度。
    """
    kwargs: Dict[str, Any] = {}
    if return_energy_type is not None:
        kwargs["returnEnergyType"] = return_energy_type
    if device is not None:
        kwargs["device"] = device
    if precision in ("single", "double"):
        kwargs["precision"] = precision
    if charge is not None:
        kwargs["charge"] = charge
    if multiplicity is not None:
        kwargs["multiplicity"] = multiplicity
    return kwargs


_MACE_LOCAL_MODEL_PATHS = {
    "mace-off24-medium": os.path.expanduser("~/.cache/mace/MACE-OFF24_medium.model"),
}


def _build_mace_potential(model_name: str):
    """
    部分 openmm-ml 版本的模型注册表里没有 mace-off24-medium 等较新的预训练名，
    但支持用 model_name="mace" + modelPath=本地缓存权重文件的方式加载同一个模型。
    这里对已知的本地缓存做名字 -> 路径映射；找不到对应缓存文件时退回原始按名加载
    （交给 openmm-ml 自己报错或解析）。
    """
    cached_path = _MACE_LOCAL_MODEL_PATHS.get(model_name)
    if cached_path and os.path.isfile(cached_path):
        return MLPotential("mace", modelPath=cached_path)
    return MLPotential(model_name)


class OrbScanner:
    def __init__(
        self,
        model_name="mace-off24-medium",
        n_order=6,
        charge=0,
        multiplicity=1,
        device="cuda",
    ):
        self.n_order = n_order
        self.charge = charge
        self.multiplicity = multiplicity
        self.device = device
        self.model_name = model_name
        self.context = None
        self.system = None
        if HAS_ORB:
            self.potential = _build_mace_potential(model_name)

    def _setup_vacuum_context(self, rdkit_mol):
        if self.context is None and HAS_ORB:
            vacuum_top = openmm.app.Topology()
            c = vacuum_top.addChain()
            res = vacuum_top.addResidue("MOL", c)
            for atom in rdkit_mol.GetAtoms():
                vacuum_top.addAtom(
                    atom.GetSymbol(),
                    openmm.app.element.Element.getByAtomicNumber(atom.GetAtomicNum()),
                    res,
                )
            self.system = self.potential.createSystem(vacuum_top, **_build_openmmml_kwargs(
                device=self.device,
                return_energy_type="interaction_energy",
                charge=self.charge,
                multiplicity=self.multiplicity,
            ))
            self.context = openmm.Context(self.system, openmm.VerletIntegrator(0.001))

    def run_torsion_scan(self, rdkit_mol, torsion_indices):
        self._setup_vacuum_context(rdkit_mol)
        angles = np.arange(-180, 180, 10)
        e_list = []
        conf = rdkit_mol.GetConformer()
        for ang in angles:
            import rdkit.Chem.rdMolTransforms as rdt

            rdt.SetDihedralDeg(conf, *torsion_indices, float(ang))
            pos_nm = (
                np.array(
                    [conf.GetAtomPosition(i) for i in range(rdkit_mol.GetNumAtoms())]
                )
                * 0.1
            )
            self.context.setPositions([openmm.Vec3(*p) for p in pos_nm])
            e_list.append(
                self.context.getState(getEnergy=True)
                .getPotentialEnergy()
                .value_in_unit(unit.kilojoules_per_mole)
            )
        return angles, np.array(e_list)


class OrbMMHybridFactory:
    def get_rotatable_torsions_rdkit(self, mol):
        import rdkit.Chem as Chem

        pattern = Chem.MolFromSmarts("[!#1;!D1]-[!#1;!D1]")
        if pattern is None:
            return []
        torsions = []
        for bond in mol.GetSubstructMatches(pattern):
            a1, a2 = bond[0], bond[1]
            for n1 in [x.GetIdx() for x in mol.GetAtomWithIdx(a1).GetNeighbors()]:
                if n1 != a2:
                    for n2 in [
                        x.GetIdx() for x in mol.GetAtomWithIdx(a2).GetNeighbors()
                    ]:
                        if n2 != a1 and mol.GetBondBetweenAtoms(a1, a2) is not None:
                            torsions.append((n1, a1, a2, n2))
        return torsions


# 用法（扫配体全部可转键二面角）：
#
#   from rdkit import Chem
#   mol = Chem.MolFromMolFile("ligand.sdf", removeHs=False)
#   torsions = OrbMMHybridFactory().get_rotatable_torsions_rdkit(mol)
#   scanner = OrbScanner(model_name="mace-off24-medium", charge=0)
#   for quad in torsions:
#       angles_deg, e_kj = scanner.run_torsion_scan(Chem.Mol(mol), quad)
#
# 注意 run_torsion_scan 会原地改 conformer，所以每条二面角传一份拷贝。
