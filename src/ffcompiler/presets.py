"""设计 §17 推荐默认 preset（逐字冻结；执行语义见 2026-11-11 计划 §4）。

表内字段与设计稿 YAML 一一对应；不可用项在编排层（cli._do_run）显式报错，
不静默替换：
- ``charge: new_charge`` 是 MVP-4（未立项）→ 报错并列出可用覆盖；
- ``bonded.hessian: orca``（high）读 .hess 文件 → 缺 --hessian-file 报错；
- ``rigid_torsion: {method: drih}``（rigid）是设计表内部不一致（DRIH 是 bonded
  方法，非二面角 provider）→ v1 不执行，原样存入 run_manifest 并标注待澄清。
"""

from __future__ import annotations

PRESETS: dict[str, dict] = {
    "fast": {
        "typing": "gaff2",
        "bonded": {"provider": "hessian", "hessian": "gfnff", "method": "mseminario"},
        "torsion": "mlp",
        "charge": "new_charge",
        "vdw": "gaff2",
    },
    "balanced": {
        "typing": "gaff2",
        "geometry": {"provider": "gfn2"},
        "bonded": {"provider": "hessian", "hessian": "gfn2", "method": "mseminario"},
        "torsion": "mlp",
        "charge": "new_charge",
        "vdw": "gaff2",
    },
    "high": {
        "geometry": {"provider": "orca"},
        "bonded": {"provider": "hessian", "hessian": "orca", "method": "mseminario"},
        "torsion": "mlp",
        "charge": "new_charge",
        "vdw": "gaff2",
    },
    "rigid": {
        "bonded": {"hessian": "gfn2", "method": "m2seminario"},
        "rigid_torsion": {"method": "drih"},
        "charge": "new_charge",
    },
}
