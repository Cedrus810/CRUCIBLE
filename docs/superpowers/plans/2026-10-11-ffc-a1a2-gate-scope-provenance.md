# A1/A2 收尾：Gate 覆盖范围与参数来源追踪 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 子集 Gate 不再能显示通过；`--resume` 不改写旧运行目录；电荷写回有中央检查、来源可追溯；重拟合不丢原兜底原因；`--reference mm` 的扫描缓存带上参考势身份。

**Architecture:** 只改现有 Gate 脚本（`scripts/gate_mvp1/2/3.py`）、`charge/apply.py`、`torsion/mlp.py`、`cli.py` 及对应单元测试。不加新模块，不改数值容差，不重跑 Gate（重跑属于 A3，要等用户提交后另报计划）。

**Tech Stack:** Python 3.12（openmm_dev，生产环境，不装不改）、pytest、numpy、OpenMM（MM 参考势）。`PYTHONPATH=src`。

**Spec:** `forcefield_compiler_design.md`（§4.1 每个参数必须带 provenance；§6.4.2 provenance 必须记录依赖的 charge/torsion 版本）；
`docs/superpowers/plans/2026-10-11-crucible-next-steps.md` §3 A1/A2；决策 0017、0018、0019。

## 与 2026-10-11 讨论结论不同的一处（请确认）

**与你同意的不同：** 我之前建议"`--reference mm` 的参考能干脆不缓存"。读了代码以后，这个办法不可行。
**原因：** 换电荷重拟合的测试（`tests/unit/test_torsion_refit.py`、`tests/unit/test_cli.py::test_derive_torsions_scan_cache_then_assign_charges_refit`）就是用 MM 参考扫描结果的缓存，来代替 MACE 缓存做快速测试的。不缓存，这两条测试和 CLI 的"derive → 换电荷 → 重拟合"链路都会断。
**改为（T42）：** 把参考势的身份写进缓存键：`mm-self` 记录扫描时 IR 的电荷，用逐项字面值，不用 hash（全局规则 §8.4）。torsion provenance 也同时记一份。重拟合时从 provenance 取扫描时记下的身份，所以仍然命中原缓存。也就是说，参考势在扫描时就冻结了，和 MACE 的行为一致。对另一套电荷的 IR 重新 `derive-torsions --reference mm` 时，缓存键对不上，就会重扫。
MACE 路径的缓存键**一字不变**（`reference=None` 时不往键里写这个字段），所以 `jobs/gate_mvp2_20261008_3/scans` 仍然有效。

## Global Constraints

- NFS 上不跑任何 git；每个任务完成后不 commit，最后合写一份 `.commit_drafts/` 草稿（全局 仓库代理规程 §1）。
- openmm_dev 是生产环境：不 pip/conda 安装、卸载、升级。
- 不新增数值阈值或收敛判据；现有容差（`TOL_CHARGE = 1e-10` 等）、决策 0017（sander 只报告）不变。
- 不用 hash/sha256 判断缓存能否复用（全局 §8.4）。
- 不删、不改 `jobs/` 下任何已有运行目录；`jobs/gate_mvp2_20261008_3/scans` 被 G2/G3 共用。
- 测试只跑与改动相关的文件，`-m "not external"`；这里涉及的测试都不调 MACE 和外部程序。
- 测试命令统一前缀：`cd /home/user/CRUCIBLE && source scripts/env.sh && PYTHONPATH=src python -m pytest`

## Review Focus

1. 只跑部分 provider 或分子（G1 `--providers gfnff`、G2/G3 `--molecules ethanol`）时，即使各判据都满足，`passed` 也必须为 false，并且 `scope == "subset"`（T38）。
2. 用 `--resume` 跑完后，旧运行目录 `scans/` 下每个文件的字节和 mtime 都不变；`.tmp.npz` 半成品不会被复制（T39）。
3. 自定义 provider 返回的电荷多一个、少一个或含 NaN 时，`assign_charges` 报 `ChargeError`，不截断，也不写回（T40）。
4. 某根键原来是 `mlp_error` 兜底，换电荷重拟合后，报告里仍然是 `mlp_error`，而不是 `no_cached_scan`；重拟合不改写缓存文件（T41）。
5. 已有的 MACE 扫描缓存键在 T42 前后逐字节相同（T42 `test_scan_cache_key_without_reference_is_unchanged`）。

## 文件清单

| 文件 | 改动 |
|---|---|
| `scripts/gate_mvp1.py`、`gate_mvp2.py`、`gate_mvp3.py` | `build_report` 加 `full_matrix`；acceptance 加 `scope`（T38） |
| `scripts/gate_mvp2.py` | `seed_scan_cache`；`--resume` 改为复制进本次运行目录（T39） |
| `src/ffcompiler/charge/apply.py` | 形状/有限性检查；`provenance["charge"]["previous"]`（T40） |
| `src/ffcompiler/torsion/mlp.py` | 重拟合保留原兜底行（T41）；`reference` 进缓存键与 provenance（T42） |
| `src/ffcompiler/cli.py` | `--reference mm` 传 `reference`（T42） |
| `tests/unit/test_gate_mvp{1,2,3}_report.py`、`test_charge_apply.py`、`test_torsion_refit.py`、`test_torsion_mlp.py`、`test_cli.py` | 对应用例 |
| `HANDOVER.md`、`CHANGELOG.md`、`docs/superpowers/plans/2026-10-11-crucible-next-steps.md`、`.superpowers/sdd/2026-10-11-ffc-a1a2/progress.md`、`.commit_drafts/` | 文档与草稿（T43） |

---

### Task 38: Gate 子集运行不能通过（决策 0019）

**Files:**
- Modify: `scripts/gate_mvp1.py`（`build_report` 约 :118，acceptance 约 :250，`main` 约 :578）
- Modify: `scripts/gate_mvp2.py`（`build_report` :120，acceptance :199，`main` :396）
- Modify: `scripts/gate_mvp3.py`（`build_report` :50，acceptance :81，`main` :225）
- Test: `tests/unit/test_gate_mvp1_report.py`、`test_gate_mvp2_report.py`、`test_gate_mvp3_report.py`

**Interfaces:**
- Produces: 三个 `build_report(..., full_matrix: bool)`（关键字参数，必填）；acceptance 新键 `"scope": "full" | "subset"`。

- [ ] **Step 1: 改测试 helper，并写失败用例**

三个测试文件的 `_build` 默认参数里加 `full_matrix=True`（G1：`kwargs = dict(full_matrix=True, pytest_summary=...)`；G2、G3：`args = dict(full_matrix=True, ...)`）。然后在每个文件末尾加（`acc` 取法按各文件：G1/G2 是 `_build(...)[2]`，G3 是 `_build(...)[1]`）：

```python
def test_subset_run_never_passes() -> None:
    # 决策 0019：各判据全部满足，但只跑了子集 → passed=false、scope=subset
    report_md, *_, acc = _build(full_matrix=False)
    assert acc["scope"] == "subset"
    assert acc["passed"] is False
    assert all(acc["criteria"].values())
    assert "scope=subset" in report_md


def test_full_run_records_scope() -> None:
    acc = _build()[-1]
    assert acc["scope"] == "full" and acc["passed"] is True
```

G1 的 `_build` 返回 `(report_md, csv, acc)`，G3 的返回 `(report_md, acc)`；上面写成 `report_md, *_, acc` 和 `[-1]`，三个文件都能用。

- [ ] **Step 2: 跑测试，确认失败**

Run: `... pytest tests/unit/test_gate_mvp1_report.py tests/unit/test_gate_mvp2_report.py tests/unit/test_gate_mvp3_report.py -q`
Expected: 新用例 FAIL（`TypeError: build_report() got an unexpected keyword argument 'full_matrix'`）。

- [ ] **Step 3: 实现**

三个 `build_report` 的签名都加关键字参数 `full_matrix: bool`（放在 `*` 之后）。在报告标题后插入提示，并改 acceptance。以 G2 为例（G1 的报告列表名是 `report_lines`，G3 的是 `lines`，写法相同）：

```python
    lines = [REPORT_HEADER]
    if not full_matrix:
        lines.append("\n> 子集运行（scope=subset）：判据结果只供调试，不是正式通过（决策 0019）。\n")
    lines += ["\n## 每根可转键\n", ...]  # 原有内容不变
    ...
    acceptance = {
        "gate": "G2",
        "scope": "full" if full_matrix else "subset",
        "passed": full_matrix and all(criteria.values()),
        ...
    }
```

G3 的 passed 是 `full_matrix and all(criteria.values()) and not failures`。

`main` 里传入：
- G1：`full_matrix=set(args.providers) == {"gfnff", "gfn2"}`
- G2：`full_matrix=set(names) == set(NAMES)`
- G3：`full_matrix=set(args.molecules) == set(G31_NAMES)`

三个脚本的 `--molecules` / `--providers` help 文字末尾加"（子集运行 passed 恒为 false，决策 0019）"。

- [ ] **Step 4: 跑测试，确认通过**

Run: 同 Step 2。Expected: 全部 PASS（原有用例也要过）。

---

### Task 39: G2 `--resume` 不再写旧运行目录

**Files:**
- Modify: `scripts/gate_mvp2.py`（`main` :307 附近；新增 `seed_scan_cache`）
- Test: `tests/unit/test_gate_mvp2_report.py`

**Interfaces:**
- Produces: `seed_scan_cache(old_root: Path, new_root: Path, names: Iterable[str]) -> list[str]`（返回已复制的 `分子/文件名`）。

- [ ] **Step 1: 写失败用例**

```python
def test_seed_scan_cache_copies_without_touching_old(tmp_path) -> None:
    old = tmp_path / "old" / "scans"
    (old / "ethanol").mkdir(parents=True)
    (old / "dmso").mkdir()
    (old / "ethanol" / "1-2.npz").write_bytes(b"scan-a")
    (old / "ethanol" / "1-2.tmp.npz").write_bytes(b"half")
    (old / "dmso" / "0-1.npz").write_bytes(b"scan-b")
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in old.rglob("*") if p.is_file()}
    new = tmp_path / "new" / "scans"
    copied = gate.seed_scan_cache(old, new, ["ethanol"])
    assert copied == ["ethanol/1-2.npz"]
    assert (new / "ethanol" / "1-2.npz").read_bytes() == b"scan-a"
    assert not (new / "ethanol" / "1-2.tmp.npz").exists()
    assert not (new / "dmso").exists()
    after = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in old.rglob("*") if p.is_file()}
    assert after == before
```

- [ ] **Step 2: 跑测试，确认失败**

Run: `... pytest tests/unit/test_gate_mvp2_report.py -q -k seed`
Expected: FAIL（`AttributeError: module 'gate_mvp2' has no attribute 'seed_scan_cache'`）。

- [ ] **Step 3: 实现**

`scripts/gate_mvp2.py` 顶部 import 加 `import shutil`，以及 `from typing import Iterable`（如未导入）。新增：

```python
def seed_scan_cache(old_root: Path, new_root: Path, names: Iterable[str]) -> list[str]:
    """--resume：把旧运行目录 scans/<分子>/*.npz 复制进本次运行的 scans/；旧目录只读。
    键不一致的文件照常由 _scan_job 判为未命中并重扫，新结果只写本次目录。"""
    copied: list[str] = []
    for name in names:
        src_dir = Path(old_root) / name
        if not src_dir.is_dir():
            continue
        dst_dir = Path(new_root) / name
        dst_dir.mkdir(parents=True, exist_ok=True)
        for src in sorted(src_dir.glob("*.npz")):
            if src.name.endswith(".tmp.npz"):
                continue
            shutil.copy2(src, dst_dir / src.name)
            copied.append(f"{name}/{src.name}")
    return copied
```

`main` 中把

```python
        cache_root = (args.resume / "scans") if args.resume else (run_dir / "scans")
        ctx.manifest["scan_cache"] = str(cache_root)
```

改为（`names` 的定义要提到这之前）：

```python
        names = tuple(args.molecules)
        cache_root = run_dir / "scans"
        ctx.manifest["scan_cache"] = str(cache_root)
        if args.resume:
            seeded = seed_scan_cache(args.resume / "scans", cache_root, names)
            ctx.manifest["scan_cache_seeded_from"] = {"dir": str(args.resume / "scans"), "files": seeded}
```

并删掉后面原有的 `names = tuple(args.molecules)`。模块 docstring 第 11 行与 `--resume` 的 help 改成"复制旧运行目录的 scans/ 到本次目录后复用（旧目录不写）"。数据量：每个扫描文件约 13 KB，12 根键约 160 KB。

- [ ] **Step 4: 跑测试，确认通过**

Run: `... pytest tests/unit/test_gate_mvp2_report.py -q`。Expected: PASS。

---

### Task 40: 电荷写回的中央检查与来源链

**Files:**
- Modify: `src/ffcompiler/charge/apply.py`
- Test: `tests/unit/test_charge_apply.py`、`tests/unit/test_torsion_refit.py`（加强弱断言）

**Interfaces:**
- Consumes: `ChargeError`（`charge/base.py`，已有）。
- Produces: `ir.provenance["charge"]["previous"]`：上一条 `provenance["charge"]`（没有就是 `None`），逐层嵌套。

- [ ] **Step 1: 写失败用例**

`tests/unit/test_charge_apply.py` 末尾（`_Fixed`、`_load` 已在文件里）：

```python
from ffcompiler.charge.base import ChargeError
from ffcompiler.charge.loaded import LoadedCharges


@pytest.mark.parametrize("delta", [1, -1])
def test_wrong_length_from_provider_raises(delta: int) -> None:
    ir, pos, mol = _load("ethanol")
    with pytest.raises(ChargeError, match="形状"):
        assign_charges(ir, _Fixed(np.zeros(len(ir.atoms) + delta)), pos, mol)


def test_non_finite_charge_raises() -> None:
    ir, pos, mol = _load("ethanol")
    q = np.zeros(len(ir.atoms))
    q[3] = np.nan
    with pytest.raises(ChargeError, match="非有限"):
        assign_charges(ir, _Fixed(q), pos, mol)


def test_reloading_keeps_previous_charge_source() -> None:
    ir, pos, mol = _load("ethanol")
    assert "charge" not in ir.provenance
    first, _ = assign_charges(ir, _Fixed(np.zeros(len(ir.atoms))), pos, mol)
    assert first.provenance["charge"]["previous"] is None
    again, _ = assign_charges(first, LoadedCharges(), pos, mol)
    prov = again.provenance["charge"]
    assert prov["method"] == "loaded"
    assert prov["previous"]["method"] == "fixed-test"
    assert prov["previous"]["source"] == {"kind": "test"}
```

（如果实现时发现 `from_amber` 的 IR 已经带 `provenance["charge"]`，就把第一个断言改成记下它，再断言 `first.provenance["charge"]["previous"]` 等于记下的值。）

`tests/unit/test_torsion_refit.py::test_refit_with_new_charges_changes_terms`：原来的 `charge_method == "loaded"` 在重拟合之前就成立，测不出更新。改成用一个 method 不同的 provider：

```python
class _Shifted(ChargeProvider):
    method = "shifted-test"

    def __init__(self, q):
        self.q = q

    def compute(self, ir, positions_nm, mol):
        return self.q.copy()

    def source_info(self):
        return {"kind": "test"}


def test_refit_with_new_charges_changes_terms(scanned, tmp_path: Path) -> None:
    _, mlp_ir, _, pos, mol, cache = _mlp_ethanol(scanned, tmp_path)
    assert mlp_ir.provenance["torsion"]["charge_method"] == "loaded"
    q = np.array([a.partial_charge for a in mlp_ir.atoms])
    q[0] += 0.05
    q[1] -= 0.05
    new_ir, report = assign_charges(mlp_ir, _Shifted(q), pos, mol, scan_cache=cache)
    rows = new_ir.provenance["torsion"]["rotors"]
    assert [r["status"] for r in rows] == ["mlp", "mlp"]
    assert new_ir.proper_torsions != mlp_ir.proper_torsions
    assert new_ir.provenance["torsion"]["charge_method"] == "shifted-test"
    assert report.torsion_refit["rotors"] == rows
```

（文件头 import 加 `from ffcompiler.charge.base import ChargeProvider`。）

- [ ] **Step 2: 跑测试，确认失败**

Run: `... pytest tests/unit/test_charge_apply.py tests/unit/test_torsion_refit.py -q`
Expected: 长度/NaN 用例 FAIL（DID NOT RAISE）；`previous` 用例 FAIL（KeyError）；refit 用例按现有代码应当 PASS，因为它测的是本来就正确的行为，这里只是把断言加强。

- [ ] **Step 3: 实现**

`src/ffcompiler/charge/apply.py`：

```python
from ffcompiler.charge.base import ChargeError, ChargeProvider, ChargeReport, equalize_residual
...
    raw = np.asarray(provider.compute(ir, positions_nm, mol), dtype=np.float64)
    n = len(ir.atoms)
    if raw.shape != (n,):
        raise ChargeError(f"{provider.method}: provider 返回的电荷形状为 {raw.shape}，IR 有 {n} 个原子")
    if not np.all(np.isfinite(raw)):
        raise ChargeError(f"{provider.method}: provider 返回了非有限电荷（NaN/Inf）")
    charges, residual = equalize_residual(raw, ir.total_charge)
    ...
    atoms = tuple(replace(a, partial_charge=float(q)) for a, q in zip(ir.atoms, charges, strict=True))
    prov = {"method": report.method, "raw_sum": report.raw_sum, "residual": report.residual,
            "source": report.source, "previous": ir.provenance.get("charge")}
```

docstring 补一句："provider 返回的形状/有限性在这里统一检查；provenance['charge']['previous'] 保留上一次的电荷来源"。

- [ ] **Step 4: 跑测试，确认通过**

Run: `... pytest tests/unit/test_charge_apply.py tests/unit/test_charge_loaded.py tests/unit/test_torsion_refit.py -q`。Expected: PASS。

---

### Task 41: 重拟合保留原兜底原因，不改写缓存

**Files:**
- Modify: `src/ffcompiler/torsion/mlp.py`（`refit_torsions` :221–251）
- Test: `tests/unit/test_torsion_refit.py`

**Interfaces:**
- Consumes: `MLPTorsionProvider.parameterize`、`ScanError`（`torsion/scan.py`）。
- Produces: 重拟合报告与 `provenance["torsion"]["rotors"]` 中，原来不是 mlp 的键原样沿用原行。

- [ ] **Step 1: 写失败用例**

```python
from ffcompiler.torsion.scan import ScanError


def test_refit_keeps_original_fallback_reason_and_cache_bytes(tmp_path: Path) -> None:
    mdir = FIXTURES / "ethanol"
    ir, pos = from_amber(*(mdir / f"ethanol.{e}" for e in ("prmtop", "inpcrd", "sdf")))
    mol = Chem.MolFromMolFile(str(mdir / "ethanol.sdf"), removeHs=False)
    base, _ = assign_charges(ir, LoadedCharges(), pos, mol)
    calls = []

    def factory(ir_):
        calls.append(1)
        if len(calls) == 2:
            raise ScanError("boom")  # 第二根键模拟参考势出错 → mlp_error
        return OpenMMCalculator(ir_)

    cache = tmp_path / "scans"
    mlp_ir, report = MLPTorsionProvider(factory, supported=ALL_Z, scan_kwargs=FAST,
                                        scan_cache=cache).parameterize(base, pos, mol)
    reasons = {tuple(r["bond"]): r["reason"] for r in report.rotors}
    assert sorted(r["status"] for r in report.rotors) == ["fallback", "mlp"]
    failed = next(b for b, why in reasons.items() if why == "mlp_error")
    before = {p.name: p.read_bytes() for p in cache.iterdir()}

    new_ir, chg = assign_charges(mlp_ir, LoadedCharges(), pos, mol, scan_cache=cache)
    rows = {tuple(r["bond"]): r for r in chg.torsion_refit["rotors"]}
    assert rows[failed]["reason"] == "mlp_error"
    assert {tuple(r["bond"]): r["reason"] for r in new_ir.provenance["torsion"]["rotors"]}[failed] == "mlp_error"
    assert {p.name: p.read_bytes() for p in cache.iterdir()} == before
```

- [ ] **Step 2: 跑测试，确认失败**

Run: `... pytest tests/unit/test_torsion_refit.py -q -k original_fallback`
Expected: FAIL（`'no_cached_scan' == 'mlp_error'`）。

- [ ] **Step 3: 实现**

`refit_torsions` 在 `lost` 检查之后、`return` 之前加：

```python
    # 原来就不是 mlp 的键（mlp_error、scan_not_converged 等）沿用原行，不改记为 no_cached_scan
    original = {frozenset(r["bond"]): r for r in prov.get("rotors", [])}
    kept = [dict(original[frozenset(r["bond"])])
            if r["status"] != "mlp" and original.get(frozenset(r["bond"]), {}).get("status") not in (None, "mlp")
            else r
            for r in report.rotors]
    report.rotors = kept
    torsion_prov = {**new_ir.provenance["torsion"], "rotors": kept}
    new_ir = new_ir.replace_terms(provenance_update={"torsion": torsion_prov})
    return new_ir, report
```

（`replace_terms` 的 `provenance_update` 是按顶层键整体替换，`ir/schema.py:111-113`，所以要传完整的 torsion 字典。）

- [ ] **Step 4: 跑测试，确认通过**

Run: `... pytest tests/unit/test_torsion_refit.py tests/unit/test_torsion_mlp.py -q`。Expected: PASS。

---

### Task 42: `--reference mm` 的缓存键带上参考势身份

**Files:**
- Modify: `src/ffcompiler/torsion/mlp.py`（`scan_cache_key` :78，`MLPTorsionProvider.__init__` :124、`parameterize` :172、`_finish` :208，`refit_torsions` :242）
- Modify: `src/ffcompiler/cli.py`（`_do_derive_torsions` :110–117）
- Test: `tests/unit/test_torsion_mlp.py`、`tests/unit/test_cli.py`

**Interfaces:**
- Produces: `scan_cache_key(..., reference: Mapping | None = None)`；`MLPTorsionProvider(..., reference: Mapping | None = None)`；`provenance["torsion"]["reference"]`（`None`，或 `{"kind": "mm-self", "charges": [...]}`）。

- [ ] **Step 1: 写失败用例**

`tests/unit/test_torsion_mlp.py`：

```python
import json

from ffcompiler.torsion.mlp import scan_cache_key


def test_scan_cache_key_without_reference_is_unchanged() -> None:
    # 保护已有 MACE 缓存（jobs/gate_mvp2_20261008_3/scans）：reference=None 时键与旧公式逐字相同
    kw = dict(model_file={"path": "m", "size": 1, "mtime": 2.0}, numbers=[6, 1], bonds=[(1, 0)])
    old = json.dumps({"model": "M", "model_file": {"path": "m", "size": 1, "mtime": 2.0},
                      "quartet": [0, 1, 2, 3], "info": {"charge": 0, "spin": 1},
                      "scan": {"step_deg": 10.0}, "numbers": [6, 1], "bonds": [[0, 1]]}, sort_keys=True)
    assert scan_cache_key("M", (0, 1, 2, 3), {"charge": 0, "spin": 1}, {"step_deg": 10.0}, **kw) == old
    with_ref = scan_cache_key("M", (0, 1, 2, 3), {"charge": 0, "spin": 1}, {"step_deg": 10.0}, **kw,
                              reference={"kind": "mm-self", "charges": [0.1, -0.1]})
    assert json.loads(with_ref)["reference"] == {"kind": "mm-self", "charges": [0.1, -0.1]}
```

`tests/unit/test_cli.py`：

```python
def test_mm_reference_cache_is_frozen_for_refit_and_rescanned_for_other_charges(tmp_path: Path) -> None:
    import json

    import numpy as np

    from ffcompiler.ir import io as ir_io

    cache = tmp_path / "scans"
    assert main(["assign-charges", *_ethanol_args(), "--method", "loaded",
                 "--output", "openmm", "--outdir", str(tmp_path / "q0")]) == 0
    assert main(["derive-torsions", *_ethanol_args(), "--ir-in", str(tmp_path / "q0" / "ir.json"),
                 "--reference", "mm", "--step-deg", "20", "--scan-cache", str(cache),
                 "--output", "openmm", "--outdir", str(tmp_path / "t0")]) == 0
    q0 = [a.partial_charge for a in ir_io.from_json(tmp_path / "q0" / "ir.json").atoms]
    q2 = np.array(q0)
    q2[0] += 0.05
    q2[1] -= 0.05
    qfile = tmp_path / "q2.txt"
    qfile.write_text("\n".join(repr(float(x)) for x in q2))
    # 1) 换电荷重拟合：参考势在扫描时冻结，仍命中缓存
    assert main(["assign-charges", *_ethanol_args(), "--ir-in", str(tmp_path / "t0" / "ir.json"),
                 "--method", "loaded", "--charges-file", str(qfile), "--scan-cache", str(cache),
                 "--output", "openmm", "--outdir", str(tmp_path / "r2")]) == 0
    # 2) 对另一套电荷的 IR 重新 derive：参考势不同 → 重扫，键里记新电荷
    assert main(["assign-charges", *_ethanol_args(), "--method", "loaded", "--charges-file", str(qfile),
                 "--output", "openmm", "--outdir", str(tmp_path / "q2")]) == 0
    assert main(["derive-torsions", *_ethanol_args(), "--ir-in", str(tmp_path / "q2" / "ir.json"),
                 "--reference", "mm", "--step-deg", "20", "--scan-cache", str(cache),
                 "--output", "openmm", "--outdir", str(tmp_path / "t2")]) == 0
    q2_ir = [a.partial_charge for a in ir_io.from_json(tmp_path / "q2" / "ir.json").atoms]
    for path in cache.glob("*.npz"):
        with np.load(path) as d:
            assert json.loads(str(d["key"]))["reference"] == {"kind": "mm-self", "charges": q2_ir}
    prov = ir_io.from_json(tmp_path / "t2" / "ir.json").provenance["torsion"]["reference"]
    assert prov == {"kind": "mm-self", "charges": q2_ir}
```

- [ ] **Step 2: 跑测试，确认失败**

Run: `... pytest tests/unit/test_torsion_mlp.py tests/unit/test_cli.py -q -m "not external" -k "reference"`
Expected: key 用例 FAIL（`unexpected keyword argument 'reference'`）；CLI 用例 FAIL（`KeyError: 'reference'`）。

- [ ] **Step 3: 实现**

`src/ffcompiler/torsion/mlp.py`：

```python
def scan_cache_key(model: str, quartet, info: Mapping, scan_kwargs: Mapping, *,
                   model_file: Mapping | None, numbers, bonds, reference: Mapping | None = None) -> str:
    """扫描缓存键：模型名与文件身份、四元组、电荷/spin、扫描参数、原子序数与键（计划 §2.7）；
    reference：参考势依赖 IR 时（mm-self）的身份，逐项字面值。为 None 时不写入，已有 MACE 缓存键不变。"""
    payload = {"model": model, "model_file": dict(model_file) if model_file else None,
               "quartet": list(quartet), "info": dict(info),
               "scan": {k: scan_kwargs[k] for k in sorted(scan_kwargs)},
               "numbers": [int(z) for z in numbers],
               "bonds": sorted([min(int(a), int(b)), max(int(a), int(b))] for a, b in bonds)}
    if reference is not None:
        payload["reference"] = dict(reference)
    return json.dumps(payload, sort_keys=True)
```

`MLPTorsionProvider.__init__` 加关键字参数 `reference: Mapping | None = None`，存成 `self.reference = dict(reference) if reference is not None else None`；docstring 加"reference：参考势依赖 IR 时的身份（mm-self 记扫描时的电荷），进缓存键与 provenance；参考势在扫描时冻结"。
`parameterize` 里 `scan_cache_key(...)` 调用加 `reference=self.reference`。
`_finish` 的 `prov` 加 `"reference": self.reference`。
`refit_torsions` 构造 provider 时加 `reference=prov.get("reference")`。

`src/ffcompiler/cli.py` 的 mm 分支：

```python
        provider = MLPTorsionProvider(lambda ir_: OpenMMCalculator(ir_),
                                      supported=frozenset(range(1, 119)),
                                      scan_kwargs={"step_deg": args.step_deg}, model_name="mm-self",
                                      scan_cache=scan_cache,
                                      reference={"kind": "mm-self",
                                                 "charges": [a.partial_charge for a in ir.atoms]})
```

注释改为："调试/测试用：参考势取输入 IR 自身的 MM，在扫描时冻结（键里记电荷）；之后换电荷重拟合仍用这次的参考能"。

已知范围：身份里只记电荷。改 bonded/LJ 后再跑 mm-self 不会被识别为不同参考势。目前 CLI 没有"bonded → IR → torsion"这条链（`derive-bonded` 不输出 IR JSON），等工作包 C 整合管线时，再按设计 §6.4.2 一并处理。另外，对同一个 `--scan-cache` 目录换参考势重扫时，会覆盖同名的 `j-k.npz`（现有行为，本计划不改）。

- [ ] **Step 4: 跑测试，确认通过**

Run: `... pytest tests/unit/test_torsion_mlp.py tests/unit/test_torsion_refit.py tests/unit/test_cli.py -q -m "not external"`
Expected: PASS（含原有 `test_derive_torsions_scan_cache_then_assign_charges_refit`）。

---

### Task 43: 文档、台账、提交草稿

**Files:**
- Create: `.superpowers/sdd/2026-10-11-ffc-a1a2/progress.md`
- Modify: `HANDOVER.md`、`CHANGELOG.md`、`docs/superpowers/plans/2026-10-11-crucible-next-steps.md`
- Create: `.commit_drafts/20261011_<时间>_a1a2-gate-scope-provenance.md`

- [ ] **Step 1: 相关测试合跑一次**

Run: `... pytest tests/unit/test_gate_mvp1_report.py tests/unit/test_gate_mvp2_report.py tests/unit/test_gate_mvp3_report.py tests/unit/test_charge_apply.py tests/unit/test_charge_loaded.py tests/unit/test_torsion_refit.py tests/unit/test_torsion_mlp.py tests/unit/test_cli.py -m "not external" -q`
Expected: 全部 PASS。把命令和输出尾部写进台账。

- [ ] **Step 2: 更新文档**

- `HANDOVER.md`：状态里写 A1/A2 已完成（附测试结果）；待办第 3 条里已修的 Minor 划掉，未做的继续保留；写明"正式 Gate 未重跑，按 A3 等用户提交后另报计划"。
- `CHANGELOG.md`：加 Unreleased 条目，按 T38–T42 各写一行。
- `2026-10-11-crucible-next-steps.md` §3：A1/A2 各条后面标"已做（本计划 T3x）"或"未做"。
- 台账逐任务记 complete 和测试结果，偏离计划的地方记成 Ruling。

- [ ] **Step 3: 写提交草稿**

`.commit_drafts/` 新建一份草稿，内容包括 commit message、文件清单（每个文件一句话），以及一条可以直接执行的 `git add … && git commit -m "…"`，覆盖本计划和决策 0019 的全部文件（包括 `docs/decisions/0019-gate-subset-scope.md`、`docs/decisions/README.md`、本计划文件、`2026-10-11-crucible-next-steps.md`）。

---

## 本计划不做

- A3：重跑 G0–G3，以及 G1.8 理论审阅（需要用户先提交，并另报计划表）。
- A4 的低影响遗留：`--scan-cache` 错误提示、重复的 `find_rotors`、mol2 报错、`step_deg` 提前检查、线程上限时机、ASE 异常分类、甲铵测试、SDF 报错、可转键类报告。
- 同一个缓存目录换参考势时覆盖同名文件（见 T42）。
- 任何 Gate 判据的数值或范围改动。
