# Crucible

[English](README.md) | [中文](README.zh-CN.md) | **日本語**

エンジン非依存の低分子力場コンパイラです。パラメータ化済みの分子を共通の中間表現
`ForceFieldIR` に読み込み、必要に応じて量子化学 Hessian からパラメータを再導出し、
**GROMACS**・**AMBER**・**OpenMM** の各形式へネイティブに書き出します。

**バージョン 0.2.0**：MVP-0（IR と 3 種のエクスポータ）と MVP-1（Hessian からの bond/angle
パラメータ導出）が完了し、Gate G0・Gate G1 ともに合格しています。変更点は
[CHANGELOG](CHANGELOG.md) を参照してください。

## 機能

| 領域 | 内容 |
|---|---|
| IR | `ForceFieldIR`：原子、bond、angle、proper/improper 二面角、LJ と電荷、除外リストと 1-4 ペア。検証と NumPy による参照エネルギー |
| 読み込み | AMBER prmtop/inpcrd + SDF（形式電荷）。ネイティブ実装で ParmEd 不要 |
| 書き出し | GROMACS top/itp/gro、AMBER prmtop/inpcrd、OpenMM System/XML をすべてネイティブに出力 |
| Hessian | Gaussian fchk、xtb 6.7.1（GFN-FF、GFN2）。収束・電荷・品質のチェック付き |
| 結合項の導出 | mSeminario、m2Seminario、DRIH、全体 Hessian フィット。分子グラフの対称性で平均化 |
| 検証 | エンジン間のエネルギー/力比較（NumPy、OpenMM ×4、sander ×2、GROMACS）と 2 段階の受け入れ Gate |

## 動作環境

既存の科学計算環境をそのまま使います。本プロジェクトは何もインストールしません。

- Python ≥ 3.12（`numpy`、`scipy`、`openmm`、`rdkit`）
- 外部プログラム（該当する機能・テストでのみ必要）：
  - `xtb` 6.7.1：既定は `/home/user/software/xtb-dist/bin/xtb`、`FFC_XTB` で変更可能
  - GROMACS：`FFC_GMX`、`PATH`、`/home/user/gmx26.3/bin/gmx` の順に探索。mixed 精度で問題なし
  - AmberTools の `sander` / `tleap`：エンジン間比較とフィクスチャ生成
  - Sobtop 2026.1.16：参照データの再生成時のみ（`scripts/make_sobtop_refs.py`）

```bash
source scripts/env.sh
```

## クイックスタート

AMBER でパラメータ化された分子を 3 つのエンジン形式に変換：

```bash
PYTHONPATH=src python -m ffcompiler.cli convert \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --output gromacs amber openmm --outdir jobs/ethanol_conversion
```

GFN2 Hessian から bond と angle を再導出（xtb は `<outdir>/xtb/` で実行）：

```bash
PYTHONPATH=src python -m ffcompiler.cli derive-bonded \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --hessian-provider gfn2 --method mseminario \
  --output gromacs amber openmm --outdir jobs/ethanol_mseminario
```

- `--hessian-provider {gfnff,gfn2,fchk}`（`fchk` には `--hessian-file` が必要）
- `--method {mseminario,m2seminario,drih}`
- `--no-symmetrize`：対称等価な項の平均化を行わない
- `--linear-angles {error,library}`：既定は `error`。`library` は 175° を超える結合角の元のパラメータを保持

書き出したファイルが各エンジンで IR のエネルギーを再現するか確認：

```bash
PYTHONPATH=src python -m ffcompiler.cli check-equivalence \
  --prmtop ... --inpcrd ... --sdf ... --workdir jobs/ethanol_check
```

Python から使う場合：

```python
from ffcompiler.importers import from_amber
from ffcompiler.hessian.xtb import run_xtb_hessian
from ffcompiler.bonded import derive_bonded

ir, pos = from_amber("mol.prmtop", "mol.inpcrd", "mol.sdf")
h = run_xtb_hessian([a.atomic_number for a in ir.atoms], pos,
                    method="gfn2", charge=ir.total_charge, workdir="jobs/xtb")
new_ir, geom, report = derive_bonded(ir, h, method="drih", input_positions_nm=pos)
```

## 規約

- **IR の単位**：nm、rad、kJ mol⁻¹、電気素量。
- **調和項**は E = ½ k Δq² で保持し、k は Hessian の対角成分そのものです。GROMACS と
  OpenMM は k をそのまま使い、AMBER は K = k/2（E = K Δq²）に kcal・Å の換算を加えます。
- **DRIH** は内部で原子単位（Bohr、E_h/a₀²）を使います。伸縮の対角成分は
  E_h/a₀² ≈ 9.37583 × 10⁵ kJ mol⁻¹ nm⁻²、角度の対角成分は E_h/rad² = 2625.50 kJ mol⁻¹ rad⁻²
  で換算します（角度に a₀² を掛け直さない）。
- 対応する Hamiltonian：気相、カットオフなし、Lorentz–Berthelot 混合則、一様な AMBER 1-4
  スケーリング。それ以外は読み込み時に拒否します。

## 検証状況

| Gate | 実行ディレクトリ | 結果 |
|---|---|---|
| G0：書き出しの等価性 | `jobs/gate_mvp0_20261005_3/` | 合格：11 分子 × 4 構造、8 本の比較レッグ、2596 行、失敗 0 |
| G1：Hessian 手法 | `jobs/gate_mvp1_20261005_4/` | 合格：3 手法すべてが 12/12 の参照ケースで Sobtop と一致。全体フィットの自己無撞着誤差 6.0 × 10⁻⁹ |

導出した力場の振動数 RMSE（参照 Hessian 比、11 分子の中央値、cm⁻¹。報告のみで合否判定には使わない）：

| Hessian | GAFF2 | mSeminario | m2Seminario | DRIH | 全体フィット |
|---|---|---|---|---|---|
| GFN2 | 110 | 109 | 188 | 219 | 118 |
| GFN-FF | 101 | 95 | 76 | 203 | **57** |

Gate を自分で実行する場合（毎回新しい `jobs/gate_*_<日付>[_N]/` を作成し、過去の結果は残します）：

```bash
python scripts/gate_mvp0.py
python scripts/gate_mvp1.py --jobs 6   # xtb を 6 並列 × 各 4 スレッド
PYTHONPATH=src python -m pytest -m "not external" -q   # 単体テスト（外部プログラム不要）
```

## 既知の制限

- 直線に近い結合角（> 175°）は導出しません。アセトニトリルは診断用ケースのみです。
- m2Seminario は GFN-FF Hessian で負の環角力定数を与えることがあります（シクロプロパン）。
  Sobtop も同じ値を返します。`derive_bonded` は黙ってフォールバックせず
  `NonPositiveForceConstantError` を送出するので、その場合は mSeminario か DRIH を使ってください。
- mSeminario の bond 力定数は結合の記述方向に依存します（Sobtop と同じ挙動）。既定で有効な
  対称化により、その非対称性は解消されます。
- 書き出しは単一残基のみ。周期系や他の混合則は対象外です。

## ディレクトリ構成

```text
src/ffcompiler/
  ir/          ForceFieldIR のスキーマ、検証、NumPy エネルギー、JSON 入出力
  importers/   AMBER prmtop の読み込み
  exporters/   GROMACS・AMBER・OpenMM への書き出し
  formats/     prmtop / rst7 のネイティブ読み書き
  hessian/     CartesianHessian、fchk と xtb のプロバイダ
  bonded/      内部座標、mSeminario、m2Seminario、DRIH、対称化、全体フィット
  validation/  エンジン間比較、GROMACS アダプタ、振動解析、pytest 集計
scripts/       Gate、フィクスチャ生成、Sobtop 参照データ
tests/         単体テストとエンジン間テスト、フィクスチャ
docs/          計画、決定記録、理論文書
jobs/          正式実行（manifest/acceptance/report のみ追跡）
```

## ドキュメント

- [docs/theory/hessian_methods.md](docs/theory/hessian_methods.md)：導出、単位、規約
- [docs/decisions/](docs/decisions/README.md)：決定記録 0001–0015
- [docs/superpowers/plans/](docs/superpowers/plans/)：ロードマップと各フェーズの計画
- [HANDOVER.md](HANDOVER.md)：現在の検証済み状態と次の作業
- [CHANGELOG.md](CHANGELOG.md)

## ロードマップ

MVP-2 MLP 二面角 · MVP-3 電荷プロバイダ · MVP-4 学習型 BCC · MVP-5 フォールバック制御は
未着手です。リポジトリは CRUCIBLE として非公開でホストしています。
