# Crucible

[English](README.md) | [中文](README.zh-CN.md) | **日本語**

エンジン非依存の低分子力場コンパイラです。パラメータ化済みの分子を共通の中間表現
`ForceFieldIR` に読み込み、必要に応じて量子化学 Hessian からパラメータを再導出し、
**GROMACS**・**AMBER**・**OpenMM** の各形式へネイティブに書き出します。

**バージョン 0.2.0**：MVP-0（IR と 3 種のエクスポータ）、MVP-1（Hessian からの bond/angle
パラメータ導出）、MVP-2（MLP 二面角、OpenMM 基準・決定 0017）、MVP-3（電荷プロバイダ：
読み込み、AM1-BCC、NAGL）が完了し、2026-10-18 に現在のソースで G0–G3 を再実行して
すべて合格しています。2026-10-13 以降の追加：ORCA `.hess` 読み込み、全コマンドで共有する
IR + 座標のパイプライン（`--ir-in` / `--coords-in`）、`global-fit` と `check-staleness`
コマンド、Hessian フォールバックの枠組み（MVP-5、Python API のみ）、学習型 BCC の骨格
（MVP-4、学習済みモデルなし）、およびオプションの逐次二面角再フィット
（`derive-torsions --refinement-rounds {0,1,2,3}`、既定値 0 = シングルパス）。
2026-10-18 の設計書との突合（34 項目、現行の逸脱一覧）は
[ACCEPTANCE_20261018_CLOSURE](docs/history/ACCEPTANCE_20261018_CLOSURE.md) を参照。
変更点は [CHANGELOG](CHANGELOG.md) を参照してください。

## 機能

| 領域 | 内容 |
|---|---|
| IR | `ForceFieldIR`：原子、bond、angle、proper/improper 二面角、LJ と電荷、除外リストと 1-4 ペア。検証と NumPy による参照エネルギー |
| 読み込み | AMBER prmtop/inpcrd + SDF（形式電荷）。ネイティブ実装で ParmEd 不要 |
| 書き出し | GROMACS top/itp/gro、AMBER prmtop/inpcrd、OpenMM System/XML をすべてネイティブに出力 |
| Hessian | Gaussian fchk、ORCA `.hess`、xtb 6.7.1（GFN-FF、GFN2）。収束・電荷・品質のチェック付き。GFN-FF → GFN2 → DFT ファイルのフォールバック枠組み（Python API のみ） |
| 結合項の導出 | mSeminario、m2Seminario、DRIH、全体 Hessian フィット。分子グラフの対称性で平均化 |
| 二面角（MVP-2） | 回転可能結合の認識、MACE-POLAR-1-M（ASE 経由）による拘束緩和スキャン、MM 差し引き、AMBER 周期項へのフィット。失敗時は GAFF2 に戻し理由を記録。オプションの逐次再フィット回数（`ffcompiler derive-torsions`） |
| 電荷（MVP-3） | 電荷の読み込み、AM1-BCC（antechamber/sqm）、NAGL（openff-gnn-am1bcc-1.0.0）。総電荷を厳密に保存。電荷を変えるとキャッシュ済みスキャンから MLP 二面角を再フィット（`ffcompiler assign-charges`）。RESP は未対応 |
| パイプライン | IR と座標 sidecar をコマンド間で受け渡し（`--ir-in` / `--coords-in`）。フィットの依存を記録し古くなったかを検査（`ffcompiler check-staleness`） |
| 検証 | エンジン間のエネルギー/力比較（NumPy、OpenMM ×4、sander ×2、GROMACS）と 4 つの受け入れ Gate（G0–G3） |

## 動作環境

既存の科学計算環境をそのまま使います。本プロジェクトは何もインストールしません。

- Python ≥ 3.12（`numpy`、`scipy`、`openmm`、`rdkit`）。MLP 二面角と NAGL 電荷の機能には
  さらに ASE + MACE、OpenFF NAGL（モデル `openff-gnn-am1bcc-1.0.0`）
- 外部プログラム（該当する機能・テストでのみ必要）：
  - `xtb` 6.7.1：既定は `~/software/xtb-dist/bin/xtb`、`FFC_XTB` で変更可能
  - GROMACS：`FFC_GMX`、`PATH`、`~/gmx26.3/bin/gmx` の順に探索。mixed 精度で問題なし
  - AmberTools の `sander` / `tleap`：エンジン間比較とフィクスチャ生成。AM1-BCC 電荷は
    さらに `antechamber` を使用（`FFC_ANTECHAMBER` または `PATH`）
  - 二面角スキャン用 MACE-POLAR-1-M モデル：既定は `~/MLP/mace/MACE-POLAR-1-M.model`、`FFC_MACE_MODEL` で変更可能
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

- `--hessian-provider {gfnff,gfn2,fchk,orca-hess}`（`fchk` と `orca-hess` には `--hessian-file` が必要）
- `--method {mseminario,m2seminario,drih}`
- `--no-symmetrize`：対称等価な項の平均化を行わない
- `--linear-angles {error,library}`：既定は `error`。`library` は 175° を超える結合角の元のパラメータを保持

MACE-POLAR の緩和スキャンから回転可能結合の二面角をフィット（スキャンは `<outdir>/scans/` にキャッシュ）：

```bash
PYTHONPATH=src python -m ffcompiler.cli derive-torsions \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --output openmm --outdir jobs/ethanol_torsions
```

- `--step-deg`（既定 10°）、`--device`（既定 cpu）、`--scan-cache`（既定 `<outdir>/scans`）
- `--refinement-rounds {0,1,2,3}`：シングルパスの後に逐次更新を追加実行し、各結合を他の結合の
  現在のフィット項に対して再フィットします。既定値 0 はシングルパスのまま。電荷変更時
  （`assign-charges`）は記録したプロトコルをキャッシュ済みスキャンで再現します。

書き出したファイルが各エンジンで IR のエネルギーを再現するか確認：

```bash
PYTHONPATH=src python -m ffcompiler.cli check-equivalence \
  --prmtop ... --inpcrd ... --sdf ... --workdir jobs/ethanol_check
```

同じ IR 上で工程をつなぐ：各工程は `ir.json` と `coords.json`（その工程で使った／返した構造、
例えば Hessian の構造）を書き出し、次の工程は `--ir-in` / `--coords-in` で読み込みます
（原子番号・元素・原子名が一致しなければ拒否、並べ替えはしません）：

```bash
# 1) 電荷
PYTHONPATH=src python -m ffcompiler.cli assign-charges --prmtop P --inpcrd C --sdf S \
  --method am1bcc --output openmm --outdir jobs/mol/q
# 2) Hessian から bond/angle を導出（1 の続き）
PYTHONPATH=src python -m ffcompiler.cli derive-bonded \
  --ir-in jobs/mol/q/ir.json --coords-in jobs/mol/q/coords.json \
  --hessian-provider gfn2 --method mseminario --output openmm --outdir jobs/mol/b
# 3) 最終 IR と座標を書き出し
PYTHONPATH=src python -m ffcompiler.cli convert \
  --ir-in jobs/mol/b/ir.json --coords-in jobs/mol/b/coords.json \
  --output gromacs amber openmm --outdir jobs/mol/final
# フィット時の入力（電荷・二面角など）が後で変わっていれば終了コード 1
PYTHONPATH=src python -m ffcompiler.cli check-staleness --ir-in jobs/mol/b/ir.json
```

`global-fit` は QMDFF 風の全体 Hessian フィットです（Hessian の指定は同じ。`--prior`、`--alpha`）。
k = 0 の下限に当たるクラスがあると `ir_failed_lower_bound.json` を書き、書き出しはせず終了コード 1。


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

2026-10-18 に現在のソース（オプションの逐次再フィット組み込み後）で 4 つの Gate を再実行し、
すべて合格しました。数値は 2026-10-13 の実行と同じで、各 Gate の pytest は 524 + 59 すべて合格です。
Gate の合否判定は既定のシングルパスに対するもので、逐次モードは単体テストとビット単位の
再現チェック（`experiments/20261018_refinement_verify/`）で担保します。

| Gate | 実行ディレクトリ | 結果 |
|---|---|---|
| G0：書き出しの等価性 | `jobs/gate_mvp0_20261018_2/` | 合格：11 分子 × 4 構造、8 本の比較レッグ、2596 行、失敗 0 |
| G1：Hessian 手法 | `jobs/gate_mvp1_20261018_2/` | 合格：3 手法すべてが 12/12 の参照ケースで Sobtop と一致。全体フィットの自己無撞着誤差 6.0 × 10⁻⁹ |
| G2：MLP 二面角 | `jobs/gate_mvp2_20261018_2/` | 合格：12/12 の回転可能結合を MLP でフィット、フォールバック 0。sander 以外のエンジン間比較はすべて合格。既知の制約（決定 0017）：sander が π 付近のフィット位相を π に丸める（nma 二面角 2.3 × 10⁻⁴ > 2 × 10⁻⁴ kcal/mol）、報告のみ |
| G3：電荷プロバイダ | `jobs/gate_mvp3_20261018_2/` | 合格：AM1-BCC が 12 フィクスチャの電荷を再現（最大差 4.3 × 10⁻¹¹ e）。3 プロバイダ × 11 分子で総電荷保存（≤ 8.3 × 10⁻¹⁷ e）、エンジン間比較はすべて合格。電荷変更後の MLP 二面角再フィット：電荷が同じなら完全一致、NAGL 電荷でも 12/12 結合が MLP のまま |

導出した力場の振動数 RMSE（参照 Hessian 比、11 分子の中央値、cm⁻¹。報告のみで合否判定には使わない）：

| Hessian | GAFF2 | mSeminario | m2Seminario | DRIH | 全体フィット |
|---|---|---|---|---|---|
| GFN2 | 110 | 109 | 188 | 219 | 118 |
| GFN-FF | 101 | 95 | 76 | 203 | **57** |

Gate を自分で実行する場合（毎回新しい `jobs/gate_*_<日付>[_N]/` を作成し、過去の結果は残します）：

```bash
python scripts/gate_mvp0.py
python scripts/gate_mvp1.py --jobs 6   # xtb を 6 並列 × 各 4 スレッド
python scripts/gate_mvp2.py --jobs 8 --threads 4 [--resume 旧G2実行ディレクトリ]   # MACE スキャン（--resume 指定時はその scans/ を再利用）
python scripts/gate_mvp3.py --scan-cache jobs/<今回のG2実行>/scans          # 必須：今回受け入れる G2 の scans/
PYTHONPATH=src python -m pytest -m "not external" -q   # 単体テスト（外部プログラム不要）
```

## 既知の制限

- **既知の不具合（MVP-2）**：安息香酸メチルの結合 1-2 では、MLP でフィットした二面角が MACE の
  プロファイルを GAFF2 より悪く再現します（プロファイル RMSE 3.36 vs 2.70 kJ/mol、フィット RMSE 0.60）。
  12 本の回転可能結合のうちこの 1 本だけです。2026-10-17 に診断済み（`jobs/torsion_diagnosis_20261017/report.md`）：
  結合ごとのフィットで他の結合を GAFF2 のまま差し引くため、1-2 の項が硬すぎる GAFF2 のメチル（0-1）項を吸収し、
  0-1 を置き換えるとその補償が外れます。`--refinement-rounds 3` で 0.44 kJ/mol、統合残差は
  0.83 → 0.22 まで下がります（サンプル内、`jobs/torsion_b3_sequential_20261017/report.md`）。
  どの更新順序も同じ不動点に収束し、3 ラウンド後の順序による差は ≤ 1e-4 kJ/mol です
  （`jobs/torsion_b3_order_20261018/report.md`）。このモードはオプションのままです：既定は
  シングルパス、B4 ホールドアウト検証は未実施、逐次モードは Gate の判定基準に入っていません。
- **スキャン品質（MVP-2）**：nma 1-3（アミド C–N）の往復スキャンのヒステリシスは最大 57 kJ/mol
  （36 フレーム中 14 フレームが 1 kJ/mol 超、±60–120° の障壁域に集中）。これはフィットではなく
  スキャンの問題で、再スキャンするかは未決、自動フォールバック閾値も設けていません。
- 1 本の結合の二面角エネルギーはすべて 1 つの二面角四つ組に載せます（MVP-2 計画 §7 第 3 項）。
  その適用範囲はまだ確定していません。
- ORCA `.hess` 読み込みは ORCA 6.1.1 の閉殻サンプル 1 件でのみ検証しています
  （振動数の相対差 5.7 × 10⁻⁸）。
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
  ir/          ForceFieldIR のスキーマ、検証、NumPy エネルギー、JSON 入出力、座標 sidecar
  importers/   AMBER prmtop の読み込み
  exporters/   GROMACS・AMBER・OpenMM への書き出し
  formats/     prmtop / rst7 のネイティブ読み書き
  hessian/     CartesianHessian、fchk / ORCA .hess / xtb のプロバイダ、フォールバック
  bonded/      内部座標、mSeminario、m2Seminario、DRIH、対称化、全体フィット、依存の鮮度検査
  torsion/     回転可能結合、MACE-POLAR 緩和スキャン、MM 差し引き、Fourier フィット、MLP プロバイダ
  charge/      電荷プロバイダ（読み込み、AM1-BCC、NAGL）、学習型 BCC の骨格、書き戻しと二面角の再フィット
  validation/  エンジン間比較、GROMACS アダプタ、振動解析、pytest 集計
scripts/       Gate、フィクスチャ生成、Sobtop 参照データ
tests/         単体テストとエンジン間テスト、フィクスチャ
docs/          計画、決定記録、理論文書
jobs/          正式実行（manifest/acceptance/report のみ追跡）
```

## ドキュメント

- [docs/history/ACCEPTANCE_20261018_CLOSURE.md](docs/history/ACCEPTANCE_20261018_CLOSURE.md)：
  現在のバージョンの結案検収（設計書 34 項目との突合。§2 に既知の不具合・決定待ち項目・未実施の検証）
- [docs/theory/hessian_methods.md](docs/theory/hessian_methods.md)：導出、単位、規約
- [docs/decisions/](docs/decisions/README.md)：決定記録 0001–0021
- [docs/superpowers/plans/](docs/superpowers/plans/)：ロードマップと各フェーズの計画
- [HANDOVER.md](HANDOVER.md)：現在の検証済み状態と次の作業
- [CHANGELOG.md](CHANGELOG.md)

## ロードマップ

MVP-2 MLP 二面角は完了、Gate G2 合格。MVP-3 電荷プロバイダは完了、Gate G3 合格（RESP は後回し、
決定 0018）。MVP-4 学習型 BCC はプロバイダの骨格と研究提案のみで、学習は未承認。MVP-5 は
フォールバックの枠組みのみで、制御の研究は未着手です。結案報告 §2 に従い、次は B4 ホールドアウト
検証、本物の MACE 参照による C5 の正式な組み合わせ検収、MVP-4 の立案検討。その後は
molecule/typing の入口と preset/単一コマンド CLI（設計 §1/§17/§22）。既定の bonded 手法・
既定の再フィット回数・バージョン番号などの決定待ち項目はそこに列挙してあります。
リポジトリは CRUCIBLE として非公開でホストしています。
