# Crucible

[English](README.md) | [中文](README.zh-CN.md) | **日本語**

Crucible はエンジン非依存の低分子力場コンパイラです。1 回のパラメータ化の結果を共通の
中間表現 `ForceFieldIR` に保持し、**GROMACS**・**AMBER**・**OpenMM** の各形式へネイティブに
書き出します。書き出しのたびに 3 エンジン間のエネルギー/力の等価検証を自動で行います。

```text
AMBER prmtop/inpcrd + SDF
        │  読み込み（ネイティブ実装、ParmEd 不要）
        ▼
   ForceFieldIR ◄──── 座標 sidecar（全コマンド共通：--ir-in / --coords-in）
        │
        ├─ bond/angle ← mSeminario · m2Seminario · DRIH · 全体 Hessian フィット
        │                Hessian：Gaussian fchk · ORCA .hess · CP2K .out · xtb（GFN-FF / GFN2）
        ├─ 二面角     ← MACE-POLAR-1-M 拘束緩和スキャン → Fourier 項（GAFF2 フォールバック、理由を記録）
        └─ 電荷       ← 読み込み · AM1-BCC（antechamber）· NAGL
        │
        ▼
   ForceFieldIR（完全——全パラメータに provenance を付与）
        │
   ┌────┼─────────┐
   ▼    ▼         ▼
GROMACS  AMBER   OpenMM   → エンジン間エネルギー / 力の比較
```

**バージョン 0.2.0。** MVP-0～MVP-3 は完了し、Gate も通過しています（2026-10-18 に現在の
ソースで G0–G3 を再実行、すべて合格）。2026-10-13 以降の追加：ORCA `.hess` 読み込み、
CP2K `.out` 読み込み、全コマンド共有の IR + 座標パイプライン、`global-fit` と
`check-staleness` コマンド、周期ボックス対応、Hessian フォールバックの枠組み（MVP-5、
Python API のみ）、学習型 BCC の骨格（MVP-4、学習済みモデルなし）、およびオプションの
逐次二面角再フィット（`--refinement-rounds {0,1,2,3}`、既定値 0 = シングルパス）。
設計書との 34 項目の結案検収は
[ACCEPTANCE_20261018_CLOSURE](docs/history/ACCEPTANCE_20261018_CLOSURE.md)、プロジェクトの
権威ある細粒度の状態は [HANDOVER.md](HANDOVER.md)、変更点は [CHANGELOG](CHANGELOG.md)
を参照してください。

## 機能

### 1 つの IR、3 つのエンジン

`ForceFieldIR` は原子、bond、angle、proper/improper 二面角、LJ と電荷、除外リストと 1-4
ペア、および任意の周期ボックスを保持し、検証と NumPy による参照エネルギーを備えます。
3 つのエクスポータは GROMACS top/itp/gro、AMBER prmtop/inpcrd、OpenMM System/XML を
ネイティブに書き出します——ParmEd も中間フォーマットも経由しません。同じ IR を 3 エンジンに
書き出した場合、項ごとのエネルギーと原子力が一致しなければなりません。`check-equivalence`
がこの比較を行い（8 レッグ：NumPy、OpenMM ×4、sander ×2、GROMACS）、G0 Gate は
11 分子 × 4 構造で実行します。

### Hessian からの bonded パラメータ

Cartesian Hessian から bond/angle パラメータを導出する 4 手法が、共通のプロバイダ抽象と
単位系（内部は原子単位）を共有します：

- **mSeminario**——局所投影。汎用の既定候補；
- **m2Seminario**——対称性と方向の一致を重視した変体。剛体・高対称系向け；
- **DRIH**——冗長内部座標の対角成分。剛体・特殊な結合のフラグメント向け；
- **全体 Hessian フィット**（QMDFF 式）——固定部分（LJ + Coulomb + 1-4 + 二面角）を差し引いた
  上で全力定数を Hessian 全体にフィット。事前正則化付き。条件数・残差・x₀ での残留力を報告；
- **localhess**（実験的）——局所コンプライアンス手法。内部整合性のみ検証済み
  （決定 0023。Gate の判定基準には入れない）。

グラフ対称性による等価項の平均化が既定で有効（`--no-symmetrize` で無効化）。非正の力定数は
黙ってフォールバックせず例外を送出します。

### MLP からの二面角

分子グラフから回転可能結合を認識し、MACE-POLAR-1-M（ASE 経由）に対する拘束緩和スキャンを
行い、現在の力場の MM 寄与を差し引き、プロファイルを AMBER 周期項にフィットします。
GAFF2 へのフォールバックは結合ごとに理由を記録し、再フィットで黙って元のフォールバック結合を
昇格させることはありません（決定 0020）。オプションの有界逐次再フィット
（`--refinement-rounds {0,1,2,3}`）は各結合を他の結合の現在のフィットに対して再フィットします。
電荷変更時は記録されたプロトコルをキャッシュ済みスキャンからビット単位で再生します。

### 電荷プロバイダ

電荷の読み込み、AM1-BCC（antechamber/sqm）、NAGL（`openff-gnn-am1bcc-1.0.0`）が共通の
プロバイダインターフェースを共有し、総電荷を厳密に保存します。プロバイダを変更すると、
再スキャンせずキャッシュ済みスキャンから MLP 二面角を再フィットします。学習型 BCC の研究骨格
（結合反対称・電荷保存の補正）は整備済み。学習済みモデルはまだありません。RESP は後回し
（決定 0018）。

### パーセプション層と力場パラメータライブラリ

幾何ベースの結合認識層（Cordero 2008 半径 ×1.15 とデータファイルの距離ウィンドウ、線形
カットオフ付き angle/二面角生成、三斜 PBC の最近接イメージ、データ駆動の improper 中心）と、
Sobtop 形式 `LJ_param.dat` / `bonded_param.dat` / `bondcrit.dat` パラメータファイルの
読み込み/マッチャ（X ワイルドカード優先順位、IDIVF、負の周期、LJ 変換）が参照データの
ワークフローを支えます。決定 0022 と `docs/theory/` のリバースエンジニアリング報告を参照。

### Provenance と鮮度検査

全パラメータが `source` とメタデータを持ち、各フィットは入力（どの電荷、どの二面角版、どの
Hessian 識別情報）を項ごとのスナップショットとして記録します。フィットが記録した入力の
いずれかが変更されていれば `check-staleness` が終了コード 1 で失敗します。IR + 座標 sidecar は
コマンド間で幾何を忠実に受け渡します——原子の順序・元素・名前は決して並べ替えません。

## クイックスタート

既存の科学計算環境をそのまま使います。本プロジェクトは何もインストールしません。

- Python ≥ 3.12（`numpy`、`scipy`、`openmm`、`rdkit`）。MLP 二面角と NAGL 電荷の機能には
  さらに ASE + MACE、OpenFF NAGL が必要
- 外部プログラム（該当する機能・テストでのみ必要）：
  - `xtb` 6.7.1：既定は `~/software/xtb-dist/bin/xtb`、`FFC_XTB` で変更可能
  - GROMACS：`FFC_GMX`、`PATH`、`~/gmx26.3/bin/gmx` の順に探索。mixed 精度で問題なし
  - AmberTools `sander` / `tleap`：エンジン間比較とフィクスチャ生成。AM1-BCC 電荷は
    さらに `antechamber` を使用（`FFC_ANTECHAMBER` または `PATH`）
  - MACE-POLAR-1-M モデル（二面角スキャン用）：既定は
    `~/MLP/mace/MACE-POLAR-1-M.model`、`FFC_MACE_MODEL` で変更可能
  - Sobtop 2026.1.16：参照データの再生成時のみ（`scripts/make_sobtop_refs.py`）

```bash
source scripts/env.sh
```

パラメータ化済みの AMBER 分子を 3 エンジンへ変換：

```bash
PYTHONPATH=src python -m ffcompiler.cli convert \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --output gromacs amber openmm --outdir jobs/ethanol_conversion
```

GFN2 Hessian から bond/angle を再導出（xtb は `<outdir>/xtb/` で実行）：

```bash
PYTHONPATH=src python -m ffcompiler.cli derive-bonded \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --hessian-provider gfn2 --method mseminario \
  --output gromacs amber openmm --outdir jobs/ethanol_mseminario
```

- `--hessian-provider {gfnff,gfn2,fchk,orca-hess,cp2k-out}`（`fchk`、`orca-hess`、`cp2k-out`
  には `--hessian-file` が必要）
- `--method {mseminario,m2seminario,drih}`；`--method localhess` は上記の実験的手法
- `--no-symmetrize`：対称等価な項の平均化を行わない
- `--linear-angles {error,library}`：既定は `error`。`library` は 175° を超える角の元の
  パラメータを保持

MACE-POLAR 緩和スキャンから回転可能結合の二面角をフィット
（`<outdir>/scans/` にキャッシュ）：

```bash
PYTHONPATH=src python -m ffcompiler.cli derive-torsions \
  --prmtop tests/fixtures/gaff2/ethanol/ethanol.prmtop \
  --inpcrd tests/fixtures/gaff2/ethanol/ethanol.inpcrd \
  --sdf    tests/fixtures/gaff2/ethanol/ethanol.sdf \
  --output openmm --outdir jobs/ethanol_torsions
```

- `--step-deg`（既定 10°）、`--device`（既定 cpu）、`--scan-cache`（既定 `<outdir>/scans`）
- `--refinement-rounds {0,1,2,3}`：シングルパス後の追加逐次ラウンド。既定値 0 は
  シングルパスのまま

書き出したファイルが各エンジンで IR のエネルギーを再現するか確認：

```bash
PYTHONPATH=src python -m ffcompiler.cli check-equivalence \
  --prmtop ... --inpcrd ... --sdf ... --workdir jobs/ethanol_check
```

同じ IR 上で工程をつなぐ——各工程は `ir.json` + `coords.json` を書き出し、次の工程が
両方を読み込みます：

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
# フィット時の入力が後で変わっていれば終了コード 1
PYTHONPATH=src python -m ffcompiler.cli check-staleness --ir-in jobs/mol/b/ir.json
```

`global-fit` は QMDFF 風の全体 Hessian フィットです（Hessian の指定は同じ。`--prior`、
`--alpha`）。k = 0 の下限に当たるクラスがあると `ir_failed_lower_bound.json` を書き、
書き出しはせず終了コード 1。

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
  E_h/a₀² ≈ 9.37583 × 10⁵ kJ mol⁻¹ nm⁻²、角度の対角成分は
  E_h/rad² = 2625.50 kJ mol⁻¹ rad⁻² で換算します（角度に a₀² を掛け直さない）。
- 対応する Hamiltonian：気相、カットオフなし、Lorentz–Berthelot 混合則、一様な AMBER 1-4
  スケーリング。それ以外は読み込み時に拒否します。

## 検証状況

4 つの Gate は 2026-10-18（オプションの逐次再フィット組み込み後）に現在のソースで再実行され、
すべて合格しました。数値は 2026-10-13 の実行と同じで、各 Gate の pytest は 524 + 59 すべて
合格です。Gate の合否判定は既定のシングルパスに対するもので、逐次モードは単体テストと
ビット単位の再現チェック（`experiments/20261018_refinement_verify/`）で担保します。

| Gate | 実行 | 結果 |
|---|---|---|
| G0 — 書き出しの等価性 | `jobs/gate_mvp0_20261018_2/` | 合格：11 分子 × 4 構造、8 レッグ、2596 行、失敗 0 |
| G1 — Hessian 手法 | `jobs/gate_mvp1_20261018_2/` | 合格：3 手法すべてが 12/12 の参照ケースで Sobtop と一致。全体フィットの自己無撞着誤差 6.0 × 10⁻⁹ |
| G2 — MLP 二面角 | `jobs/gate_mvp2_20261018_2/` | 合格：12/12 の回転可能結合を MLP でフィット、フォールバック 0。sander 以外のエンジン間レッグはすべて合格。既知の制約（決定 0017）：sander が π 付近のフィット位相を丸める（nma 二面角 2.3 × 10⁻⁴ > 2 × 10⁻⁴ kcal/mol）、報告のみ |
| G3 — 電荷プロバイダ | `jobs/gate_mvp3_20261018_2/` | 合格：AM1-BCC が 12 フィクスチャの電荷を再現（最大差 4.3 × 10⁻¹¹ e）。3 プロバイダ × 11 分子で総電荷保存（≤ 8.3 × 10⁻¹⁷ e）、エンジン間比較はすべて合格。電荷変更後の MLP 二面角再フィット：電荷が同じなら完全一致、NAGL 電荷でも 12/12 結合が MLP のまま |

導出した力場の参照 Hessian に対する振動数 RMSE（11 分子の中央値、cm⁻¹。報告のみで
判定には使わない）：

| Hessian | GAFF2 | mSeminario | m2Seminario | DRIH | 全体フィット |
|---|---|---|---|---|---|
| GFN2 | 110 | 109 | 188 | 219 | 118 |
| GFN-FF | 101 | 95 | 76 | 203 | **57** |

Gate を自分で実行する場合（毎回新しい `jobs/gate_*_<日付>[_N]/` を作成し、過去の結果は残す）：

```bash
python scripts/gate_mvp0.py
python scripts/gate_mvp1.py --jobs 6   # xtb を 6 並列 × 各 4 スレッド
python scripts/gate_mvp2.py --jobs 8 --threads 4 [--resume 旧G2ディレクトリ]   # MACE スキャン
python scripts/gate_mvp3.py --scan-cache jobs/<今回のG2実行>/scans          # 必須
PYTHONPATH=src python -m pytest -m "not external" -q   # 単体テスト（外部プログラム不要）
```

## 既知の制限

- **既知の不具合（MVP-2）**：安息香酸メチルの結合 1-2 では、MLP でフィットした二面角が
  MACE のプロファイルを GAFF2 より悪く再現します（プロファイル RMSE 3.36 対 2.70 kJ/mol、
  フィット RMSE 0.60）。12 本の回転可能結合のうちこの 1 本だけです。2026-10-17 に診断済み
  （`jobs/torsion_diagnosis_20261017/report.md`）：結合ごとのフィットで他の結合を GAFF2 の
  まま差し引くため、1-2 の項が硬すぎる GAFF2 のメチル（0-1）項を吸収し、0-1 を置き換えると
  その補償が外れます。`--refinement-rounds 3` で 0.44 kJ/mol、統合残差は 0.83 → 0.22 まで
  下がります（サンプル内、`jobs/torsion_b3_sequential_20261017/report.md`）。どの更新順序も
  同じ不動点に収束し、3 ラウンド後の順序による差は ≤ 1e-4 kJ/mol です
  （`jobs/torsion_b3_order_20261018/report.md`）。オプションのまま：既定はシングルパス、
  ホールドアウト検証（B4）は未実施、逐次モードは Gate の判定基準に入っていません。
- **スキャン品質（MVP-2）**：nma 1-3（アミド C–N）の往復スキャンのヒステリシスは最大
  57 kJ/mol（36 フレーム中 14 フレームが 1 kJ/mol 超、±60–120° の障壁域に集中）。これは
  フィットではなくスキャンの問題で、再スキャンするかは未決、自動フォールバック閾値も
  設けていません。
- 1 本の結合の二面角エネルギーはすべて 1 つの二面角四つ組に載せます（MVP-2 計画 §7 第 3 項）。
  その適用範囲はまだ確定していません。
- ORCA `.hess` 読み込みは ORCA 6.1.1 の閉殻サンプル 1 件でのみ検証しています
  （振動数の相対差 5.7 × 10⁻⁸）。
- 直線に近い角（> 175°）は導出しません。アセトニトリルは診断用ケースのみです。
- m2Seminario は GFN-FF Hessian で負の環角力定数を与えることがあります（シクロプロパン）。
  Sobtop も同じ値を返します。`derive_bonded` は黙ってフォールバックせず
  `NonPositiveForceConstantError` を送出するので、その場合は mSeminario か DRIH を使ってください。
- mSeminario の bond 力定数は結合の記述方向に依存します（Sobtop と同じ挙動）。既定で有効な
  対称化により、その非対称性は解消されます。
- 書き出しは単一残基のみ。IR は任意の周期ボックスを保持できます（`box_nm`、決定 0023）：
  NumPy エネルギーは最近接イメージ + 結合グラフのアンラップ、GROMACS/OpenMM は完全な
  三斜ボックスを書き出し、AMBER prmtop は直交系のみ。周期系のエンジン間レッグはまだ
  check-equivalence に入っていません。

## リポジトリ構成

```text
src/ffcompiler/
  ir/          ForceFieldIR スキーマ、検証、NumPy エネルギー、JSON 入出力、座標 sidecar
  perception/  幾何 bond/angle/dihedral/improper 認識、Cordero 半径、PBC ミラー
  fflib/       Sobtop 形式 LJ/bonded/bondcrit ライブラリの読み込みとワイルドカード照合
  importers/   AMBER prmtop の読み込み
  exporters/   GROMACS・AMBER・OpenMM の書き出し
  formats/     prmtop / rst7 のネイティブ読み書き
  hessian/     CartesianHessian、fchk / ORCA .hess / CP2K .out / xtb プロバイダ、フォールバック
  bonded/      内部座標、mSeminario、m2Seminario、DRIH、localhess、対称化、全体フィット、鮮度検査
  torsion/     回転可能結合、MACE-POLAR 緩和スキャン、MM 差し引き、Fourier フィット、MLP プロバイダ
  charge/      電荷プロバイダ（読み込み、AM1-BCC、NAGL）、学習型 BCC の骨格、書き戻しと二面角の再フィット
  validation/  エンジン間比較、GROMACS アダプタ、振動解析、pytest 集計
scripts/       Gate、フィクスチャ生成、Sobtop 参照データ
tests/         単体テストとエンジン間テスト、フィクスチャ
docs/          計画、決定記録、理論文書
jobs/          正式実行（manifest/acceptance/report のみ追跡）
```

## ドキュメント

- [forcefield_compiler_design.md](forcefield_compiler_design.md)——設計書（検収の基準。
  次フェーズの節 §3.1、§6.4.4–6.4.6、§10.4–10.5 は「設計済み・未実装」と明記）
- [docs/history/ACCEPTANCE_20261018_CLOSURE.md](docs/history/ACCEPTANCE_20261018_CLOSURE.md)——
  現在のバージョンの結案検収（設計 34 項目との突合。§2 に既知の不具合・決定待ち項目・
  未実施の検証）
- [docs/theory/hessian_methods.md](docs/theory/hessian_methods.md)——導出、単位、規約
- [docs/decisions/](docs/decisions/README.md)——決定記録 0001–0024
- [docs/superpowers/plans/](docs/superpowers/plans/)——ロードマップと MVP 計画
- [HANDOVER.md](HANDOVER.md)——現在の検証済み状態と次の作業
- [CHANGELOG.md](CHANGELOG.md)

## ロードマップ

**完了し Gate 済み**：MVP-0（IR + エクスポータ）、MVP-1/1b（Hessian 手法 + 全体フィット）、
MVP-2（MLP 二面角）、MVP-3（電荷プロバイダ）——G0–G3 が現在のソースで合格。

**骨格のみ、研究は未着手/未完了**：MVP-4 学習型 BCC（骨格 + 研究提案のみ、学習は未承認）、
MVP-5 フォールバック制御（階段は整備済み、不確実性駆動の制御は未着手）。

**設計済み、未実装**（2026-10-25 設計改訂、決定 0024——設計書 §3.1、§6.4.4–6.4.6、
§10.4–10.5）：

- **MVP-6 分子エントリ**——SDF を入力し骨格 IR を得る：antechamber を外部 GAFF2 typer として
  使用（typing のみ）、GAFF2 ライブラリ値を事前値/フォールバックに、構造の黙って生成はしない。
  単一コマンド CLI と preset（§17/§22）への道を開く。
- **MVP-7 スキャン品質**——マルチスタート緩和スキャンと格子点ごとの分岐診断。分岐は診断のみに
  使いフィット対象にはしない。nma 1-3 の再スキャンが最初の適用例。
- **MVP-8 平衡構造ジョイントフィット**——q⁰ をフィット変数に加え勾配マッチングを追加
  （§6.4.3-1 が確保していたインターフェース）。オプトインで、現行の全体フィットに
  ビット単位で退化。
- **MVP-9 多配座ジョイントフィット + 交互更新の収束**——同じ目的関数で C ≥ 1、q⁰ 共有、
  ホールドアウト配座検証、および charge → torsions → global-fit 交互更新の ε 許容の
  完了判定。
- TorsionDrive はオプトインの 2 次元連成スキャン用に確保。本番環境へインストールするかは
  明示的な未決事項。

**未実施の検証**：B4 二面角再フィットのホールドアウト検証、本物の MACE 参照による C5 の
正式な組み合わせ検収、Level 3–4 ベンチマーク（設計 §20）。

**決定待ち**：既定の bonded 手法、既定の再フィット回数、バージョン番号、TorsionDrive の
インストール——HANDOVER.md と結案報告 §2 を参照。
