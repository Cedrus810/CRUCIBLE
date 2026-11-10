# 决策 0026：外部工具位置由用户自己的配置文件给出（tools.toml + init-config + check-env）

2026-10-28。计划：`docs/superpowers/plans/2026-10-28-tool-config.md`。

## 背景

10-25 的发布提交为了脱敏，把工具的默认路径改成了相对 HOME（`~/software/xtb-dist`、
`~/gmx26.3`、`~/MLP/mace/...`）。这样做默认每个用户的目录布局都一样，换一个账号就会出问题：
10-28 用另一个账号跑 Gate 时，G1 的 xtb 和一个 CLI 测试都因为 `~` 解析到别人的家目录而失败。
另外，各工具的查找方式不统一：有的支持环境变量，有的只认 PATH，`scripts/gate_mvp1.py`
完全不认 `FFC_XTB`。

## 用户原话（2026-10-28）

- 「是这样 本质是为了脱敏。。。。」
- 「但是 这么说 给别人交付的时候 应该给人家一个config 然后他们用户可以自己写位置啊」
- 「所以搞一个 config 就很好用啊。。。叫其他人自己写？」
- 「然后我们允许用户缺项然后他自己可以缺功能」
- 对计划及两处补充（缺 sander 时跳过、check-env 输出功能表）：「可以 走起」

## 裁定

1. 新模块 `ffcompiler/toolconfig.py`。每个工具按以下顺序查找：环境变量 `FFC_<KEY>` → 配置文件
   `[tools]` → PATH → 旧的相对 HOME 默认位置（保留它，原来的用法不写配置也照常能跑）。
   已有的环境变量名 `FFC_XTB`、`FFC_GMX`、`FFC_ANTECHAMBER`、`FFC_MACE_MODEL` 保持不变。
2. 配置文件放在 `FFC_CONFIG` 指定的路径，没指定时用 `$XDG_CONFIG_HOME/crucible/tools.toml`
   （默认 `~/.config/crucible/tools.toml`），不入库。仓库里只有 `crucible.example.toml`，
   路径全部是 `/path/to/...` 占位。不支持放在仓库目录下的配置文件：那样得改 `.gitignore`，
   一旦漏加，真实路径就会被提交出去。
3. 配置文件里的路径：`~` 会展开；相对路径按配置文件所在目录解析，不按当前工作目录。
4. 环境变量或配置文件里写了路径、但这个路径不可用时，直接报错（ToolConfigError），不悄悄换成
   PATH 上的另一个版本。环境变量设成空字符串，当作没设。
5. xtb 不查 PATH（MVP-1 计划 §1 固定使用 6.7.1，沿用原来的规则）。
6. 缺某个工具时，只有用到它的功能不可用；报错信息会列出试过的每一处位置，并提示运行
   `init-config`。CLI `check-equivalence` 缺 sander 时跳过 sander 两路，记入 `report.skipped`；
   `compare()` 的 `skip_missing_sander` 参数默认是 False，所以 Gate（全套验收）缺 sander 仍然不通过。
7. `ffcompiler init-config`：把能自动找到的工具先填进去，找不到的写成注释掉的占位行；
   文件已存在时，不加 `--force` 就不覆盖。`ffcompiler check-env`：打印工具表（状态、来源、
   路径、xtb/gmx 的版本）和"可用功能"表；只在配置写错时返回 1，缺工具只报告、返回 0。
8. `scripts/env.sh` 的 conda 根目录、环境名、modulefiles 目录和 gmx 模块名都可以用
   `FFC_MAMBA_ROOT`、`FFC_ENV_NAME`、`FFC_MODULEFILES`、`FFC_GMX_MODULE` 覆盖，默认值不变。
9. `scripts/make_sobtop_refs.py` 只在 `tools/sobtop/` 还没复制时才需要 Sobtop 源目录
   （原来每次都检查源目录）。

## 验证（2026-10-28）

只靠配置文件（非原作者账号、不改 HOME、不设 `FFC_<工具>`）：`pytest -m "not external"` 603 passed，全量 external 60 passed；G0–G3（`jobs/gate_mvp*_20261028_3/`）全部通过，report.md 与决策 0025 修复后那一轮逐字相同。

## 不做的事

- 运行记录（manifest、provenance、Gate 报告）里会写下工具的绝对路径，这是可追溯性要求的；
  入库前怎么脱敏另外决定。
- NAGL 模型跟着 Python 环境走（从 `openff.nagl_models` 包里取），不进配置。ORCA、Gaussian、CP2K
  本项目只读它们的输出文件，不调用程序，也不进配置。
