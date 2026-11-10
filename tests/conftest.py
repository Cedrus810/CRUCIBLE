"""tests 包标记文件。

openmm_dev 的 site-packages 里有一个无关的 ``tests`` 包；为了让 pytest 收到本仓库的
tests 而不是那个包，``tests/`` 及每个子目录都要有 ``__init__.py``（计划 §1.6）。
"""

import pytest


@pytest.fixture(autouse=True)
def _isolate_tool_config(monkeypatch, tmp_path_factory, request):
    """单测不读开发者本机的 ~/.config/crucible/tools.toml（决策 0026）；external 测试要用真工具，不隔离。"""
    if request.node.get_closest_marker("external") is None:
        empty = tmp_path_factory.mktemp("toolcfg") / "tools.toml"
        empty.write_text("[tools]\n", encoding="utf-8")
        monkeypatch.setenv("FFC_CONFIG", str(empty))
