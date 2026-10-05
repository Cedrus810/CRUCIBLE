"""tests 包标记文件。

openmm_dev 的 site-packages 里有一个无关的 ``tests`` 包；为了让 pytest 收到本仓库的
tests 而不是那个包，``tests/`` 及每个子目录都要有 ``__init__.py``（计划 §1.6）。
"""
