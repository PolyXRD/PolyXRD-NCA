"""仓库骨架自检（T0.1 验收）。"""
from __future__ import annotations


def test_version() -> None:
    from nca.version import APP_NAME, __version__

    assert __version__ == "0.1.0"
    assert APP_NAME == "PolyXRD NCA"


def test_package_importable() -> None:
    import nca
    import nca.core

    assert nca is not None and nca.core is not None
