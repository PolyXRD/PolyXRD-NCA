"""pytest 公共 fixture —— 数据根、实样路径与 Qt 离屏环境（T0.1 / T6.1b）。

约定：
- DATA_ROOT 指向家族只读数据根 D:/Project/XRD（COD 源、实样谱）。
- 只读引用，测试不写回数据根。
- GUI 测试跑在 **Qt 离屏平台**（offscreen）上，不弹窗、不抢焦点。
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# 让 tests/ 直接 import nca（无需安装包）
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

DATA_ROOT = Path("D:/Project/XRD")
SAMPLE_811 = DATA_ROOT / "XRData/csuHJW/txt/Data040：NCM：811_raw.txt"
CIF_DIR = ROOT / "nca" / "data" / "cif"


@pytest.fixture(scope="session")
def data_root() -> Path:
    """家族只读数据根。"""
    return DATA_ROOT


@pytest.fixture(scope="session")
def sample_811() -> Path:
    """真实 NCA_80_15_05 实样谱（岛津导出两列文本）。"""
    if not SAMPLE_811.exists():
        pytest.skip(f"实样缺失: {SAMPLE_811}")
    return SAMPLE_811


@pytest.fixture(scope="session")
def cif_dir() -> Path:
    """内置 CIF 目录。"""
    return CIF_DIR


# --------------------------------------------------------------------------- #
# Qt 离屏环境（T6.1b）
# --------------------------------------------------------------------------- #

#: 离屏平台名（必须在 QApplication 创建**之前**设入环境变量才生效）
QT_PLATFORM = "offscreen"


@pytest.fixture(scope="session")
def qapp():
    """会话级 QApplication（离屏）。

    QApplication 整个进程只能有一个，重复创建会崩；故用 session 作用域。
    平台强制 offscreen，测试机上没有显示服务器也能跑。
    """
    os.environ.setdefault("QT_QPA_PLATFORM", QT_PLATFORM)
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture(scope="session", autouse=True)
def _close_qt_windows_at_session_end():
    """会话结束时关掉所有顶层窗口 —— 否则 QThread 会带着运行中的线程被销毁。

    **这是个真缺陷，不是测试洁癖**：`TaskRunner` 的 QThread 生命周期绑在窗口上，
    窗口被 GC 而线程还在跑时，Qt 打印 `QThread: Destroyed while thread is still
    running` 并**直接 abort 进程**（实测 `test_gui_export` / `test_gui_main_window`
    全部用例通过后进程以 127 退出）。

    修法两层：应用侧 `MainWindow.closeEvent` 调 `flow.shutdown()`（已有），
    测试侧在这里统一收网 —— 只针对真的建过 QApplication 的会话，纯核心测试
    不会被拉进 Qt。
    """
    yield
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        return
    for widget in app.topLevelWidgets():
        # 测试替身可能没有 shutdown（如 test_tch_init 的 FakeRunner）——
        # 收尾阶段不能让替身的缺口变成「用例失败」，否则掩盖真正的问题。
        try:
            widget.close()
        except Exception:  # noqa: BLE001 —— 会话收尾，异常必须吞掉
            pass
    app.processEvents()
