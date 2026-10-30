"""T8.2 · 多型号泛化 + 真实前驱体实样验证。

**数据现状如实登记**：家族数据根（`D:/Project/XRD/XRData/`）里**没有任何
NCA 实测谱**。唯一可用的三元层状正极实样是 `csuHJW/Data040：NCM：811`
（NCM811，Ni:Co:Mn = 8:1:1）—— 它与 NCA 同为 R-3m #166 层状、晶胞几乎
相同，因此**只能作同构锚点**：能证明管线在该骨架上收敛、晶胞/尺寸/占位
量级合理，**不能**证明模型化学正确（TM 位是 Al 还是 Mn，XRD 分不出，
但散射因子不同 → 系统性残差）。这一点由 Rwp 的 0.6pp 代价量化，见 C 路。
另有：

- `csuHJW/Data043：Ni(OH)2` —— 真实 **(Ni,Co,Mn)(OH)₂ 三元氢氧化物前驱体**
  实样（Ni:Co:Mn = 5:3:2），且**随附 Maud 2.55 的独立精修结果**
  （`Data043：Ni(OH)2.par`：P-3m1，a=3.130435 Å、c=4.619966 Å、
  z_O=0.24166）。这是本卡最有价值的一项 —— 真实数据 + 外部独立参考值。
- `ShiJiKong-shXRD/Si#B##.../sim_123*.txt` —— Maud 的 **Si 模拟谱**，
  与 NCA 无关（已核实文件头 `_maud_comment Simulated_data`），排除。
- 其余非 NCA 实样（尖晶石/LFP/岩盐/锐钛矿/ZnO 等）已由 T5.2 样本库
  承担负样本判据，本卡不重复。

因此本卡采用**三路验证**（偏离原「多型号实样」设想，已登记）：

| 路 | 对象 | 证明什么 |
|---|---|---|
| A | 80:15 / 82:13 / 85:10 / 90:05 四配比**合成谱** | 管线不止对 80:15 成立（配比泛化） |
| B | 真实前驱体 Data043 | 判定门正确分流 + 精修晶胞与 Maud 独立值一致 |
| C | 真实 NCM811（Data040） | 同构锚点：T8.1 的改动没有把实样锚点带偏 |

**为什么配比泛化只能用合成谱**：实验室 XRD 对 NCA 各配比的晶胞几乎相同
（内置四份 CIF 晶胞逐位相同，a=2.8645、c=14.1610），峰强差也只有
0.03–0.14 %（见 `test_model_intensity_differences_are_sub_percent`），
Ni:Co:Al 比例须由 ICP/EDS 确定（见 `docs/STD_FREEZE.md` §1 边界声明）。
所以「配比」在本卡里的可验证含义是**四份 CIF 各自能收敛、且回收晶胞一致**
—— 而不是「不同 Ni 含量能被 XRD 分开」，后者在物理上不成立，不能当判据。
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from nca.core import constants as C
from nca.core.dtypes import Pattern
from nca.core.identify import identify
from nca.core.patterns import load_pattern
from nca.core.pipeline import default_options_for, run_pipeline
from nca.core.profile.cif import CifStructure
from nca.core.profile.forward import ProfileParams, compute_phase_profile
from nca.core.solve import SolveOptions, builtin_structure, solve

#: 型号列表（内置 CIF 名）—— NCA 四档常见配比
MODELS = ("NCA_80_15_05", "NCA_82_13_05", "NCA_85_10_05", "NCA_90_05_05")

#: 合成谱真值：把 CIF 晶胞沿 a、c 各偏一点（远大于回收精度，防"原样返回"）
_A_MUL = 1.004
_C_MUL = 0.996

#: 真实实样路径（只读引用家族数据根）
_DATA = "D:/Project/XRD/XRData"
PRECURSOR_REAL = f"{_DATA}/csuHJW/Data043/Data043：Ni(OH)2.txt"
NCA_80_15_05_REAL = f"{_DATA}/csuHJW/txt/Data040：NCM：811_raw.txt"

#: Data043 随附 Maud 2.55 独立精修结果（`Data043：Ni(OH)2.par`）
MAUD_A = 3.130435
MAUD_C = 4.619966


def _y_deg_for(d_nm: float, wavelength: float = 1.5406) -> float:
    """晶粒尺寸 → TCH 洛伦兹尺寸项 Y（度），口径同 `test_roundtrip`。"""
    return math.degrees(1.0) * C.SCHERRER_K * wavelength / (d_nm * 10.0)


def _synth_model(model: str, *, seed: int = 7) -> tuple[Pattern, float, float]:
    """按型号 CIF 合成一份带 Poisson 噪声的谱，返回 (谱, a_gt, c_gt)。

    真值晶胞 = CIF 晶胞 × (1.004, 0.996) —— 起点故意错开，避免"输入等于
    输出"的空验证。
    """
    st0 = builtin_structure(model)
    a0, c0 = st0.cell[0], st0.cell[2]
    a_gt, c_gt = a0 * _A_MUL, c0 * _C_MUL
    st = CifStructure(f"{model}_gt", (a_gt, a_gt, c_gt, 90.0, 90.0, 120.0),
                      st0.atoms)
    st._variant = 771234          # noqa: SLF001 —— 缓存盐值，防 id 复用撞车
    prm = ProfileParams(U=0.03, V=-0.012, W=0.010, Y=_y_deg_for(120.0), X=0.0)
    tt = np.arange(12.0, 78.0, 0.03)
    prof = compute_phase_profile(st, tt, prm, 12.0, 78.0)
    bg = 40.0 + 30.0 * (tt - 12.0) / 66.0
    y = prof / prof.max() * 6000.0 + bg
    y = np.random.default_rng(seed).poisson(y).astype(np.float64)
    return Pattern(two_theta=tt, intensity=y, meta={}), a_gt, c_gt


@pytest.fixture(scope="module")
def model_fits() -> dict:
    """四型号各跑一次 stage1（峰形/晶胞/背景），module 级缓存。"""
    out = {}
    for model in MODELS:
        pat, a_gt, c_gt = _synth_model(model)
        res = solve(pat, builtin_structure(model), SolveOptions())
        out[model] = (res, a_gt, c_gt,
                      100.0 * math.sqrt(pat.intensity.size
                                        / pat.intensity.sum()))
    return out


# ---------------------------------------------------------------------------
# A · 型号泛化（合成谱）
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("model", MODELS)
def test_cell_a_recovered_per_model(model, model_fits) -> None:
    """逐型号：a 回收误差 ≤ 5 mÅ（实测 ≤ 0.11 mÅ，留 40 倍裕度）。"""
    res, a_gt, _c_gt, _floor = model_fits[model]
    assert abs(res.a - a_gt) <= 0.005, \
        f"{model}: a={res.a:.5f} vs 真值 {a_gt:.5f}"


@pytest.mark.parametrize("model", MODELS)
def test_cell_c_recovered_per_model(model, model_fits) -> None:
    """逐型号：c 回收误差 ≤ 10 mÅ（实测 ≤ 0.87 mÅ）。"""
    res, _a_gt, c_gt, _floor = model_fits[model]
    assert abs(res.c - c_gt) <= 0.010, \
        f"{model}: c={res.c:.5f} vs 真值 {c_gt:.5f}"


@pytest.mark.parametrize("model", MODELS)
def test_rwp_at_noise_floor_per_model(model, model_fits) -> None:
    """逐型号：Rwp 应贴着合成谱的 Poisson 统计下限（≈9.6 %）。

    **不能照搬真实样的 3 % 经验值** —— 弱背景合成谱（峰 6000 计数、
    背景 40–70 计数）的统计下限就在 9.6 % 量级。
    """
    res, _a_gt, _c_gt, floor = model_fits[model]
    assert res.rwp < 1.6 * floor, \
        f"{model}: Rwp={res.rwp:.3f} % vs 统计下限 {floor:.3f} %"


def test_all_models_converge_to_same_cell(model_fits) -> None:
    """四型号都能收敛，且**晶胞逐位相同** —— 这是 NCA 家族的物理事实。

    内置四份 NCA CIF 共用同一 LiNiO₂ 骨架晶胞（a=2.8645、c=14.161），差别
    只在 TM 位 Ni:Co 占位。XRD 分不开 80:15 / 82:13 / 85:10 / 90:05 ——
    Ni:Co:Al 比例须由 ICP/EDS 给出（见 `docs/STD_FREEZE.md` §1）。本断言
    把这个边界固化成事实：谁若给某档 CIF 单独改晶胞，这里立刻红。
    """
    assert len(model_fits) == len(MODELS)
    cells = {m: (model_fits[m][0].a, model_fits[m][0].c) for m in MODELS}
    a0, c0 = cells[MODELS[0]]
    for m, (a, c) in cells.items():
        assert abs(a - a0) < 1e-3, f"{m} 的 a={a:.5f} 与 {MODELS[0]} 的 {a0:.5f} 不同"
        assert abs(c - c0) < 5e-3, f"{m} 的 c={c:.5f} 与 {MODELS[0]} 的 {c0:.5f} 不同"


def test_model_intensity_differences_are_sub_percent() -> None:
    """四档配比的峰强差**确有其事，但小到 XRD 分不开** —— 边界声明的量化版。

    用同一组峰形参数正演四份谱（归一后取最大绝对差）：

    | 对比 | max\\|ΔI\\| |
    |---|---|
    | 80:15 vs 82:13 | 0.00028 |
    | 80:15 vs 85:10 | 0.00070 |
    | 80:15 vs 90:05 | 0.00140 |

    原因：TM 位平均散射因子 f̄ = 0.8×28 + 0.15×27 + 0.05×13 ≈ 27.1（80:15）
    对 0.9×28 + 0.05×27 + 0.05×13 ≈ 27.2（90:05）—— 差 0.4 %，落在任何
    实验室 XRD 的计数统计里。**所以型号必须由 ICP/EDS 定，不能由 XRD 定**。

    两条断言各管一头：差值必须 > 0（四份 CIF 不是彼此复制），且必须
    < 0.5 %（不能吹成"XRD 可区分型号"）。谁给某档 CIF 改散射因子/占位
    改出 1 % 以上的差异，这条会红。
    """
    tt = np.arange(15.0, 70.0, 0.02)
    prm = ProfileParams(U=0.02, V=-0.01, W=0.008)
    profs = {}
    for m in MODELS:
        y = compute_phase_profile(builtin_structure(m), tt, prm, 15.0, 70.0)
        profs[m] = y / max(y.max(), 1e-12)
    devs = {}
    for i, m1 in enumerate(MODELS):
        for m2 in MODELS[i + 1:]:
            dev = float(np.max(np.abs(profs[m1] - profs[m2])))
            devs[(m1, m2)] = dev
            assert dev > 1e-5, f"{m1} 与 {m2} 的谱逐点相同 —— CIF 未真正区分"
            assert dev < 5e-3, \
                f"{m1} 与 {m2} 的峰强差 {dev:.5f} 过大，与「XRD 不可分型号」矛盾"
    # 配比距离越远，差异越大（单调性自检，防"差异来自噪声"）
    assert devs[("NCA_80_15_05", "NCA_90_05_05")] > \
        devs[("NCA_80_15_05", "NCA_82_13_05")]


# ---------------------------------------------------------------------------
# B · 真实前驱体实样（Data043，(Ni,Co,Mn)(OH)₂，523 配比）
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def precursor_identify():
    """真实前驱体实样跑判定门。"""
    if not Path(PRECURSOR_REAL).exists():
        pytest.skip(f"实样缺失: {PRECURSOR_REAL}")
    return identify(load_pattern(PRECURSOR_REAL))


@pytest.fixture(scope="module")
def precursor_fit():
    """以三元氢氧化物 CIF 直接跑管线（绕过 GUI 判定门，验证精修本身）。

    管线对「结构里没有 Li 位点」的设计行为是**降级跳过 stage2**，
    不抛异常（见 `pipeline.run_pipeline` docstring）。
    """
    if not Path(PRECURSOR_REAL).exists():
        pytest.skip(f"实样缺失: {PRECURSOR_REAL}")
    st = builtin_structure("precursor_NCA_OH")
    return run_pipeline(load_pattern(PRECURSOR_REAL), st,
                        default_options_for(st))


def test_real_precursor_not_mistaken_for_noise(precursor_identify) -> None:
    """**T8.2 抓出的缺陷回归**：真实前驱体谱不得被判成「含大量噪声」。

    修前：寻峰器不看突出度，把 9 个真实峰切成 62 个噪声肩 → 峰数超
    `C_MAX_OBS_PEAKS` → 判定门在「噪声」分支直接返回 other，用户拿到的
    结论是"谱图含大量噪声"，与事实相反。
    """
    assert precursor_identify.metrics["n_peaks"] <= C.C_MAX_OBS_PEAKS, \
        f"峰数 {precursor_identify.metrics['n_peaks']} 超上限"
    assert not any("噪声" in r for r in precursor_identify.reasons), \
        f"仍被判成噪声: {precursor_identify.reasons}"


def test_real_precursor_routes_to_hydroxide(precursor_identify) -> None:
    """真实前驱体应分流到 `precursor_hydroxide`（**不是** nca、不是 other）。

    `ok=False` 是对的：前驱体是 P-3m1 氢氧化物，不能拿层状 R-3m 的 NCA
    正极 CIF 去精修。界面据此提示"这是前驱体"，而不是拒之门外。
    """
    assert precursor_identify.ok is False
    assert precursor_identify.kind == C.KIND_PRECURSOR_HYDROXIDE, \
        f"kind={precursor_identify.kind}｜{'；'.join(precursor_identify.reasons)}"


@pytest.mark.slow
def test_real_precursor_cell_matches_maud(precursor_fit) -> None:
    """**本卡最有价值的一项**：精修晶胞与 Maud 2.55 的独立结果一致。

    参考值取自实样随附的 `Data043：Ni(OH)2.par`（Maud 用 Popa 各向异性
    尺寸/应变模型、双线 Kα1/Kα2、背景多项式 + 手工节点，与我们的模型
    完全不同）—— 两条独立路线给出同一晶胞，才说明晶胞精修可信。
    """
    res = precursor_fit
    a, c = res.cell[0], res.cell[2]
    dev_a = abs(a - MAUD_A) / MAUD_A * 100.0
    dev_c = abs(c - MAUD_C) / MAUD_C * 100.0
    assert dev_a <= 0.3, f"a={a:.5f} vs Maud {MAUD_A}（偏差 {dev_a:.3f} %）"
    assert dev_c <= 0.3, f"c={c:.5f} vs Maud {MAUD_C}（偏差 {dev_c:.3f} %）"


@pytest.mark.slow
def test_real_precursor_skips_occupancy_without_breaking(precursor_fit) -> None:
    """无 Li 位点 → stage2 降级为 None，管线不中断（设计要求）。"""
    res = precursor_fit
    assert res.stage2 is None, "前驱体不应产出 Li 占位结果"
    assert res.occ_li is None and res.z_o is None
    assert res.structure is not None
    assert res.rwp > 0.0


@pytest.mark.slow
def test_real_precursor_rwp_consistent_with_noise(precursor_fit) -> None:
    """Rwp 应与该谱的统计下限同量级（实测 ≈2.1 × 下限）。

    这份谱很弱（最强峰 3389 计数、5–80° 共 3751 点），且前驱体是片状
    晶体、各向异性展宽明显 —— 单一 Caglioti 模型拟不到 Maud 那种
    Popa 各向异性模型的水平。断言按「统计下限的倍数」而非绝对百分比，
    避免把模型能力边界误判成缺陷。
    """
    pat = load_pattern(str(Path(PRECURSOR_REAL)))
    floor = 100.0 * math.sqrt(pat.intensity.size / pat.intensity.sum())
    assert precursor_fit.rwp < 2.6 * floor, \
        f"Rwp={precursor_fit.rwp:.3f} % vs 统计下限 {floor:.3f} %"


# ---------------------------------------------------------------------------
# C · 真实 NCA_80_15_05 锚点（T8.1 改动不得把已冻结的实样带偏）
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def nca811_fit():
    if not Path(NCA_80_15_05_REAL).exists():
        pytest.skip(f"实样缺失: {NCA_80_15_05_REAL}")
    st = builtin_structure("NCA_80_15_05")
    return run_pipeline(load_pattern(NCA_80_15_05_REAL), st, default_options_for(st))


@pytest.mark.slow
def test_real_811_rwp_anchor(nca811_fit) -> None:
    """同构锚点 Rwp ≤ 3.9 %（NCA 模型精修 NCM811 实样，实测 3.7452 %）。

    与 NCM 版的 3.1151 % 相比高 0.63pp —— 这不是退化，是**化学错配的
    预期代价**：实样 TM 位是 Mn，模型给的是 Al（散射因子 25 vs 13），
    8:1:1 与 8:1.5:0.5 的峰强分布本就不同。三条判据（C3 双峰比、Occ、D）
    仍各自成立，说明结构骨架对得上，只是模型化学不对。
    """
    assert nca811_fit.rwp <= 3.9, f"Rwp={nca811_fit.rwp:.4f} %"


@pytest.mark.slow
def test_real_811_microstructure_anchor(nca811_fit) -> None:
    """同构锚点：晶粒尺寸 105 nm 附近（实测 104.8，区间 95–116）。

    NCM 模型下是 155 nm；换 NCA 模型后样品宽化的一部分被结构性残差
    吸收到峰形参数里，Y 变小、D 随之变小。区间按 ±10 % 取。
    """
    d = nca811_fit.crystallite.D_nm
    assert d is not None and math.isfinite(d)
    assert 95.0 <= d <= 116.0, f"D={d:.1f} nm"


@pytest.mark.slow
def test_real_811_occupancy_anchor(nca811_fit) -> None:
    """同构锚点：Li 占位 0.987 ±0.5 %（NCM 模型实测 0.98689）。"""
    occ = nca811_fit.occ_li
    assert occ is not None
    assert 0.978 <= occ <= 0.995, f"occ_Li={occ:.5f}"


# ---------------------------------------------------------------------------
# D · 双候选择优的边界解让位（T8.2 抓出的缺陷）
# ---------------------------------------------------------------------------

def _cand(y: float, rwp: float):
    """构造 `_pick_candidate` 需要的最小候选（只有 rwp 与 Y 被读到）。"""
    from types import SimpleNamespace

    return SimpleNamespace(
        rwp=rwp, result=SimpleNamespace(params=SimpleNamespace(Y=y)))


def test_pick_candidate_prefers_interior_over_boundary() -> None:
    """边界解（Y=0）Rwp 仅略优时，应取内部解 —— 本卡实测场景。

    实测数据：合成 NCA_82_13_05 真值 Y=0.06043；谷扫描候选 Y=0.06058 /
    Rwp 12.294 %，普通候选 Y=0.0（D=inf）/ Rwp 12.190 %。按 Rwp 会挑中
    边界解，把 130 nm 的样品报成"超出检出上限"。
    """
    from nca.core.pipeline import _pick_candidate

    boundary = _cand(0.0, 12.190)
    interior = _cand(0.06058, 12.294)
    assert _pick_candidate([boundary, interior]) is interior
    assert _pick_candidate([interior, boundary]) is interior


def test_pick_candidate_keeps_boundary_when_interior_much_worse() -> None:
    """内部解明显更差（超过 5 % 容差）时仍取边界解 —— 规则不能盖过数据。"""
    from nca.core.pipeline import _pick_candidate

    boundary = _cand(0.0, 12.0)
    interior = _cand(0.06, 13.0)
    assert _pick_candidate([boundary, interior]) is boundary


def test_pick_candidate_plain_rwp_when_best_is_interior() -> None:
    """胜出者本身就是内部解时，行为与 T8.1 的纯 Rwp 取优逐位一致。"""
    from nca.core.pipeline import _pick_candidate

    good = _cand(0.05, 9.1)
    worse = _cand(0.06, 9.4)
    assert _pick_candidate([good, worse]) is good


def test_pick_candidate_both_boundary() -> None:
    """两个候选都停在边界 → 无内部解可让，仍按 Rwp 取优。"""
    from nca.core.pipeline import _pick_candidate

    a = _cand(0.0, 10.2)
    b = _cand(0.0, 10.1)
    assert _pick_candidate([a, b]) is b
