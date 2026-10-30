"""第二相槽位与多相定量测试（T5.3 验收）。

**核心判据：Hill–Howard 定量往返**
合成 NCA_80_15_05 + Li₂CO₃ 混合谱（按 `w_p ∝ S_p·Z_p·M_p·V_p` 自洽给 scale）→
多相精修 → Hill–Howard 换算 → 恢复 90/10（容差 ±2 wt%）。

合成口径的推导（不可改，否则定量往返不成立）：
精修的线性层用 lstsq 解出各相 `S_p`（profile_p 已知，无噪声时精确等于
生成时给的值）。`S_p ∝ 衍射晶胞数 ∝ w_p/(Z_p·M_p·V_p)`，故
`w_p ∝ S_p·Z_p·M_p·V_p` —— Hill–Howard 成立。

**★酸检（RutAna 继承）**：ZMV 关系必须正确。TiO₂ 锐钛型(Z=4)/
金红石型(Z=2)：(ZMV)_锐/(ZMV)_金 = (4×136.3)/(2×62.42) ≈ 4.368。
"""
from __future__ import annotations

import numpy as np
import pytest

from nca.core.multiphase import (
    MAX_IMPURITIES,
    ROLE_IMPURITY,
    ROLE_MAIN,
    PhaseSlot,
    builtin_impurities,
    make_plan,
    quantify,
    suggest_impurities,
    unexplained_peaks,
)
from nca.core.profile.forward import ProfileParams, compute_phase_profile
from nca.core.refine.nonlinear import RefineConfig, refine
from nca.core.refine.qpa import (
    format_formula,
    formula_mass,
    infer_z,
    qpa_weights,
    z_consistency,
)
from nca.core.solve import builtin_structure

GRID = np.arange(10.0, 80.0, 0.02)
PRM = ProfileParams(zero=0.0, U=0.03, V=-0.01, W=0.008, eta=0.6,
                    anode="Cu", ka_ratio=0.5)
SEED = 20260701

#: 主相化学式**不手填** —— 由 CIF 自带的 `_chemical_formula_sum` 决定。
#: 内置 NCA_80_15_05 是 COD 1520789（Arai et al., J. Power Sources 90 (2000) 76），
#: 实为锂富集镍欠缺型 Li1.03Ni0.77Co0.1Mn0.1O2，与教科书配方
#: LiNi0.8Co0.1Mn0.1O2 不同（TM 比 Ni:Co:Mn = 0.77:0.1:0.1 ≈ 8:1:1，
#: 但 Li/TM = 1.06 > 1）。手填教科书配方会让 Z 反算失败、摩尔质量算错。
MAIN_FORMULA = None


def _main_slot() -> PhaseSlot:
    st = builtin_structure("NCA_80_15_05")
    return PhaseSlot(name="NCA_80_15_05", structure=st, formula=MAIN_FORMULA,
                     role=ROLE_MAIN)


def _imp(name: str) -> PhaseSlot:
    for s in builtin_impurities():
        if s.name == name:
            return s
    raise KeyError(name)


def _synth_mixed(slots_weights, seed: int = SEED):
    """按 Hill–Howard 自洽合成多相谱。返回 (Pattern, 设计质量分数列表)。"""
    from nca.core.dtypes import Pattern

    rows = []
    for slot, w in slots_weights:
        s = slot.scale_for_weight(w)
        rows.append(s * compute_phase_profile(slot.structure, GRID, PRM))
    total = np.sum(rows, axis=0)
    total = total / max(total.max(), 1e-12) * 10000.0
    rng = np.random.default_rng(seed)
    y = total + 200.0 + rng.normal(0.0, np.sqrt(np.maximum(total, 1.0)))
    return Pattern(two_theta=GRID.copy(), intensity=y)


# ---------------------------------------------------------------------------
# QPA 基础（★酸检）
# ---------------------------------------------------------------------------

def test_qpa_acid_check_zmv_relation() -> None:
    """★酸检：(ZMV)_锐/(ZMV)_金 = (4×136.3)/(2×62.42) ≈ 4.368。"""
    z_a, v_a, m_a = 4, 1.0, 136.3      # 锐钛型：Z=4, V=1（仅比值用）
    z_r, v_r, m_r = 2, 1.0, 62.42      # 金红石型：Z=2
    ratio = (z_a * v_a * m_a) / (z_r * v_r * m_r)
    assert ratio == pytest.approx(4.368, rel=1e-3)
    # 同 scale 时锐钛相质量分数远高于金红石（不是 50%）
    ws = qpa_weights([("锐", 1.0, 4, 1.0, 136.3), ("金", 1.0, 2, 1.0, 62.42)])
    assert ws[0].weight == pytest.approx(0.8137, abs=1e-3)


def test_formula_mass() -> None:
    assert formula_mass("TiO2") == pytest.approx(79.865, abs=1e-3)
    assert formula_mass("Li2CO3") == pytest.approx(73.889, abs=1e-2)
    assert formula_mass("NiO") == pytest.approx(74.692, abs=1e-2)
    assert formula_mass("C") == pytest.approx(12.011, abs=1e-3)


def test_infer_z_for_project_phases() -> None:
    """各相的 Z（晶胞内化学式单元数）必须正确。"""
    # 不给化学式 —— 由 CIF 自带配方反算
    assert infer_z(builtin_structure("NCA_80_15_05")) == 3
    for name, z in (("Li2CO3", 4), ("NiO", 4), ("LiAlO2", 4)):
        slot = _imp(name)
        assert slot.zmv()[0] == z, f"{name} 的 Z 应为 {z}"


def test_formula_comes_from_cif_not_handwritten() -> None:
    """★锁定内置 CIF 的真实配方（不可用别人的配方替代）。

    内置 NCA_80_15_05 由 COD 1520789（LiNiO₂ 骨架）改占位派生，
    `_chemical_formula_sum` = `Al0.05 Co0.15 Li1 Ni0.8 O2`，
    与 NCA 常见工业配比 80:15:5 一致，Z=3 反算自洽。

    反例用 **NCM 教科书配方** `LiNi0.8Co0.1Mn0.1O2`：本 CIF 里根本没有 Mn，
    Co 占位 0.45 也对不上 0.1×3=0.3，反算必然抛错（而不是悄悄算出一个错的 Z）。
    """
    st = builtin_structure("NCA_80_15_05")
    assert st.formula is not None, "NCA_80_15_05 CIF 应带 _chemical_formula_sum"
    assert format_formula(st.formula) == "LiNi0.8Co0.15Al0.05O2"
    assert infer_z(st) == 3
    # 别家配方必须被拒绝（而不是悄悄算出一个错的 Z）
    with pytest.raises(ValueError):
        infer_z(st, "LiNi0.8Co0.1Mn0.1O2")


def test_z_consistency_flags_declared_mismatch() -> None:
    """CIF 自述 Z 与反算 Z 不一致时要能报出来（内置前驱体即属此列）。"""
    st = builtin_structure("NCA_80_15_05")
    info = z_consistency(st, infer_z(st))
    assert info["z_computed"] == 3
    assert info["consistent"] is True
    assert info["note"] is None
    # 自述值不符时给提示
    bad = z_consistency(st, 5)
    assert bad["consistent"] is False
    assert "不符" in (bad["note"] or "")


def test_qpa_weights_rejects_invalid() -> None:
    with pytest.raises(ValueError):
        qpa_weights([])
    with pytest.raises(ValueError):
        qpa_weights([("x", -1.0, 1, 1.0, 1.0)])


# ---------------------------------------------------------------------------
# 槽位
# ---------------------------------------------------------------------------

def test_builtin_impurities_present() -> None:
    names = {s.name for s in builtin_impurities()}
    assert {"Li2CO3", "NiO", "LiAlO2"} <= names
    for s in builtin_impurities():
        assert s.role == ROLE_IMPURITY
        z, vol, mass = s.zmv()
        assert z >= 1 and vol > 0 and mass > 0


def test_plan_respects_impurity_limit() -> None:
    main = _main_slot()
    # 默认路径：取内置候选的前 MAX 个
    plan = make_plan(main.structure, MAIN_FORMULA)
    assert len(plan.impurities) <= MAX_IMPURITIES
    assert plan.slots[0].role == ROLE_MAIN
    # 显式给 ≤ MAX 个 —— 通过
    plan2 = make_plan(main.structure, MAIN_FORMULA,
                      impurities=list(builtin_impurities())[:MAX_IMPURITIES])
    assert len(plan2.impurities) == MAX_IMPURITIES
    # 显式给 > MAX 个 —— **必须报错，不得静默截断**
    with pytest.raises(ValueError):
        make_plan(main.structure, MAIN_FORMULA,
                  impurities=[main] * (MAX_IMPURITIES + 1))


def test_plan_phases_align_with_slots() -> None:
    main = _main_slot()
    plan = make_plan(main.structure, MAIN_FORMULA,
                     impurities=[_imp("Li2CO3")])
    assert len(plan.phases()) == len(plan.enabled_slots) == 2


# ---------------------------------------------------------------------------
# 定量往返（核心判据）
# ---------------------------------------------------------------------------

def test_quantify_roundtrip_two_phase() -> None:
    """NCA_80_15_05 90% + Li₂CO₃ 10% → 精修 → 定量 → 恢复 90/10。"""
    main = _main_slot()
    imp = _imp("Li2CO3")
    pat = _synth_mixed([(main, 0.90), (imp, 0.10)])

    cfg = RefineConfig(refine_zero=True, refine_uvw=True, refine_eta=True,
                       n_bg=8, tmin=12.0, tmax=79.0, bg_mode="cheb",
                       refine_asym=False, refine_bfactor=False,
                       refine_cell=False, max_iter=25)
    res = refine([main.structure, imp.structure], pat, PRM, cfg)
    assert res.success
    assert len(res.scales) == 2

    ws = quantify([main, imp], list(res.scales))
    assert ws[0].weight + ws[1].weight == pytest.approx(1.0, abs=1e-9)
    assert ws[0].weight == pytest.approx(0.90, abs=0.02), \
        f"主相质量分数 {ws[0].weight:.4f} 偏离设计值 0.90"
    assert ws[1].weight == pytest.approx(0.10, abs=0.02), \
        f"杂质相质量分数 {ws[1].weight:.4f} 偏离设计值 0.10"


def test_quantify_roundtrip_three_phase() -> None:
    """主相 + 2 个杂质相：80/12/8。"""
    main = _main_slot()
    imp1 = _imp("Li2CO3")
    imp2 = _imp("NiO")
    pat = _synth_mixed([(main, 0.80), (imp1, 0.12), (imp2, 0.08)])

    cfg = RefineConfig(refine_zero=True, refine_uvw=True, refine_eta=True,
                       n_bg=8, tmin=12.0, tmax=79.0, bg_mode="cheb",
                       refine_asym=False, refine_bfactor=False,
                       refine_cell=False, max_iter=25)
    res = refine([main.structure, imp1.structure, imp2.structure], pat,
                 PRM, cfg)
    ws = quantify([main, imp1, imp2], list(res.scales))
    assert ws[0].weight == pytest.approx(0.80, abs=0.03)
    assert ws[1].weight == pytest.approx(0.12, abs=0.03)
    assert ws[2].weight == pytest.approx(0.08, abs=0.03)


def test_scale_for_weight_is_inverse_of_hill_howard() -> None:
    """`scale_for_weight` 与 Hill–Howard 严格互逆（纯代数自洽）。"""
    main = _main_slot()
    imp = _imp("NiO")
    s_main = main.scale_for_weight(0.75)
    s_imp = imp.scale_for_weight(0.25)
    ws = qpa_weights([("m", s_main, *main.zmv()), ("i", s_imp, *imp.zmv())])
    assert ws[0].weight == pytest.approx(0.75, abs=1e-9)
    assert ws[1].weight == pytest.approx(0.25, abs=1e-9)


def test_quantify_rejects_length_mismatch() -> None:
    main = _main_slot()
    with pytest.raises(ValueError):
        quantify([main, _imp("NiO")], [1.0])


# ---------------------------------------------------------------------------
# 残差诊断
# ---------------------------------------------------------------------------

def test_unexplained_peaks_find_impurity_signal() -> None:
    """单相精修 NCA_80_15_05（含 Li₂CO₃）的残差里应能找到未解释峰。"""
    main = _main_slot()
    imp = _imp("Li2CO3")
    pat = _synth_mixed([(main, 0.90), (imp, 0.10)])

    cfg = RefineConfig(refine_zero=True, refine_uvw=True, refine_eta=True,
                       n_bg=8, tmin=12.0, tmax=79.0, bg_mode="cheb",
                       refine_asym=False, refine_bfactor=False,
                       refine_cell=False, max_iter=25)
    res = refine([main.structure], pat, PRM, cfg)
    peaks = unexplained_peaks(pat.two_theta, pat.intensity, res.ycalc)
    assert peaks, "单相模型未留下任何未解释峰"
    assert all(p.height > 0 for p in peaks)


def test_suggest_impurities_identifies_li2co3() -> None:
    """残差诊断应把 Li₂CO₃ 列为首要建议。"""
    main = _main_slot()
    imp = _imp("Li2CO3")
    pat = _synth_mixed([(main, 0.90), (imp, 0.10)])

    cfg = RefineConfig(refine_zero=True, refine_uvw=True, refine_eta=True,
                       n_bg=8, tmin=12.0, tmax=79.0, bg_mode="cheb",
                       refine_asym=False, refine_bfactor=False,
                       refine_cell=False, max_iter=25)
    res = refine([main.structure], pat, PRM, cfg)
    peaks = unexplained_peaks(pat.two_theta, pat.intensity, res.ycalc)
    sugg = suggest_impurities(peaks)
    assert sugg, "未给出任何杂质相建议"
    assert sugg[0].name == "Li2CO3", f"首位建议是 {sugg[0].name}"
    d = sugg[0].as_dict()
    assert d["hits"] >= 1


def test_no_false_impurity_on_pure_sample() -> None:
    """★负控制：纯单相样品不应报出任何杂质相。

    残差诊断必须"宁可不报，不可乱报" —— 否则界面会给用户一份假的杂质清单。
    """
    main = _main_slot()
    pat = _synth_mixed([(main, 1.00)])

    cfg = RefineConfig(refine_zero=True, refine_uvw=True, refine_eta=True,
                       n_bg=8, tmin=12.0, tmax=79.0, bg_mode="cheb",
                       refine_asym=False, refine_bfactor=False,
                       refine_cell=False, max_iter=25)
    res = refine([main.structure], pat, PRM, cfg)
    peaks = unexplained_peaks(pat.two_theta, pat.intensity, res.ycalc)
    sugg = suggest_impurities(peaks)
    assert not sugg, \
        f"纯单相样品误报杂质: {[(s.name, s.hits) for s in sugg]}"
