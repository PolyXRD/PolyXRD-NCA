"""物相判定门（T5.1）。

**目的**：导入谱后先判断「这是不是 NCA（或其前驱体）」，不是就停下并提示。
本模块**不跑精修**，只做峰位层面的快速结构类型判定 —— 必须能在 GUI 里
秒级返回。

三重判据（STD_FREEZE §1，阈值唯一出处 `nca/core/constants.py`）
--------------------------------------------------------------
| 判据 | 指标 | 阈值 |
|---|---|---|
| C1 | 前 12 强峰与 R-3m 理论峰的 2θ 均方根偏差 | ≤ 0.10° |
| C2 | (006)/(012) 与 (108)/(110) 双峰可分辨 | Δ2θ ≥ 0.10° |
| C3 | I(003)/I(104) | ∈ [0.5, 3.0] |

**C1 的「理论峰」是标定后的（R-3m 晶胞快速标定）** —— 这是本卡的关键口径，
STD_FREEZE v3 已补记：

1. CIF 自带的晶胞是文献值，与实测差 0.2–0.3%（实测 a=2.87185 vs CIF 2.86450），
   而 2θ 偏差随 tanθ 放大：18.8° 处 −0.04°，44.5° 处 −0.10°，77° 处 −0.23°。
   直接拿 CIF 晶胞比，C1 的 RMS 实测 0.164°，**恒超 0.10° 阈值**，判定门形同虚设。
2. 判定门先做一次**快速晶胞标定**（匹配 → LM 拟合 a,c → 重匹配，3 轮），
   再用标定后的理论峰跑 C1。实测 RMS 0.164° → **0.0138°**，标定得
   a=2.87158 / c=14.19411，与 40 秒全谱精修的锚点 (2.87185 / 14.19502)
   相差 0.01%。
3. 标定检验的是**结构类型**（层状 R-3m 的峰位序列特征），不检验晶胞绝对值 ——
   这正是 STD_FREEZE §1 边界声明要表达的：NCA 四档配比 / NCM / LiNiO₂ 的晶胞
   几乎相同，实验室 XRD 分不开，Ni:Co:Al 比例须由 ICP/EDS 确定。

**分流**：C1–C3 全过 → `nca`；否则在候选池里找**结构可信**的匹配，按其空间群
归类为 `precursor_hydroxide` / `precursor_carbonate`；都不匹配 → `other`。

**结构可信度门槛（三道，缺一不可）**
----------------------------------
晶胞快速标定只有 (a, c) 两个自由度 —— **什么谱都能拟合出一个像样的峰位
RMS**。只按 RMS 挑最佳候选会翻车（T5.1 实测缺陷，已修）：

| 样品 | 旧口径（只看 RMS） | 修正后 |
|---|---|---|
| 岩盐 NiO | 判成 `precursor_carbonate`（RMS 0.069°） | 无可信候选 → `other` |
| 尖晶石 LiMn₂O₄ | 判成 `nca` | 无可信候选 → `other` |
| 金红石 TiO₂ | 判成 `nca`（RMS 0.079°） | 无可信候选 → `other` |
| 纯噪声谱 | 判成 `precursor_hydroxide` | 无可信候选 → `other` |

三道门槛：
1. **峰位**：命中 ≥ `C1_MIN_MATCH` 且 RMS ≤ `C1_RMS_MAX`；
2. **晶胞**：标定后相对原 CIF 晶胞的偏差 ≤ `CALIB_CELL_TOL_PCT`（3%）——
   标定能把尖晶石的 a=8.25 Å 拉到 2.80 Å，这道护栏专治此类；
3. **强度指纹**：整谱 Pearson 相关 ≥ `C_FINGERPRINT_MIN`（0.85）——
   真结构 0.92–0.98，非结构 −0.02–0.55，判别裕度充足。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

from nca.core import constants as C
from nca.core import fingerprint
from nca.core.calibrate import (
    MatchPair,
    calibrate,
    calibrate_pool,
    cell_drift_pct,
    match,
    rms,
)
from nca.core.criteria import (
    HKL_003,
    HKL_006,
    HKL_012,
    HKL_104,
    HKL_108,
    HKL_110,
    c1,
    c2,
    c3,
)
from nca.core.dtypes import Pattern
from nca.core.peaks import (
    WAVELENGTH_KA1,
    ObservedPeak,
    find_peaks,
    prepare,
    reflections,
)
from nca.core.profile.cif import CifStructure

__all__ = [
    "HKL_003",
    "HKL_006",
    "HKL_012",
    "HKL_104",
    "HKL_108",
    "HKL_110",
    "IdentifyResult",
    "MatchPair",
    "category_of",
    "default_pool",
    "identify",
]

#: 兼容别名 —— 保 `tests/`、`tools/` 与历史内部调用点的旧名不破
_c1 = c1
_c2 = c2
_c3 = c3
_calibrate = calibrate
_calibrate_pool = calibrate_pool
_cell_drift_pct = cell_drift_pct
_match = match
_rms = rms

#: 空间群国际表号 → 类别
#:
#: NCA 版在 NCM 三分类之外补两类（都是 NCA 场景里必须能被单独指认的相）：
#: - 92 / 14 / 63 = 含铝副相：γ-LiAlO₂(P4₁2₁2) / 三水铝石(P2₁/n) / 勃姆石(Cmcm)；
#: - 225 = 岩盐相 NiO —— NCA 表面重构与欠氧化前驱体的标志。
#: 注意 166 同时被 γ-LiAlO₂ 的高压 R-3m 变体命中，故层状放行还要过
#: `NCA_FAMILY_A_RANGE` 的 a 区间护栏（见 constants 注释）。
_SG_CATEGORY = {
    166: C.KIND_NCA,                     # R -3 m：层状氧化物
    164: C.KIND_PRECURSOR_HYDROXIDE,     # P -3 m 1：类 β-Ni(OH)2 氢氧化物
    167: C.KIND_PRECURSOR_CARBONATE,     # R -3 c：方解石型碳酸盐
    92: C.KIND_AL_PHASE,                 # P 41 21 2：γ-LiAlO2
    14: C.KIND_AL_PHASE,                 # P 1 21/n 1：三水铝石 Al(OH)3
    63: C.KIND_AL_PHASE,                 # C m c m：勃姆石 AlOOH
    225: C.KIND_ROCKSALT,                # F m -3 m：NiO 岩盐
}

#: 判定门内置参考结构（STD_FREEZE 的「R-3m 理论峰」出处）
REFERENCE_NAME = "NCA_80_15_05"


@dataclass(frozen=True)
class IdentifyResult:
    """判定门结论。

    - `ok`：True = **可以进入 NCA 精修流程**（仅当 kind == nca 且 C1–C3 全过）。
    - `kind`：结构类型（nca / precursor_hydroxide / precursor_carbonate / other），
      即使 `ok` 为 False 也给出，用于界面措辞（「这是前驱体，请走前驱体流程」）。
    - `reasons`：逐条中文理由（界面直接显示）。
    - `metrics`：量化明细（可进结果卡 / 日志）。
    """

    ok: bool
    kind: str
    reasons: tuple[str, ...] = ()
    metrics: dict = field(default_factory=dict)

    def summary(self) -> str:
        head = "通过" if self.ok else "未通过"
        return f"[{head}] {self.kind}｜" + "；".join(self.reasons)


# ---------------------------------------------------------------------------
# 候选池
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def default_pool() -> tuple[CifStructure, ...]:
    """内置 CIF 全池。

    分四组（顺序即优先级：同名匹配时先命中的进 `_reference`）：
    1. **层状 NCA**（4 档配比 + 端元 LiNiO₂）—— 参考结构与族属判据的对象；
    2. **对照层状相**：LiCoO₂（a≈2.814，族区间外）、脱锂层状相；
    3. **前驱体**：NCA 氢氧化物/碳酸盐 + β-Ni(OH)₂ / NiCO₃ 骨架；
    4. **NCA 副相**：γ-LiAlO₂、三水铝石、勃姆石、岩盐 NiO、残碱 Li₂CO₃。

    第 4 组进池的意义不是"能被判成 NCA"，而是让判定门**说得出它是什么**：
    铝偏析、岩盐表面相、残碱都有明确的处置建议，笼统报 other 对产线无用。
    """
    from nca.core.solve import builtin_structure

    names = (
        # 层状 NCA
        "NCA_80_15_05", "NCA_82_13_05", "NCA_85_10_05", "NCA_90_05_05",
        "LiNiO2",
        # 对照层状相
        "LiCoO2", "layered_delith",
        # 前驱体
        "precursor_NCA_OH", "precursor_NiOH2",
        "precursor_NCA_CO3", "precursor_NiCO3",
        # NCA 副相
        "LiAlO2", "AlOH3_gibbsite", "AlOOH_boehmite",
    )
    out = []
    for n in names:
        try:
            out.append(builtin_structure(n))
        except (FileNotFoundError, ValueError):
            continue
    # 杂质相（NiO 岩盐 / Li₂CO₃ 残碱）单独加载：它们在 impurity 子目录，
    # 由 multiphase 的加载器统一入口读，避免两处各自拼路径。
    try:
        from nca.core.multiphase import _load_impurity_cif

        for stem in ("NiO", "Li2CO3"):
            try:
                out.append(_load_impurity_cif(stem))
            except (FileNotFoundError, ValueError):
                continue
    except ImportError:
        pass
    return tuple(out)


def category_of(struct: CifStructure) -> str:
    """按空间群归类（CIF 未带空间群时归 other）。"""
    from nca.core.symmetry import lookup

    sg = getattr(struct, "spacegroup", None)
    if not sg:
        return C.KIND_OTHER
    try:
        info = lookup(sg)
    except KeyError:
        return C.KIND_OTHER
    if info.number is None:
        return C.KIND_OTHER
    return _SG_CATEGORY.get(int(info.number), C.KIND_OTHER)


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def identify(pat: Pattern, cif_pool=None, *, anode: str = "Cu",
             ka2_ratio: float = 0.5) -> IdentifyResult:
    """判定门主流程。

    参数
    ----
    pat : 原始谱（可含背景与 Kα2 —— 判定门内部自行预处理）
    cif_pool : 候选结构列表；None = 内置 14 份

    返回语义的约定
    --------------
    - `ok=True` **当且仅当** kind == nca 且 C1–C3 全过 → 可以进入 NCA 精修。
    - `kind` 是**结构类型**：nca / precursor_hydroxide / precursor_carbonate /
      other。`ok=False` 时 kind 仍会给出（例如 LiCoO₂ 是层状 R-3m，kind=nca
      但 ok=False），界面据此措辞（"是层状氧化物但晶胞不是 NCA，是否仍按 NCA
      精修？"），而不是一律拒之门外。
    """
    lam = WAVELENGTH_KA1
    net = prepare(pat, ka2_ratio=ka2_ratio, anode=anode)
    obs = find_peaks(net, rel_min=C.PEAK_REL_MIN)
    metrics: dict = {"n_peaks": len(obs),
                     "peak_rel_min": C.PEAK_REL_MIN}

    if len(obs) < C.C_MIN_OBS_PEAKS:
        return IdentifyResult(False, C.KIND_OTHER,
                              (f"可用峰仅 {len(obs)} 个（不足 "
                               f"{C.C_MIN_OBS_PEAKS}），谱图可能不是 XRD 粉末谱。",),
                              metrics)
    if len(obs) > C.C_MAX_OBS_PEAKS:
        return IdentifyResult(False, C.KIND_OTHER,
                              (f"寻得 {len(obs)} 个峰，超过上限 "
                               f"{C.C_MAX_OBS_PEAKS} —— 谱图含大量噪声，"
                               f"不是可用的 XRD 粉末谱。",),
                              metrics)

    pool = tuple(cif_pool) if cif_pool else default_pool()
    ref = _reference(pool)
    if ref is None:
        return IdentifyResult(False, C.KIND_OTHER,
                              ("候选池中没有可用的 R-3m 参考结构。",), metrics)

    # --- 快速路径：只用内置参考结构标定，判据全过就收工（绝大多数情况）---
    st_ref = ref
    cell, pairs = _calibrate(ref, obs, lam)
    c1_ok, c1 = _c1(obs, _refl_for_cell(st_ref, cell, lam))
    c2_ok, c2 = _c2(obs, cell, lam) if c1_ok else (False, None)
    c3_ok, c3 = _c3(obs, cell, lam) if c1_ok else (False, None)

    # --- 完整路径：快速路径未全过 → 用池中最佳标定重试 ---
    pool_cal = []
    if not (c1_ok and c2_ok and c3_ok):
        pool_cal = _calibrate_pool(pool, obs, lam)
        if pool_cal:
            st_ref, cell, pairs, _rms0 = pool_cal[0]
            c1_ok, c1 = _c1(obs, _refl_for_cell(st_ref, cell, lam))
            c2_ok, c2 = _c2(obs, cell, lam) if c1_ok else (False, None)
            c3_ok, c3 = _c3(obs, cell, lam) if c1_ok else (False, None)

    metrics["c1"] = c1
    metrics["c2"] = c2 if c2 else {"pass": None, "evaluated": False,
                                   "reason": "C1 未通过，未评估"}
    metrics["c3"] = c3 if c3 else {"pass": None, "evaluated": False,
                                   "reason": "C1 未通过，未评估"}
    metrics["calibrated_cell"] = {"a": round(cell[0], 5), "c": round(cell[2], 5)}
    metrics["calib_hits"] = len(pairs)
    metrics["calib_rms"] = round(_rms(pairs), 5) if pairs else None
    metrics["reference"] = st_ref.name

    reasons: list[str] = []
    if c1_ok:
        reasons.append(f"C1 通过：{c1['hits']} 个主峰峰位 RMS "
                       f"{c1['rms_deg']:.4f}° ≤ {C.C1_RMS_MAX}°。")
    else:
        reasons.append(f"C1 未通过：命中 {c1['hits']} 个主峰，RMS "
                       f"{c1['rms_deg'] if c1['rms_deg'] is not None else '—'}°"
                       f"（阈值 {C.C1_RMS_MAX}°，至少需命中 "
                       f"{C.C1_MIN_MATCH} 个）。")

    c2_ok = c3_ok = False
    if c1_ok:
        c2_ok, c2 = _c2(obs, cell, lam)
        c3_ok, c3 = _c3(obs, cell, lam)
        metrics["c2"], metrics["c3"] = c2, c3
        reasons.append("C2 通过：(006)/(012) 与 (108)/(110) 双峰均可分辨。"
                       if c2_ok else "C2 未通过：层状特征双峰未能分辨。")
        reasons.append(f"C3 通过：I(003)/I(104) = {c3['ratio']:.3f}。"
                       if c3_ok else
                       f"C3 未通过：I(003)/I(104) = "
                       f"{c3['ratio'] if c3['ratio'] is not None else '—'}"
                       f"（应落在 {C.C3_RATIO_MIN}–{C.C3_RATIO_MAX}）。")
    else:
        metrics["c2"] = {"pass": None, "evaluated": False,
                          "reason": "C1 未通过，未评估"}
        metrics["c3"] = {"pass": None, "evaluated": False,
                          "reason": "C1 未通过，未评估"}

    if c1_ok and c2_ok and c3_ok:
        a_lo, a_hi = C.NCA_FAMILY_A_RANGE
        a_val = cell[0]
        reasons.append(f"标定晶胞 a={a_val:.4f} Å、c={cell[2]:.4f} Å。")
        if not (a_lo <= a_val <= a_hi):
            reasons.append(f"晶胞族属判据未通过：标定 a={a_val:.4f} Å 不在 NCA "
                           f"族区间 {a_lo}–{a_hi} Å —— 疑似 LiCoO₂ "
                           f"(a≈2.814 Å)、γ-LiAlO₂ 高压 R-3m 变体 (a≈2.800 Å) "
                           f"等其他层状氧化物，或脱锂体系；"
                           f"该谱不按 NCA 成品放行（层状结构类型仍予保留）。")
            metrics["nca_family_a"] = {"value": round(a_val, 5),
                                       "range": [a_lo, a_hi], "pass": False}
            return IdentifyResult(False, C.KIND_NCA, tuple(reasons), metrics)
        metrics["nca_family_a"] = {"value": round(a_val, 5),
                                   "range": [a_lo, a_hi], "pass": True}
        reasons.append(f"晶胞族属判据通过：a 落在 NCA 族 {a_lo}–{a_hi} Å。")
        reasons.append("提示：NCA / NCM / LiNiO₂ 同为层状 R-3m 且晶胞几乎相同，"
                       "实验室 XRD 无法区分 Ni:Co:Al 与 Ni:Co:Mn 的比例 —— "
                       "本软件判定的是**结构类型**，元素配比须由 ICP/EDS 确定。")
        return IdentifyResult(True, C.KIND_NCA, tuple(reasons), metrics)

    reasons.append("判据未全部通过，转入结构类型分流。")
    trusted = _best_match(pool_cal, obs, lam, net, pat.two_theta)
    metrics["trusted_matches"] = trusted
    if not trusted:
        reasons.append("候选池中无任何结构能以可信方式匹配该谱"
                       f"（峰位 RMS ≤ {C.C1_RMS_MAX}°、晶胞偏差 ≤ "
                       f"{C.CALIB_CELL_TOL_PCT}%、指纹相关 ≥ "
                       f"{C.C_FINGERPRINT_MIN}）。")
        return IdentifyResult(False, C.KIND_OTHER, tuple(reasons), metrics)

    best = trusted[0]
    kind = best["category"]
    reasons.append(f"最佳可信匹配「{best['name']}」：峰位 RMS "
                   f"{best['rms_deg']:.4f}°、命中 {best['hits']} 峰、"
                   f"指纹相关 {best['corr']:.3f}、标定 a={best['a']:.4f} Å。")
    if kind == C.KIND_PRECURSOR_HYDROXIDE:
        reasons.append("→ 氢氧化物前驱体（P-3m1），不是层状氧化物成品；"
                       "请改用前驱体流程（晶胞 + 晶粒尺寸，不含 Li 占位）。")
    elif kind == C.KIND_PRECURSOR_CARBONATE:
        reasons.append("→ 碳酸盐前驱体（R-3c），不是层状氧化物成品。")
    elif kind == C.KIND_AL_PHASE:
        reasons.append("→ 含铝副相（γ-LiAlO₂ / Al(OH)₃ / AlOOH）—— "
                       "铝在共沉淀或烧结中偏析成独立相，不是 NCA 层状主相。"
                       "建议核查铝源加入方式与前驱体洗涤/干燥工艺。")
    elif kind == C.KIND_ROCKSALT:
        reasons.append("→ 岩盐相 NiO（Fm-3m）—— NCA 表面重构或前驱体欠氧化的标志，"
                       "不是层状 R-3m 主相。")
    elif kind == C.KIND_NCA:
        reasons.append("→ 层状 R-3m 锂过渡金属氧化物，但与标准 NCA 参考的"
                       "晶胞/峰位不完全吻合（可能是 LiCoO₂、脱锂或富锂体系）。")
    else:
        reasons.append("→ 判定为非 NCA 样品。")
    return IdentifyResult(False, kind, tuple(reasons), metrics)


def _reference(pool: tuple[CifStructure, ...]) -> CifStructure | None:
    """判定门的 R-3m 参考结构：优先内置 NCA_80_15_05，否则池中任一 R-3m。"""
    for st in pool:
        if st.name == REFERENCE_NAME:
            return st
    for st in pool:
        if category_of(st) == C.KIND_NCA:
            return st
    return None


def _refl_for_cell(ref: CifStructure, cell, lam: float):
    """按标定后晶胞重算反射表（强度随晶胞微变，hkl 不变）。"""
    scaled = CifStructure(ref.name, tuple(cell), ref.atoms)
    return reflections(scaled, lam, C.IDENTIFY_TMIN, C.IDENTIFY_TMAX)


def _theoretical_clusters(ref, tol: float) -> int:
    """窗口内**可分辨**的理论峰簇数（覆盖率口径的分母）。

    为什么不能直接用 `len(reflections(...))`：反射表里 Kα1 与 Kα2 是两条独立
    记录（岩盐 NiO 的 (111) 记成 43.38° 与 43.40° 两行），直接计数会把同一
    条衍射线算成两个反射 —— 实测 NiO 得 7 条「理论反射」，而实际只剩 4 条
    可分辨线，覆盖率被压到 4/7=0.57，门槛失效。

    这里按 `C1_MATCH_TOL_DEG` 聚类（与标定/匹配用的容差同口径）后计数，
    保证「分母」与「分子（命中数）」是同一个分辨尺度下的量。
    """
    pos = sorted(float(r.two_theta) for r in ref)
    if not pos:
        return 0
    n = 1
    last = pos[0]
    for p in pos[1:]:
        if p - last > tol:
            n += 1
            last = p
    return n


def _best_match(pool_cal, obs: list[ObservedPeak], lam: float, net: Pattern,
                grid: np.ndarray) -> list[dict]:
    """在**已标定**的候选里挑出结构可信的匹配（按峰位 RMS 升序）。

    三道门槛（缺一不可）：
    1. 峰位：`命中 ≥ C1_MIN_MATCH` **或** `反射覆盖率 ≥ C1_MIN_COVERAGE`
       （后者专为高对称少峰相而设，见 `constants.C1_MIN_COVERAGE`），
       且 RMS ≤ `C1_RMS_MAX`；
    2. 晶胞：标定后相对原 CIF 晶胞的偏差 ≤ `CALIB_CELL_TOL_PCT`
       —— 标定只有 2 个自由度，没有这道护栏，什么谱都能被"拟合"成候选；
    3. 强度指纹：整谱 Pearson 相关 ≥ `C_FINGERPRINT_MIN`
       —— 只看峰位会把"强行标定"的假匹配当真（T5.1 实测缺陷）。
       指纹本身在 `nca/core/fingerprint.py`（峰宽自适应，T8.2）。
    """
    out: list[dict] = []
    for st, cell, pairs, pair_rms in pool_cal:
        if pair_rms > C.C1_RMS_MAX:
            continue
        n_theo = _theoretical_clusters(_refl_for_cell(st, cell, lam),
                                       C.C1_MATCH_TOL_DEG)
        coverage = (len(pairs) / n_theo) if n_theo else 0.0
        if len(pairs) < C.C1_MIN_MATCH and coverage < C.C1_MIN_COVERAGE:
            continue
        da, dc = _cell_drift_pct(cell, st.cell)
        if max(da, dc) > C.CALIB_CELL_TOL_PCT:
            continue
        corr = fingerprint.corr(st, cell, net, grid)
        if corr < C.C_FINGERPRINT_MIN:
            continue
        out.append({"name": st.name, "category": category_of(st),
                    "rms_deg": round(pair_rms, 5), "hits": len(pairs),
                    "n_theoretical": n_theo, "coverage": round(coverage, 3),
                    "corr": round(corr, 4),
                    "a": round(cell[0], 5), "c": round(cell[2], 5),
                    "cell_drift_pct": [round(da, 2), round(dc, 2)]})
    out.sort(key=lambda r: r["rms_deg"])
    return out
