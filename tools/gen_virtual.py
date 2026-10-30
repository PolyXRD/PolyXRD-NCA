"""生成含尺寸/应变宽化的虚拟粉末谱（T2.2 验收工具）。

用法
----
```bash
# 生成 D = 50 / 100 / 200 nm 的虚拟 NCA_80_15_05 谱（默认输出到 tests/data/virtual/）
D:/Project/XRD/python310/python.exe tools/gen_virtual.py --size 50 100 200

# 同时施加微观应变，并指定输出目录
D:/Project/XRD/python310/python.exe tools/gen_virtual.py --size 100 --strain 0.001 \
    --out tests/data/virtual
```

输出为**两列文本**（2θ / 强度），与 `nca.core.patterns.load_pattern` 可直接读的
格式一致 —— 生成的谱可以直接喂给 `solve()` 做往返验证。

**口径**（与 `docs/STD_FREEZE.md` §3 一致）：
- 尺寸项 `Y = (180/π)·Kλ/D`（度，常数）
- 应变项 `X = (180/π)·4ε`（度，常数）
- 峰形用 TCH 合成（`nca.core.profile.size_strain`），η 逐峰由 TCH 给出
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from nca.core.crystallite import D_DETECT_MAX_NM
from nca.core.profile.forward import ProfileParams, compute_phase_profile
from nca.core.profile.size_strain import SCHERRER_K, tch_combine, tch_eta
from nca.core.profile.struct import compute_peaks
from nca.core.solve import builtin_structure

#: Cu Kα1 波长（Å）
LAMBDA_KA1 = 1.540598
#: Cu Kα2 波长（Å）与相对强度比 I(Kα2)/I(Kα1)
#:
#: ★ 合成谱**必须**含 Kα2（NCA 版修正）。判定门 `peaks.prepare()` 按 Cu 靶
#: 无条件做 Rachinger 剥离；若合成谱只有 Kα1，剥离算法仍会减掉一个位移副本，
#: 在每个强峰的高角侧留下等间距伪肩峰（实测岩盐 NiO：真峰 4 条被切成 14 条，
#: 肩峰间距 0.22°→0.44° 随 2θ 增大，正是 Kα2 位移随 tanθ 放大的特征）。
#: 真实 Cu 靶实验室谱本就带 Kα2，合成谱跟上才是正确建模。
LAMBDA_KA2 = 1.544426
KA2_RATIO = 0.5
#: 默认仪器高斯宽化（Caglioti）
DEFAULT_U, DEFAULT_V, DEFAULT_W = 0.03, -0.01, 0.008


def y_from_size_nm(d_nm: float) -> float:
    """尺寸 → TCH 洛伦兹常数 Y（度）。"""
    d_ang = d_nm * 10.0
    return (180.0 / math.pi) * SCHERRER_K * LAMBDA_KA1 / d_ang


def x_from_strain(eps: float) -> float:
    """应变 → TCH 洛伦兹常数 X（度）。"""
    return (180.0 / math.pi) * 4.0 * eps


def synth_tch(struct, grid: np.ndarray, d_nm: float, eps: float,
              u: float, v: float, w: float) -> np.ndarray:
    """按 TCH 峰形正演一条合成谱（含 Cu Kα1 + Kα2 双线）。"""
    y_ang = y_from_size_nm(d_nm)
    x_ang = x_from_strain(eps)
    peaks, _ = compute_peaks(
        struct, wavelength=LAMBDA_KA1,
        two_min=float(grid[0]), two_max=float(grid[-1]),
        wavelengths=[(LAMBDA_KA1, 1.0), (LAMBDA_KA2, KA2_RATIO)])
    prof = np.zeros_like(grid, dtype=np.float64)
    for pos, height, _d in peaks:
        if height <= 0:
            continue
        theta = math.radians(pos / 2.0)
        t = math.tan(theta)
        fg = math.sqrt(max(u * t * t + v * t + w, 1e-8))
        fl = x_ang * t + y_ang / math.cos(theta)
        fwhm = max(tch_combine(fg, fl), 1e-4)
        eta = tch_eta(fg, fl)
        # 伪 Voigt：单位面积归一化后乘以积分强度
        area = fwhm * (1.06447 + 0.50633 * eta) / 2.0
        pv = (eta / (1.0 + eta)) * np.exp(-4.0 * math.log(2.0)
                                         * ((grid - pos) / fwhm) ** 2)
        lor = (1.0 / (1.0 + eta)) / (1.0 + 4.0 * ((grid - pos) / fwhm) ** 2)
        prof += height / max(area, 1e-12) * (pv + lor)
    return prof


def write_pattern(path: Path, grid: np.ndarray, y: np.ndarray,
                  header: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# {header}"]
    lines += [f"{t:.4f} {i:.4f}" for t, i in zip(grid, y)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="生成含尺寸/应变的虚拟粉末谱")
    ap.add_argument("--size", type=float, nargs="+", required=True,
                    help="晶粒尺寸 D（nm），可给多个")
    ap.add_argument("--strain", type=float, default=0.0,
                    help="微观应变 ε（无量纲），默认 0")
    ap.add_argument("--structure", default="NCA_80_15_05", help="结构名（内置 CIF）")
    ap.add_argument("--out", default="tests/data/virtual", help="输出目录")
    ap.add_argument("--tmin", type=float, default=10.0)
    ap.add_argument("--tmax", type=float, default=90.0)
    ap.add_argument("--step", type=float, default=0.02)
    ap.add_argument("--seed", type=int, default=20260701)
    args = ap.parse_args(argv)

    struct = builtin_structure(args.structure)
    grid = np.arange(args.tmin, args.tmax, args.step)
    out_dir = ROOT / args.out if not Path(args.out).is_absolute() else Path(args.out)
    rng = np.random.default_rng(args.seed)

    for d_nm in args.size:
        if d_nm <= 0:
            raise SystemExit(f"尺寸必须为正: {d_nm}")
        note = ""
        if d_nm > D_DETECT_MAX_NM:
            note = f" [超出检出上限 {D_DETECT_MAX_NM:g} nm]"
        prof = synth_tch(struct, grid, d_nm, args.strain,
                         DEFAULT_U, DEFAULT_V, DEFAULT_W)
        prof = prof / max(prof.max(), 1e-12) * 10000.0
        y = prof + 200.0 + rng.normal(0.0, np.sqrt(np.maximum(prof, 1.0)))
        name = f"virtual_{args.structure}_D{int(d_nm)}"
        if args.strain:
            name += f"_eps{args.strain:g}"
        header = (f"PolyXRD-NCA 虚拟谱 | 结构 {args.structure} | "
                  f"D={d_nm:g} nm | ε={args.strain:g} | "
                  f"X={x_from_strain(args.strain):.5f}° "
                  f"Y={y_from_size_nm(d_nm):.5f}°{note}")
        path = out_dir / f"{name}.txt"
        write_pattern(path, grid, y, header)
        print(f"写出 {path}  (D={d_nm:g} nm, ε={args.strain:g}){note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
