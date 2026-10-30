"""构建判定门样本库（T5.2 · NCA 版）。

从家族只读数据根挑选**真实**谱作负样本/锚点，缺失的类别用**合成**谱补齐，
并生成 `tests/data/manifest.json`（每个样本的来源、期望结论）。

用法
----
```bash
D:/Project/XRD/python310/python.exe tools/build_sample_lib.py
```

**设计原则**
- **真实优先**：尽量用家族实测谱（含仪器宽化、织构、背景、非晶等真实世界的
  所有麻烦），合成谱只用于补齐没有实测数据的类别。
- **数据根只读**：只复制到项目内 `tests/data/`，不写回 `D:/Project/XRD`。
- **期望值写进 manifest**：判定门的行为由 manifest 断言，不靠人记。

★ **NCA 的关键数据现状（必须如实记录）**
--------------------------------------
家族数据根里**没有 NCA 实测谱**。`Data040：NCM：811` 是与 NCA 同构的层状
R-3m 实样（同为 #166，a/c 相差 <0.2%），本库把它作为**同构锚点**使用，
在 manifest 的 note 里写明真实来源与"XRD 不能区分 Ni:Co:Mn 与 Ni:Co:Al"。
NCA 专属配比的正样本全部是**合成谱**（由内置 CIF 正演），这是当前唯一的
可复现手段，也是 v0.1.0 的已知局限（见 docs/VALIDATION.md §局限）。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "tools") not in sys.path:
    sys.path.insert(0, str(ROOT / "tools"))

from gen_virtual import (DEFAULT_U, DEFAULT_V, DEFAULT_W,  # noqa: E402
                         synth_tch, x_from_strain, y_from_size_nm)
from nca.core.solve import builtin_structure  # noqa: E402

DATA_ROOT = Path("D:/Project/XRD")
LIB = ROOT / "tests" / "data"
SEED = 20260701

#: 合成谱的默认 2θ 网格（贴近真实实验室谱）
GRID = np.arange(10.0, 80.0, 0.02)


def _synth_text(struct: CifStructure, d_nm: float, eps: float = 0.0,
                seed: int = SEED) -> str:
    """合成一条两列文本谱（含背景 + 计数噪声）。"""
    prof = synth_tch(struct, GRID, d_nm, eps, DEFAULT_U, DEFAULT_V, DEFAULT_W)
    prof = prof / max(prof.max(), 1e-12) * 10000.0
    rng = np.random.default_rng(seed)
    y = prof + 200.0 + rng.normal(0.0, np.sqrt(np.maximum(prof, 1.0)))
    lines = [f"# PolyXRD-NCA 合成谱 | {struct.name} | D={d_nm:g} nm | "
             f"ε={eps:g} | X={x_from_strain(eps):.5f}° "
             f"Y={y_from_size_nm(d_nm):.5f}° | seed={seed}"]
    lines += [f"{t:.4f} {i:.4f}" for t, i in zip(GRID, y)]
    return "\n".join(lines) + "\n"


def _copy_real(src_rel: str, dst: Path) -> None:
    """从数据根复制真实谱（只读数据根，只写项目内）。

    统一转成**两列文本**（`#` 头 + 两列数值），不保留原扩展名：
    原文件有 .RAW/.TXT/.mdi 等多种格式，扩展名一变读谱器的格式嗅探就失效。
    """
    from nca.core.patterns import load_pattern

    src = DATA_ROOT / src_rel
    if not src.exists():
        raise FileNotFoundError(f"数据根缺少: {src}")
    pat = load_pattern(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# 来源: {src_rel}",
             f"# 区间: {pat.two_theta[0]:.4f}-{pat.two_theta[-1]:.4f}° "
             f"({pat.two_theta.size} 点)"]
    lines += [f"{t:.4f} {i:.4f}"
              for t, i in zip(pat.two_theta, pat.intensity)]
    dst.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# 样本定义
# ---------------------------------------------------------------------------
#: (类别, 文件名, 数据根相对路径, 说明, expect_ok, expect_kind)
#: expect_kind 取值：nca / precursor_hydroxide / precursor_carbonate /
#:                   al_phase / rocksalt / other
REAL_SAMPLES = [
    # --- 层状 R-3m 真实锚点（NCM811 实样；与 NCA 同构） ---
    ("positive/nca", "layered_anchor_real.txt",
     "XRData/csuHJW/txt/Data040：NCM：811_raw.txt",
     "真实层状 R-3m 实样（NCM811）。**家族数据根无 NCA 实测谱**，本样本作为"
     "同构锚点：NCA 与 NCM 同为 #166、晶胞相差 <0.2%，实验室 XRD 不可区分",
     True, "nca"),
    # --- 尖晶石 ---
    ("negative/spinel", "alcoo_650c_real.txt",
     "XRData/csuHJW/data005：AlCoO-650C.txt", "AlCoO₂ 尖晶石", False, "other"),
    ("negative/spinel", "alcoo_750c_real.txt",
     "XRData/csuHJW/data004：AlCoO-750C.txt", "AlCoO₂ 尖晶石", False, "other"),
    ("negative/spinel", "magnetite_iucr.txt",
     "XRData/IUCr/MAGNETIT.RAW", "Fe₃O₄ 尖晶石", False, "other"),
    # --- 橄榄石（LFP / FePO₄）---
    ("negative/lfp", "lifepo4_smz.txt",
     "XRData/smz geshi/LiFePO4#01##20250520-100220_100.txt",
     "LiFePO₄ 橄榄石 Pnma", False, "other"),
    ("negative/lfp", "fepo4_smz.txt",
     "XRData/smz geshi/FePO4#01##20250520-132913_100.txt",
     "FePO₄ 橄榄石型", False, "other"),
    # --- 岩盐型（Fm-3m）---
    ("negative/rocksalt", "fluorite_iucr.txt",
     "XRData/IUCr/FLUORITE.RAW", "CaF₂，Fm-3m 萤石型", False, "other"),
    # --- 前驱体 ---
    ("negative/precursor", "ni_oh2_csu.txt",
     "XRData/csuHJW/Data043/Data043：Ni(OH)2.txt",
     "真实三元氢氧化物前驱体（Ni:Co:Mn=5:3:2，随附 Maud 精修 a=3.130435/"
     "c=4.619966）。★ 本工具能正确指认为 precursor_hydroxide —— 因为 NCA "
     "候选池含 NCA 专属氢氧化物（precursor_NCA_OH，同一 P-3m1 骨架）；"
     "该样本 `ok=False`（前驱体不是可精修的 NCA 成品），kind 判对是诊断价值所在",
     False, "precursor_hydroxide"),
    # --- 其他非层状 / 混合相 ---
    ("negative/other", "anatase_csu.txt",
     "XRData/csuHJW/txt/Data033：单相锐钛矿.txt", "单相锐钛矿", False, "other"),
    ("negative/other", "y2o3_csu.txt",
     "XRData/csuHJW/txt/Data038：Y2O3.txt", "Y₂O₃ 方铍石型（立方）",
     False, "other"),
    ("negative/other", "zno_caco3_csu.txt",
     "XRData/csuHJW/Data041/Data041：ZnO+CaCO3.txt",
     "ZnO + CaCO₃（CaCO₃ 为 R-3c）", False, "other"),
    ("negative/other", "limno2_si_csu.txt",
     "XRData/csuHJW/txt/Data030：LiMnO2+Si.txt", "正交 LiMnO₂ + Si",
     False, "other"),
    ("negative/other", "zrb_zrb2_csu.txt",
     "XRData/csuHJW/txt/Data031：ZrB-ZrB2.txt", "六方硼化物", False, "other"),
]

def _nio_structure():
    """岩盐 NiO 结构（读真实 COD 1010093，不手搓原子表）。

    ★ 为什么必须读 CIF 而不是手写：早先版本手写
    `[("Ni",0,0,0),("O",.5,.5,.5)]` + 立方胞 a=4.18 —— 那是**原胞描述塞进
    惯用立方胞**，只给 2 个原子。结构因子把 NaCl 型本该消光的 (100)/(110)
    等全算了出来（实测正演出 11 条线而不是 4 条），合成谱根本不像岩盐，
    判定门自然认不出。NCM 版把这个 bug 掩盖了（对该样本期望就是 other），
    换成"期望 rocksalt"后才暴露。
    """
    from nca.core.multiphase import _load_impurity_cif

    return _load_impurity_cif("NiO")


def build() -> list[dict]:
    manifest: list[dict] = []

    for folder, name, source, note, exp_ok, exp_kind in REAL_SAMPLES:
        dst = LIB / folder / name
        _copy_real(source, dst)
        manifest.append({
            "file": str((LIB / folder / name).relative_to(LIB)).replace("\\", "/"),
            "category": folder.split("/")[1],
            "origin": "real",
            "source": source,
            "note": note,
            "expect_ok": exp_ok,
            "expect_kind": exp_kind,
        })

    # --- 合成谱 ---
    # (类别, 文件名, 结构, D_nm, ε, 说明, expect_ok, expect_kind)
    synth_specs = [
        # 正样本：四档 NCA 配比 + LiNiO₂ 端元
        ("positive/nca", "nca_80_15_05_syn.txt", "NCA_80_15_05", 150.0, 0.0,
         "主配比 NCA 80/15/5（派生自 COD 1520789）合成", True, "nca"),
        ("positive/nca", "nca_82_13_05_syn.txt", "NCA_82_13_05", 200.0, 0.0,
         "NCA 82/13/5 合成", True, "nca"),
        ("positive/nca", "nca_85_10_05_syn.txt", "NCA_85_10_05", 150.0, 0.0,
         "NCA 85/10/5 合成", True, "nca"),
        ("positive/nca", "nca_90_05_05_syn.txt", "NCA_90_05_05", 120.0, 0.0,
         "高镍 NCA 90/5/5 合成", True, "nca"),
        ("positive/nca", "linio2_syn.txt", "LiNiO2", 120.0, 0.0005,
         "LiNiO₂ 端元（COD 1533770）合成，含 0.05% 应变", True, "nca"),
        # 层状但不是 NCA 族：LiCoO₂（a=2.814 在族区间外）
        ("negative/lco", "lco_d120nm_syn.txt", "LiCoO2", 120.0, 0.0,
         "LiCoO₂（COD 1533825）合成 —— 层状 R-3m 但 a 不在 NCA 族区间",
         False, "nca"),
        ("negative/lco", "lco_d200nm_syn.txt", "LiCoO2", 200.0, 0.0,
         "LiCoO₂ 合成（另一晶粒尺寸）", False, "nca"),
        # 岩盐
        ("negative/rocksalt", "nio_d80nm_syn.txt", None, 80.0, 0.0,
         "岩盐 NiO (Fm-3m, a=4.177) 合成 —— NCA 表面重构标志相",
         False, "rocksalt"),
        # 前驱体
        ("negative/precursor", "nico3_d150nm_syn.txt", "precursor_NiCO3",
         150.0, 0.0, "碳酸盐前驱体骨架（COD 1548823，R-3c）合成",
         False, "precursor_carbonate"),
        ("negative/precursor", "nih2_d150nm_syn.txt", "precursor_NiOH2",
         150.0, 0.0, "氢氧化物前驱体骨架（COD 1548811，P-3m1）合成",
         False, "precursor_hydroxide"),
        ("negative/precursor", "nca_oh_d150nm_syn.txt", "precursor_NCA_OH",
         150.0, 0.0, "NCA 氢氧化物前驱体（80:15:5）合成", False,
         "precursor_hydroxide"),
        # 含铝副相（NCA 专属）
        ("negative/al_phase", "lialo2_d120nm_syn.txt", "LiAlO2", 120.0, 0.0,
         "γ-LiAlO₂（COD 4002658，四方 P4₁2₁2）合成 —— 铝偏析副相",
         False, "al_phase"),
    ]
    for (folder, name, struct_name, d_nm, eps, note,
         exp_ok, exp_kind) in synth_specs:
        struct = _nio_structure() if struct_name is None else builtin_structure(struct_name)
        dst = LIB / folder / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(_synth_text(struct, d_nm, eps), encoding="utf-8")
        manifest.append({
            "file": str((LIB / folder / name).relative_to(LIB)).replace("\\", "/"),
            "category": folder.split("/")[1],
            "origin": "synthetic",
            "source": note,
            "note": note,
            "expect_ok": exp_ok,
            "expect_kind": exp_kind,
        })

    path = LIB / "manifest.json"
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2),
                    encoding="utf-8")
    return manifest


def main() -> int:
    manifest = build()
    print(f"样本库已构建：{len(manifest)} 个样本 -> {LIB}")
    cats: dict[str, int] = {}
    for m in manifest:
        cats[m["category"]] = cats.get(m["category"], 0) + 1
    for k in sorted(cats):
        print(f"  {k:<12} {cats[k]}")
    pos = sum(1 for m in manifest if m["expect_ok"])
    print(f"  正样本 {pos} / 负样本 {len(manifest) - pos}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
