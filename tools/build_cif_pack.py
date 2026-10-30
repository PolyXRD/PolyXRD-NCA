"""内置 CIF 打包器（T0.3 · NCA 版）。

从 `D:/Project/XRD/sync-COD/COD/` 取 COD 原始 CIF 落盘，并按 NCA 名义配比
派生出层状 NCA（80/15/5、82/13/5、85/10/5、90/5/5）与三元氢氧化物 / 碳酸盐
前驱体，最后生成 INDEX.json。

NCA 与 NCM 的关键差异（本文件落地要点）
----------------------------------------
1. **层状相 TM 位换 Al**：Al³⁺ 取代 Mn 进入 3b 位，散射因子相差一个量级
   （Z: Al 13 vs Mn 25），|F|² 与各峰积分强度比值随之改变 —— 这也是本工具
   与 PolyXRD-NCA 家族其他成员共用同一套精修内核、却需要独立 CIF 库的原因。
2. **Al 相关副相必须入库**：γ-LiAlO₂（P4₁2₁2）、三水铝石 Al(OH)₃（P2₁/n）、
   勃姆石 AlOOH（Cmcm）。NCA 的 Al 在共沉淀/烧结中易偏析成独立相，这些相
   的峰若不入候选池，判定门会把它们当"未解释峰"。
3. **γ-LiAlO₂ 不是层状 R-3m**（P4₁2₁2，a=5.1965）—— 但 COD 里另有一条
   高压 R-3m 的 LiAlO₂（9009880，a=2.8003），与 NCA 同空间群，属于判定门
   必须能排除的假阳性（靠 `constants.NCA_FAMILY_A_RANGE` 的 a 区间护栏）。

设计要点
--------
1. **原样复制**的 CIF 只做书写规范化（去 esd 括号、原子 loop 重排 7 列），
   数值一个不动，可溯源、可校验 sha256。
2. **派生**的 CIF 只改 `_chemical_formula_sum` 与 TM 位占位，保留原始对称操作
   （parse_cif 依赖 `_symmetry_equiv_pos_as_xyz` 展开，缺了会少算原子）。
3. 派生文件头写入 `# PolyXRD derived from COD <id>`，明确非实测结构。
4. 全部走幂等写入：重复运行结果一致。

用法::

    D:/Project/XRD/python310/python.exe tools/build_cif_pack.py
    D:/Project/XRD/python310/python.exe tools/build_cif_pack.py --check   # 只校验不写
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COD_ROOT = Path("D:/Project/XRD/sync-COD/COD")
OUT_DIR = ROOT / "nca" / "data" / "cif"
IMP_DIR = OUT_DIR / "impurity"
INDEX = OUT_DIR / "INDEX.json"

#: 原样复制清单（主目录）：目标名 -> (COD 编号, 说明)
COPY_LIST: dict[str, tuple[int, str]] = {
    "LiNiO2.cif": (1533770, "LiNiO₂ 端元（真实测定，NCA 的 Ni 极限参考）"),
    "LiCoO2.cif": (1533825, "LiCoO₂（真实测定；判定门负样本对照，a≈2.814 在 NCA 族区间外）"),
    "layered_delith.cif": (
        4002445,
        "层状 R-3m 脱锂态参考（真实测定 Li₀.₆₃；用于充电态/脱锂峰位对照，"
        "非 NCA 专用配比）",
    ),
    "precursor_NiOH2.cif": (1548811, "β-Ni(OH)₂ 氢氧化物前驱体骨架 P-3m1（真实测定）"),
    "precursor_NiCO3.cif": (1548823, "NiCO₃ 碳酸盐前驱体骨架 R-3c（真实测定）"),
    "LiAlO2.cif": (
        4002658,
        "γ-LiAlO₂ P4₁2₁2（真实测定）—— NCA 表面/残铝副相，不是层状 R-3m",
    ),
    "AlOH3_gibbsite.cif": (1200016, "三水铝石 Al(OH)₃ P2₁/n（真实测定）—— 铝源残留相"),
    "AlOOH_boehmite.cif": (9012247, "勃姆石 AlOOH Cmcm（真实测定）—— 铝源脱水中间相"),
}

#: 杂质候选（impurity 子目录）：目标名 -> (COD 编号, 说明)
IMPURITY_LIST: dict[str, tuple[int, str]] = {
    "NiO.cif": (1010093, "岩盐 NiO Fm-3m（真实测定）—— NCA 表面重构/前驱体欠氧化相"),
    "Li2CO3.cif": (9008283, "残碱 Li₂CO₃ C2/c（真实测定）"),
}

#: 派生清单：目标名 -> (模板 COD, 新化学式, TM 位占位, 说明)
#:
#: 占位口径：Li 层继承模板的 0.98 Li + 0.02 Ni（2% Li/Ni 反位），
#: **TM 位显式给 Li 0.02** 把反位配平回去 —— 这样总 Li = 1.00、
#: 总 Ni = 目标值，化学式与占位严格自洽（Hill–Howard 定量要用）。
DERIVE_LIST: dict[str, tuple[int, str, dict[str, float], str]] = {
    "NCA_80_15_05.cif": (
        1520789, "Al0.05 Co0.15 Li1 Ni0.8 O2",
        {"Ni": 0.78, "Co": 0.15, "Al": 0.05, "Li": 0.02},
        "NCA 80/15/5 主配比派生（TM 位 78:15:5 + 2% Li/Ni 反位 -> 总 Ni 0.80、Li 1.00）",
    ),
    "NCA_82_13_05.cif": (
        1520789, "Al0.05 Co0.13 Li1 Ni0.82 O2",
        {"Ni": 0.80, "Co": 0.13, "Al": 0.05, "Li": 0.02},
        "NCA 82/13/5 派生（TM 位 80:13:5 + 2% 反位 -> 总 Ni 0.82、Li 1.00）",
    ),
    "NCA_85_10_05.cif": (
        1520789, "Al0.05 Co0.1 Li1 Ni0.85 O2",
        {"Ni": 0.83, "Co": 0.10, "Al": 0.05, "Li": 0.02},
        "NCA 85/10/5 派生（TM 位 83:10:5 + 2% 反位 -> 总 Ni 0.85、Li 1.00）",
    ),
    "NCA_90_05_05.cif": (
        1520789, "Al0.05 Co0.05 Li1 Ni0.9 O2",
        {"Ni": 0.88, "Co": 0.05, "Al": 0.05, "Li": 0.02},
        "NCA 90/5/5 高镍派生（TM 位 88:5:5 + 2% 反位 -> 总 Ni 0.90、Li 1.00）",
    ),
    "precursor_NCA_OH.cif": (
        1548811, "Al0.05 Co0.15 H2 Ni0.8 O2",
        {"Ni": 0.80, "Co": 0.15, "Al": 0.05},
        "NCA 氢氧化物前驱体（80:15:5，β-Ni(OH)₂ 骨架 P-3m1）",
    ),
    "precursor_NCA_CO3.cif": (
        1548823, "C Al0.05 Co0.15 Ni0.8 O3",
        {"Ni": 0.80, "Co": 0.15, "Al": 0.05},
        "NCA 碳酸盐前驱体（80:15:5，方解石型骨架 R-3c）",
    ),
}

#: 元素 -> 离子标识（写进 _atom_site_type_symbol）
ION = {"Li": "Li+1", "Ni": "Ni+3", "Co": "Co+3", "Mn": "Mn+3", "Al": "Al+3",
       "O": "O-2", "H": "H1+", "C": "C4+"}

#: 过渡金属元素集合（用于定位 TM 位，Li 必须排除）
TM_ELEMS = {"Ni", "Co", "Mn", "Fe", "Al", "Cr", "V", "Ti", "Cu", "Zn"}


def _elem(type_symbol: str) -> str:
    """从 'Ni+3' / 'O-2' / 'Al+3' 里取元素符号。"""
    out = []
    for ch in type_symbol:
        if ch.isalpha():
            out.append(ch)
        else:
            break
    return "".join(out)


def cod_path(cod_id: int) -> Path:
    """COD 三级布局路径：1/00/00/1000003.cif"""
    s = str(cod_id)
    rel = Path(s[0]) / s[1:3] / s[3:5] / f"{s}.cif"
    return COD_ROOT / rel


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def split_atom_block(lines: list[str]) -> tuple[int, int, int]:
    """定位原子位 loop 的数据区起止 [start, end)，返回 (header_idx, start, end)。

    返回的 header_idx 指向 `_atom_site_label` 所在行。
    """
    h = None
    for i, ln in enumerate(lines):
        if ln.strip() == "_atom_site_label":
            h = i
            break
    if h is None:
        raise ValueError("未找到 _atom_site_label")
    j = h
    while j < len(lines) and lines[j].strip().startswith("_atom_site"):
        j += 1
    start = j
    k = start
    while k < len(lines):
        s = lines[k].strip()
        if not s or s.startswith(("loop_", "_", "#", "data_")):
            break
        k += 1
    return h, start, k


def _strip_esd(tok: str) -> str:
    """去掉 CIF 的 esd 括号写法：'3.1268(1)' -> '3.1268'，'0' 原样。"""
    t = tok.strip().strip("'\"")
    if "(" in t and t.endswith(")"):
        head = t.split("(")[0]
        if head not in ("", "-", "+"):
            return head
    return t


_CELL_TAGS = (
    "_cell_length_a", "_cell_length_b", "_cell_length_c",
    "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma",
    "_cell_volume", "_cell_formula_units_Z",
)
_ATOM_COLS = ("label", "type_symbol", "fract_x", "fract_y", "fract_z",
              "occupancy", "U_iso_or_equiv")


def normalize_cif_text(text: str) -> str:
    """把 CIF 规范化为家族解析器可直接吃的形式（只改书写形式，不改数值）。

    1. 去掉晶胞/原子数值的 esd 括号（`3.1268(1)` -> `3.1268`）；
    2. 把 `_atom_site_*` 循环重排为 7 列标准集，丢弃解析器不读的列。

    原因：家族解析器用 `float()` 直读，遇到 esd 括号会返回 None ——
    晶胞退化成 1.0、带 esd 的原子被**静默丢弃**。COD 原始件普遍带 esd，
    必须在这里一次性规范化。
    """
    lines = text.splitlines()
    fixed: list[str] = []
    for ln in lines:
        s = ln.strip()
        hit = next((t for t in _CELL_TAGS if s.startswith(t)), None)
        if hit and len(s.split(None, 1)) > 1:
            tag, val = s.split(None, 1)
            fixed.append(f"{tag:<34}{_strip_esd(val)}")
        else:
            fixed.append(ln)

    try:
        h, start, end = split_atom_block(fixed)
    except ValueError:
        return "\n".join(fixed) + "\n"
    tags = [ln.strip() for ln in fixed[h:start]]
    rows = [_split_row(ln) for ln in fixed[start:end] if ln.strip()]

    def idx(key: str) -> int | None:
        t = f"_atom_site_{key}"
        return tags.index(t) if t in tags else None

    i_lbl, i_sym = idx("label"), idx("type_symbol")
    i_x, i_y, i_z = idx("fract_x"), idx("fract_y"), idx("fract_z")
    i_occ, i_u = idx("occupancy"), idx("U_iso_or_equiv")
    if None in (i_x, i_y, i_z):
        return "\n".join(fixed) + "\n"

    def cell(row: list[str], k: int | None, default: str) -> str:
        if k is None or k >= len(row):
            return default
        return _strip_esd(row[k])

    new_rows = []
    for r in rows:
        label = r[i_lbl] if i_lbl is not None and i_lbl < len(r) else "X"
        sym = r[i_sym] if i_sym is not None and i_sym < len(r) else label
        new_rows.append(" ".join([
            label, sym,
            cell(r, i_x, "0"), cell(r, i_y, "0"), cell(r, i_z, "0"),
            cell(r, i_occ, "1"), cell(r, i_u, "0"),
        ]))

    head = [f"_atom_site_{c}" for c in _ATOM_COLS]
    rebuilt = fixed[:h] + head + new_rows + fixed[end:]
    return "\n".join(rebuilt) + "\n"


def _split_row(line: str) -> list[str]:
    """切分 CIF 数据行（引号内视为一个 token）。"""
    s = line.strip()
    if "'" in s or '"' in s:
        try:
            import shlex

            return shlex.split(s)
        except ValueError:
            return s.split()
    return s.split()


def read_cell_and_formula(lines: list[str]) -> tuple[str, float, float]:
    """读取 (空间群, a, c)；缺失给 0。"""
    sg, a, c = "", 0.0, 0.0
    for ln in lines:
        s = ln.strip()
        if s.startswith("_symmetry_space_group_name_H-M"):
            sg = s.split(None, 1)[1].strip().strip("'\"") if len(s.split(None, 1)) > 1 else ""
        elif s.startswith("_cell_length_a"):
            a = _num(s)
        elif s.startswith("_cell_length_c"):
            c = _num(s)
    return sg, a, c


def _num(tag_line: str) -> float:
    raw = tag_line.split(None, 1)[1] if len(tag_line.split(None, 1)) > 1 else "0"
    raw = raw.strip().split("(")[0]
    try:
        return float(raw)
    except ValueError:
        return 0.0


def derive(template: Path, formula: str, tm_occ: dict[str, float], note: str) -> str:
    """以模板为骨架，替换 TM 位占位，返回新 CIF 文本。"""
    lines = template.read_text(encoding="utf-8", errors="replace").splitlines()
    h, start, end = split_atom_block(lines)
    tags = [ln.strip() for ln in lines[h:start]]
    rows = [ln.split() for ln in lines[start:end]]
    if not rows:
        raise ValueError(f"模板原子行为空: {template}")

    def col(name: str) -> int:
        return tags.index(f"_atom_site_{name}")

    i_lbl, i_sym = col("label"), col("type_symbol")
    i_x, i_y, i_z = col("fract_x"), col("fract_y"), col("fract_z")
    i_occ = col("occupancy")

    # 找 TM 位：占位最大的**过渡金属**行（Li 层占位也高，绝不能当锚点）
    tm_rows = [r for r in rows if _elem(r[i_sym]) in TM_ELEMS]
    if not tm_rows:
        raise ValueError(f"模板无 TM 位: {template}")
    anchor = max(tm_rows, key=lambda r: float(r[i_occ].split("(")[0]))
    x, y, z = anchor[i_x], anchor[i_y], anchor[i_z]

    # 丢掉锚点位上**全部**原有原子（含该位的少量 Li —— 由名义配比重建）
    kept = [r for r in rows if (r[i_x], r[i_y], r[i_z]) != (x, y, z)]

    extra = []
    for n, (el, occ) in enumerate(tm_occ.items()):
        row = list(anchor)
        row[i_lbl] = f"{el}TM{n + 1}"
        row[i_sym] = ION.get(el, el)
        row[i_x], row[i_y], row[i_z] = x, y, z
        row[i_occ] = f"{occ:g}"
        extra.append(row)

    body = [" ".join(r) for r in kept + extra]

    out: list[str] = []
    for ln in lines[:start]:
        if ln.strip().startswith("_chemical_formula_sum"):
            out.append(f"_chemical_formula_sum            '{formula}'")
        else:
            out.append(ln)
    out.extend(body)
    out.extend(lines[end:])

    head = [
        f"# PolyXRD derived from COD {template.stem}",
        f"#   {note}",
        "#   仅改 _chemical_formula_sum 与 TM 位占位；对称操作与其余原子保持模板原值。",
        "#   晶胞起点取模板值，精修/标定会自行修正（NCA 实测多在 a≈2.87、c≈14.19）。",
        "#   非实验测定结构，供精修初值使用。",
        "",
    ]
    return "\n".join(head + out) + "\n"


def _copy_entry(name: str, cid: int, note: str, dst_dir: Path) -> dict | None:
    """原样复制一份 COD CIF（规范化书写），返回索引条目。"""
    src = cod_path(cid)
    if not src.exists():
        print(f"[缺失] COD {cid} -> {src}", file=sys.stderr)
        return None
    dst = dst_dir / name
    raw = src.read_text(encoding="utf-8", errors="replace")
    dst.write_text(normalize_cif_text(raw), encoding="utf-8")
    sg, a, c = read_cell_and_formula(dst.read_text(encoding="utf-8").splitlines())
    return {"name": dst.stem, "file": name, "cod_id": cid, "spacegroup": sg,
            "a": a, "c": c, "derived": False, "normalized": True,
            "sha256": sha256(dst),
            "cod_sha256": hashlib.sha256(raw.encode("utf-8")).hexdigest(),
            "note": note}


def _derive_entry(name: str, cid: int, formula: str,
                  tm_occ: dict[str, float], note: str, dst_dir: Path) -> dict | None:
    src = cod_path(cid)
    if not src.exists():
        print(f"[缺失] 模板 COD {cid}", file=sys.stderr)
        return None
    dst = dst_dir / name
    text = normalize_cif_text(derive(src, formula, tm_occ, note))
    dst.write_text(text, encoding="utf-8")
    sg, a, c = read_cell_and_formula(text.splitlines())
    return {"name": dst.stem, "file": name, "cod_id": cid, "spacegroup": sg,
            "a": a, "c": c, "derived": True, "normalized": True,
            "formula": formula, "sha256": sha256(dst),
            "cod_sha256": hashlib.sha256(src.read_bytes()).hexdigest(),
            "note": note}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="构建 NCA 内置 CIF 包")
    ap.add_argument("--check", action="store_true", help="只校验，不写文件")
    args = ap.parse_args(argv)

    if args.check:
        return _check_only()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    IMP_DIR.mkdir(exist_ok=True)
    entries: list[dict] = []
    imps: list[dict] = []

    for name, (cid, note) in COPY_LIST.items():
        e = _copy_entry(name, cid, note, OUT_DIR)
        if e is None:
            return 1
        entries.append(e)

    for name, (cid, formula, tm_occ, note) in DERIVE_LIST.items():
        e = _derive_entry(name, cid, formula, tm_occ, note, OUT_DIR)
        if e is None:
            return 1
        entries.append(e)

    for name, (cid, note) in IMPURITY_LIST.items():
        e = _copy_entry(name, cid, note, IMP_DIR)
        if e is None:
            return 1
        imps.append(e)

    entries.sort(key=lambda e: (e["derived"], e["name"]))
    INDEX.write_text(json.dumps(
        {"cod_root": str(COD_ROOT), "count": len(entries),
         "entries": entries, "impurities": imps},
        ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"内置 CIF: {len(entries)} 份（原样 {len(COPY_LIST)}，派生 {len(DERIVE_LIST)}）"
          f" + 杂质 {len(imps)} 份")
    for e in entries + imps:
        flag = "派生" if e["derived"] else "原样"
        print(f"  [{flag}] {e['file']:<26} COD {e['cod_id']:<9} "
              f"{e['spacegroup']:<12} a={e['a']:.4f} c={e['c']:.4f}")
    return 0


def _check_only() -> int:
    """只校验：INDEX 与磁盘文件 sha256 一致。"""
    if not INDEX.exists():
        print("缺 INDEX.json", file=sys.stderr)
        return 1
    idx = json.loads(INDEX.read_text(encoding="utf-8"))
    bad = 0
    for e in idx["entries"] + idx.get("impurities", []):
        d = OUT_DIR if (OUT_DIR / e["file"]).exists() else IMP_DIR
        p = d / e["file"]
        if not p.exists():
            print(f"[缺文件] {e['file']}")
            bad += 1
        elif sha256(p) != e["sha256"]:
            print(f"[sha256 不符] {e['file']}")
            bad += 1
    print(f"{'OK' if not bad else 'FAIL'}  {idx['count']} 份主件 + "
          f"{len(idx.get('impurities', []))} 份杂质，异常 {bad}")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
