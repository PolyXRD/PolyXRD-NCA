# PolyXRD NCA

**A dedicated tool for XRD powder diffraction refinement of ternary lithium-ion cathode materials (NCA, nickel-cobalt-aluminium) and their precursors.**

Member of the PolyXRD family: [PolyXRD](https://github.com/PolyXRD) (data analysis) ·
PolyXRD-Converter (format conversion) · PolyXRD-SynthPattern (pattern synthesis) ·
PolyXRD-RutAna (TiO₂ polymorph ratio) · PolyXRD-DoCryst (degree of crystallinity) ·
PolyXRD-Graphite (degree of graphitization) · PolyXRD-NCM (NCM ternary) · **NCA (this project)**

## What it is

One window, right from the start: the pattern on the left (pattern reading + smoothing / background subtraction / Kα2 stripping / peak finding / single-peak fitting),
and built-in CIFs on the right (four NCA composition ratios 80:15 / 82:13 / 85:10 / 90:05, comparison layered phases,
precursors, aluminium-based secondary phases, impurities), with support for importing your own CIFs.

After importing XRD data, it **first decides "is it layered NCA"**: if not, a dialog stops and explains the reason;
if yes, it performs full-pattern Rietveld refinement starting from the CIF and outputs

- Unit cell a / b / c / α / β / γ, crystal system / point group / space group / Laue group
- Crystallite size (Scherrer, coherent diffraction domain) and strain
- Li occupancy (3a site) and the z coordinate of O, Li/Ni cation mixing
- Rwp / Rp / Rexp / RB / GOF / χ², and the parameter correlation matrix
- Pattern and result export (PNG / CSV / JSON)

> Status: **v0.1.0** (feature-complete, 431+ tests). For tasks and criteria see
> `docs/TASKS.md`, `docs/STD_FREEZE.md`.

## Relationship with PolyXRD-NCM (important)

This project is **ported from PolyXRD-NCM v0.1.0 as its skeleton**: the two are structurally isomorphic (R-3m #166 layered
lithium transition-metal oxide), the only difference being whether the 3b site is **Al** or **Mn**. During the port **not a single criterion threshold
was changed**; what changed is the data side (CIF library / candidate pool / sample library / real anchors).

Therefore:

- **Mechanistic conclusions** (identification gate, unit-cell recovery, occupancy convergence, round-trip ground-truth reversibility) are shared by both versions;
- **Absolute numerical values** must be re-anchored under the NCA convention; the NCM numbers cannot be used as the NCA baseline;
- The family data root **contains no measured NCA pattern**; the only usable real ternary layered cathode sample is NCM811
  (`XRData/csuHJW/Data040`), which is used as an **isomorphic anchor**, with the cost quantified and recorded in
  `docs/STD_FREEZE.md` §0.

## Capability boundaries (also stated in the UI)

1. The identification gate only guarantees "it is a layered R-3m lithium transition-metal oxide structure type". **The unit cells of all NCA compositions /
   NCM / LiNiO₂ are almost identical, with peak-intensity differences of 0.03–0.14 %, which XRD cannot distinguish** —
   the Ni:Co:Al ratio must come from ICP/EDS.
2. The "crystallite size" from XRD is the **coherent diffraction domain size**, not the secondary-particle size from SEM.
3. The absolute accuracy of Li occupancy is about ±0.02–0.03; it is suitable for **relative comparison within the same batch**, not as an absolute ground truth.

## Documentation

| Document | Content |
|---|---|
| `docs/PROJECT_PLAN.md` | Scheme, feasibility measurements, phased roadmap |
| `docs/TASKS.md` | Atomic task cards and acceptance commands |
| `docs/STD_FREEZE.md` | Criteria and convention freeze (**§0 is the NCA/NCM difference register, read it first**) |
| `docs/CIF_SOURCE.md` | COD sources and derivation rules of the built-in CIFs |
| `docs/VALIDATION*.md` | Real-data validation records (including NCA-convention re-anchoring) |
| `docs/UI_GUIDE.md` | UI structure, wording and interaction conventions |
| `docs/GSAS2_STATUS.md` | Current status of GSAS-II cross-validation |

## Running

```bash
D:/Project/XRD/python310/python.exe -m nca.gui                  # GUI
D:/Project/XRD/python310/python.exe -m nca.cli.nca --help       # command line
```

## Running tests

```bash
cd /d/Project/PolyXRD-Tools/PolyXRD-NCA
D:/Project/XRD/python310/python.exe -m pytest -q                  # full suite
D:/Project/XRD/python310/python.exe -m pytest -q -m "not slow"    # skip full-pattern refinement
```

## Packaging

```bat
build.bat             :: = build_exe.bat + build_setup.bat, artifacts collected into dist\
tidy_dist.bat dry     :: list only (in Git Bash you must write dry, not /dry)
```

`dist\` keeps only two deliverables: `PolyXRD-NCA-Internal-vX.Y.Z.exe` and
`PolyXRD-NCA-vX.Y.Z-Setup.exe`.

## License

Dual license (Academic / Commercial), see `LICENSE` and `LICENSE-CN`.
For commercial licensing contact: sshztx@outlook.com · https://github.com/PolyXRD
