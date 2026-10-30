# PolyXRD-NCA 项目方案与分步路线图

> 三元锂电正极材料（NCA）XRD 粉末衍射精修专用工具 · 项目立项文档
>
> 根目录：`D:\Project\PolyXRD-Tools\PolyXRD-NCA`
> 文档：`docs/PROJECT_PLAN.md`（本文）、`docs/TASKS.md`（原子任务卡）、`docs/CIF_SOURCE.md`（CIF 来源与出处）
> 创建日期：2026-10-29 · 版本 v0.1
>
> ⚠ **本文是立项期的方案文档（以 NCM 骨架讨论），数字与型号清单为当时口径。**
> 落地后的最终口径以 `docs/STD_FREEZE.md` §0/§7（NCA/NCM 差异登记）、
> `docs/CIF_SOURCE.md` v2（CIF 构成）、`docs/VALIDATION*.md` §v6（NCA 锚点）为准。

---

## 0. 一句话定位

打开即一个界面：左侧谱图（读谱 + 平滑/扣背景/去 Kα2/寻峰/单峰拟合），右侧内置
层状 NCA 三元正极材料的 COD CIF（NCA 四档配比 + 对照层状相 + 前驱体 + 铝系副相）。
导入 XRD 数据后
先做「是不是层状 NCA」判定——不是则弹窗停止；是则以 CIF 为起点做全谱 Rietveld
精修，输出晶胞参数、晶粒尺寸、晶系/点群/空间群/劳厄群、Li 占位与 Rwp/RB/GOF，
支持导出谱图与结果。家族风格（logo/splash/三语/注册/安装包/版权）继承 PolyXRD。

---

## 1. 可行性结论（已实测，非设想）

用家族现有内核（`PolyXRD-RutAna` B 法 Rietveld 引擎）**零改动**直接拟合真实
三元层状正极数据（`D:\Project\XRD\XRData\csuHJW\Data040：NCM：811_raw.txt`，3500 点，
10–80°，步长 0.02°；**该样品是 NCM811，家族数据根里没有 NCA 实测谱**，
它被当作同构锚点使用，代价见 `docs/STD_FREEZE.md` §0），实测结果：

| 指标 | 实测值 | 结论 |
|---|---|---|
| 理论峰位命中 | 11 个主峰全中，偏差 ≤0.03°（文献晶胞未精修） | 结构模型正确 |
| (006)/(012)、(108)/(110) 分裂 | 64.52° / 64.94° 清晰双峰 | 层状有序度判据存在 |
| I(003)/I(104) | 1.31（>1.2） | 阳离子混排低 |
| 整谱精修 Rwp | **7.87%**（Rp 4.33%，47 s / 940 迭代） | 引擎可用 |
| 晶胞回收 | a=2.8718，c=14.1942，c/a=4.9426 | 与文献吻合 |
| 暴露缺口 ① | Caglioti U 收敛到 0、经验衰减 B=3.32 | **缺尺寸-应变宽化项**，被迫用高角衰减补偿 |
| 暴露缺口 ② | `rbragg = None` | **RB 未实现** |

复现脚本在 `_spike/`（`probe_r3m.py`、`spike_refine.py`），未改动任何家族仓库。

**结论：可行，且内核复用度高达 60–70%。主要新增工作在四块：①尺寸-应变峰形
（TCH）②占位与 O 坐标精修 ③RB/GOF 与空间群元数据 ④"是不是层状 NCA"判定门。**

---

## 2. 需求裁定（本轮已拍板）

| # | 决策点 | 裁定 | 备注 |
|---|---|---|---|
| 1 | 内置 CIF 来源 | **COD 实测数据**（`sync-COD` 已同步） | 候选清单见 §3 |
| 2 | 判定门严格度 | **初版只挡明显非层状相**（稳妥） | 判据见 §4 |
| 3 | 第二相槽位 | **加**（Li₂CO₃/NiO/碳等杂质相） | +3–4 人日，见 P5 |
| 4 | GSAS-II 交叉验证 | **接**（开发期用，不打包） | 复用主程序桥；**当前未安装，需先装**，见 §6 |

---

## 3. 内置 CIF：COD 来源与候选（已实查 `inorg_manifest.sqlite`）

`sync-COD` 数据根 `D:\Project\XRD\sync-COD`，`inorg_manifest.sqlite` 的 `files`
表（`cod_id, lib, rel, formula, ...`，1,959,421 行）已按 `is_inorg=1 AND lib="COD"`
实查。

**实查结论（决定了整个 CIF 策略）**：COD 里**没有任何 Ni-Co-Al 层状氧化物**
条目（`Ni+Co+Al+O` 组合零命中），也没有 Ni-Co-Al 氢氧化物/碳酸盐。因此
NCA 四档配比只能"派生"，可用的**真实测定结构**只有下面这些骨架与对照：

| 用途 | 文件 | COD 编号 | 化学式（formula_sum） | a / c（Å） |
|---|---|---|---|---|
| **骨架（四档派生母本）** | `NCA_*.cif` ×4 | `1520789` | Co0.1 Li1.03 Mn0.1 Ni0.77 O2 | 2.8645 / 14.161 |
| **Ni 端元** | `LiNiO2.cif` | `1533770` | Li Ni O2 | 2.8755 / 14.1806 |
| **判定门负样本** | `LiCoO2.cif` | `1533825` | Co Li O2 | 2.8140 / 14.0480 |
| **脱锂态参考** | `layered_delith.cif` | `4002445` | Li0.63… O2 | 2.8220 / 14.4630 |
| **铝偏析相** | `LiAlO2.cif` | `4002658` | Li Al O2 | 5.1965 / 6.3464（P4₁2₁2） |
| **铝源残留** | `AlOH3_gibbsite.cif` | `1200016` | Al H3 O3 | 8.6840 / 9.7360（P2₁/n） |
| **铝源脱水相** | `AlOOH_boehmite.cif` | `9012247` | Al O OH | 2.8510 / 3.7360（Cmcm） |
| **氢氧化物前驱体骨架** | `precursor_NiOH2.cif` | `1548811` | Ni H2 O2 | 3.1268 / 4.6060（P-3m1） |
| **碳酸盐前驱体骨架** | `precursor_NiCO3.cif` | `1548823` | Ni C O3 | 4.6117 / 14.7350（R-3c） |
| **杂质（残碱）** | `impurity/Li2CO3.cif` | `9008283` | Li2 C O3 | 8.3593 / 6.1975（C2/c） |
| **杂质（岩盐相）** | `impurity/NiO.cif` | `1010093` | Ni O | 4.1684（Fm-3m） |

**原子占位示例（COD 1520789，逐行核验）**：Li 层（z=0.5）= Li 0.98 + Ni 0.02；
TM 层（z=0）= Ni 0.75 + Co 0.10 + Mn 0.10 + Li 0.05；O 在 6c，z=0.259。即**该
CIF 自带 Li/Ni 混排与完整的「LiMO₂」化学计量**，是理想的派生母本。注意它本身
是 NCM811 的实测结构 —— 派生时把 TM 位的 Mn 换成 Al 并按配比重排（见
`docs/CIF_SOURCE.md` §2），骨架与对称操作不动。

### 3.1 派生策略（T0.3 落地）

- **NCA 四档配比**（80:15 / 82:13 / 85:10 / 90:05）：以 COD 1520789 为母本，
  只改 TM 位元素与占位、只改 `_chemical_formula_sum`，头部标
  `# PolyXRD derived from COD 1520789`。
  配平口径：Li 层继承 `Li 0.98 + Ni 0.02`，TM 位显式给 `Li 0.02` → 总 Li = 1.00、
  总 Ni = 目标值（Hill–Howard 定量要用）。
- **前驱体**：Ni 端元骨架（1548811 / 1548823）+ 名义配比 80:15:5 改写 TM 位。
- **四档晶胞故意完全相同（2.8645 / 14.161）**：NCA 各配比的真实晶胞差在
  0.01 Å 量级、峰强差 0.03–0.14 %，XRD 分不开；给每档编不同晶胞等于制造
  "看起来能区分型号"的假象。
- **删掉的条目**：NCM 版的 111 / 111_alt / 523 / 622 / 811_alt / 富锂等全部移除
  —— 留着只会让判定门把 NCA 谱错配到别人的晶胞上。

### 3.2 CIF 归口与许可

- COD 数据许可：**CC0（公有领域）**，可随软件分发、可再加工，但需在「关于/致谢」
  页保留 COD 归属与条目出处（`CIF_SOURCE.md` 逐条登记 COD 编号 + 原始文献 DOI）。
- 自建/改写 CIF 需在文件头注释来源与改写内容，保证可溯源。

---

## 4. 判定门设计（初版：只挡明显非层状相）

### 4.1 输入与流程

读谱 → 自动寻峰（+去 Kα2 后）→ 提取前 12 个强峰 2θ → 与 R-3m 理论峰位序列比对
→ 三重判据（全过才算"是层状 NCA"）→ 不过则弹窗列出未过项与理由，**停止精修**。

### 4.2 三重判据（初版阈值，P0 冻结，P5 用正负样本标定）

| 判据 | 指标 | 初版阈值 | 依据 |
|---|---|---|---|
| C1 峰位序列匹配 | 前 12 强峰与 R-3m 理论峰的 2θ 均方根偏差 | RMS ≤ 0.10°（含 003/101/006/012/104/015/107/018/110/113 主干峰） | 实测 11 峰偏差 ≤0.03°，阈值留 3 倍余量 |
| C2 特征分裂 | (006)/(012) 双峰（≈38°）与 (108)/(110) 双峰（≈64.5°）存在 | 各自能检出两个可分辨峰（Δ2θ ≈ 0.15–0.4°） | 层状有序度核心判据 |
| C3 峰强比 | I(003)/I(104) | 0.5 < R < 3.0（宽松带，不判混排，只挡非层状） | 层状 NCA 典型 1.1–1.6；尖晶石/岩盐相会显著偏离 |

> C3 只做"是否层状"的粗判，**不做**混排程度判定——Li/Ni 混排是精修输出项，
> 不在此处卡。

### 4.3 初版「明显非层状相」白名单/黑名单（稳妥策略）

| 类 | 判据 | 动作 |
|---|---|---|
| LCO（LiCoO₂） | a ≈ 2.814 vs NCA a ≈ 2.86–2.88，主峰整体右移 ~0.5–0.8° | 拦下，提示"疑似钴酸锂/其他层状，非 NCA 三元" |
| 尖晶石（LiMn₂O₄ 等 Fd-3m） | 峰位序列对不上 R-3m；出现 (311)/(222) 特征峰 | 拦下 |
| 橄榄石（LFP，Pnma） | 峰位序列对不上；d 值体系差异大 | 拦下 |
| 岩盐/无序相 | C2 特征分裂消失（(006)/(012) 合并） | 拦下，提示"可能已发生岩盐化/过度放电" |

> 边界声明（写入界面措辞）：NCA 四档配比 / NCM / LiNiO₂ 的晶胞几乎相同
> （a 差 <0.01 Å，峰强差 0.03–0.14 %），实验室 XRD **分不开**；判定门只保证
> "是层状 R-3m 锂过渡金属氧化物结构类型"，**不保证 Ni:Co:Al 具体比例**
> ——那需要 ICP/EDS。

### 4.4 前驱体判定（独立子门）

前驱体（氢氧化物 P-3m1 / 碳酸盐 R-3c）与成品 NCA 是不同结构，**不能共用同一
精修模板**。导入后先按空间群指纹分流：
- 若峰位匹配 P-3m1（层状氢氧化物）→ 走前驱体模板与判定；
- 若匹配 R-3c（碳酸盐）→ 走前驱体模板与判定；
- 若匹配 R-3m → 走 NCA 精修主线。
（前驱体的具体 CIF 与判定阈值在 P0 的 T0.3 一并落地，本文只立分流框架。）

---

## 5. 技术架构与复用边界

```
PolyXRD-NCA/
├── nca/
│   ├── core/                 # 纯算法，无 GUI 依赖（可 pytest）
│   │   ├── patterns.py       # 读谱（vendor 自 RutAna patterns.py）
│   │   ├── preprocess.py     # 平滑/扣背景/去 Kα2（vendor 自 RutAna + PolyXRD）
│   │   ├── peakfind.py       # 寻峰 + pseudo-Voigt 单峰拟合（vendor 自 PolyXRD）
│   │   ├── identify.py       # ★ 判定门（§4，新建）
│   │   ├── profile/          # CIF 解析/结构因子/前向（vendor 自 RutAna profile/）
│   │   │   ├── cif.py  struct.py  shape.py  forward.py
│   │   │   └── size_strain.py    # ★ TCH 尺寸-应变峰形（新建，替代纯 Caglioti）
│   │   ├── refine/           # 精修引擎（vendor 自 RutAna refine/）
│   │   │   ├── linear.py  lm.py  nonlinear.py  rfactors.py
│   │   │   ├── occ.py            # ★ 3a 位 Li/Ni 占位 + O(z) 精修（新建）
│   │   │   └── rbragg.py         # ★ Bragg R 因子（新建，填 rfactors 空缺）
│   │   ├── crystallite.py    # 晶粒尺寸（TCH 积分宽度 → Scherrer/Williamson-Hall）
│   │   ├── symmetry.py       # ★ 空间群→点群/劳厄群/晶系映射表（新建，230 群静态表）
│   │   └── report.py         # 结果卡 + 谱图/结果导出（改造）
│   ├── gui/                  # PySide6 单窗口（左谱图+工具条 / 右 CIF 面板）
│   │   ├── main_window.py  canvas.py  cif_panel.py  resultcard.py
│   │   └── identify_dialog.py    # 判定不过弹窗
│   ├── cli/                  # CLI（无 GUI 头跑判定+精修，供 GSAS-II 交叉验证）
│   ├── i18n/                 # 三语（vendor 自 RutAna，改词条）
│   ├── edition.py            # 版次（注册版/内部评估版）
│   └── resources/            # logo / splash / icon
├── tools/                    # gen_virtual / bench / cross_gsas2（交叉验证）
├── scripts/                  # gsas2_bridge（vendor 自主程序，改 NCA 精修模式）
├── tests/                    # pytest + golden 基线
├── docs/                     # PROJECT_PLAN / TASKS / CIF_SOURCE / VALIDATION
├── requirements.txt  LICENSE*  README.md  *.spec  build*.bat
```

**复用原则**（继承 RutAna 的 vendor 约定）：
- 算法模块以 **vendor 副本**随本项目发布，不在运行时 import 兄弟仓库；
- 改动处加注释 `# PolyXRD-NCA 修改`，原样处保留 `# vendored from ...` 与版权。

**回归基线（冻结）**：`tests/golden/nca811_baseline.json` 锁定
`Rwp=3.1607% / a=2.87185 / c=14.19502`（v6 重锚；NCM 版为 3.0409 / 同一晶胞），
作为后续所有改动（TCH/占位/RB）的逐位回归锚点。

---

## 6. GSAS-II 交叉验证方案

### 6.1 前置：GSAS-II 当前**未安装**

主程序已有 `scripts/gsas2_bridge.py` 子进程桥 + `_find_gsas2_python()` 定位逻辑，
但本机 `POLYXRD_GSAS2_PYTHON` 为空、常见目录（`D:\GSASII` / `C:\GSASII` /
`ProgramData\gsas2main`）均不存在。**P0 第一步就是安装 GSAS-II**（官方
gsas2main 安装器，自带 Python），否则交叉验证无从谈起。

### 6.2 交叉验证口径（只用于开发期，不打包进 exe）

| 项 | 做法 |
|---|---|
| 调用方式 | JSON 子进程（`gsas2_bridge.py`，复用主程序模式，扩展 `refine` 支持占位/O 坐标） |
| 比对对象 | 同一样品、同一 CIF，两套引擎分别精修 → 比对晶胞 a/c、Li 占位、Rwp/GOF |
| 验收判据 | 晶胞 a/c 相对偏差 ≤0.05%；Rwp 差 ≤2 个百分点；Li 占位差 ≤0.03 |
| 样本 | 合成 round-trip（已知真值）+ 真实三元层状锚点（Data040）+ 真实前驱体（Data043） |
| 排除打包 | `*.spec` 里 exclude GSAS-II 依赖；CLI 通过环境变量 `POLYXRD_GSAS2_PYTHON` 定位，未装则跳过交叉验证并明确提示 |

---

## 7. 分步实现与验收（阶段表）

### 阶段总览

| 阶段 | 内容 | 人日 | 可验收产出 |
|---|---|---|---|
| P0 | 立项冻结 + 环境（GSAS-II 装）+ CIF 归口 | 1.5–2 | CIF_SOURCE.md、GSAS-II 就绪、golden 基线 |
| P1 | 内核搬运 + 基线复现 | 2–3 | CLI 复现 7.87%，golden 绿 |
| P2 | TCH 尺寸-应变峰形 | 3–5 | 合成数据晶粒尺寸回收 ≤5%，B 因子回归物理合理 |
| P3 | 占位 + O 坐标精修 + 相关性矩阵 | 3–5 | Li 占位 round-trip ≤0.03，参数相关性可导出 |
| P4 | RB/GOF + 空间群元数据 + 结果卡 | 2–3 | RB 实现、230 群映射表、结果卡齐全 |
| P5 | 判定门 + 第二相槽位 + 前驱体分流 | 3–5 | 正负样本判定正确率 100%，第二相可精修 |
| P6 | 单窗口 GUI | 5–7 | 左谱图/右 CIF 面板，全流程可点 |
| P7 | 导出 + 家族化 | 2–3 | logo/splash/i18n/注册/安装包 |
| P8 | 验证 + GSAS-II 交叉 + 回归 | 4–6 | VALIDATION 报告，交叉比对表 |
| **合计** | | **25–39 人日** | 全投入 3–4 周 v1.0 |

### 7.1 P0 · 立项冻结与环境（1.5–2 人日）

- **T0.1** 仓库骨架：`requirements.txt`（numpy 2.2.6 / scipy 1.15.3 /
  matplotlib 3.10.9 / PySide6 6.11.2 / pytest）、`pytest.ini`、包目录 `nca/`、
  `version.py`、`README.md`、双许可 `LICENSE`/`LICENSE-CN`。
  - 验收：`python -m pytest -q` 空跑绿。
- **T0.2** 标准与判据冻结：`docs/STD_FREEZE.md` 冻结判定门阈值（§4.2）、R 因子口径、
  尺寸-应变模型选择（TCH 伪 Voigt vs 单 Lorentz 应变）、Li 占位精修边界。
  - 验收：文档逐字冻结，阈值唯一出处。
- **T0.3** CIF 归口：`docs/CIF_SOURCE.md` 登记 §3 候选 COD 编号 + DOI；
  派生 NCA 四档配比（COD 1520789 母本改写 TM 位）与前驱体（氢氧化物 P-3m1 /
  碳酸盐 R-3c 骨架 + 名义配比）；
  把选定 CIF 落盘到 `nca/data/cif/`。
  - 验收：每份 CIF 可被 `parse_cif` 解析且原子位/占位/空间群与 COD 原文一致。
- **T0.4** 安装 GSAS-II 并打通 `gsas2_bridge.py`（用 NCA_80_15_05 跑一次 Le Bail 冒烟）。
  - 验收：`POLYXRD_GSAS2_PYTHON` 指向成功，桥返回 `ok:true`。

### 7.2 P1 · 内核搬运 + 基线复现（2–3 人日）

- **T1.1** vendor 读谱/预处理/寻峰/CIF/前向/精修模块（§5 清单）复制进 `nca/core/`，
  保留版权头。
- **T1.2** CLI 入口 `nca/cli/nca.py`：`--sample --cif --out.json --png`，
  复现 §1 的 7.87%。
- **T1.3** 冻结 golden 基线 + `tools/bench.py` 回归脚本。
  - 验收：`python -m pytest` 全绿，`Rwp=7.87±0.1 / a=2.8718 / c=14.1942` 逐位一致。

### 7.3 P2 · TCH 尺寸-应变峰形（3–5 人日）

- **T2.1** `size_strain.py`：Thompson-Cox-Hastings（TCH）伪 Voigt，FWHM 由
  `FWHM = (U·tan²θ + V·tanθ + W + Z/cos²θ)^½ + (X·tanθ + Y/cosθ)` 合成，
  晶粒尺寸与应变由 X/Y（Lorentz 洛伦兹项）反演；替代/叠加纯 Caglioti。
- **T2.2** `crystallite.py`：由积分宽度反演晶粒尺寸（Scherrer 球径 + 可选
  Williamson-Hall 分离尺寸/应变）。
- **T2.3** 合成数据 round-trip：`tools/gen_virtual.py` 生成已知晶粒尺寸（如
  50/100/200 nm）的虚拟 NCA 谱，验证回收误差 ≤5%。
  - 验收：`tests/test_size_strain.py` 回收误差 ≤5%，真实数据 B 因子回到物理合理
    值（<0.5），Rwp 不劣化（≤8%）。

### 7.4 P3 · 占位与 O 坐标精修 + 相关性矩阵（3–5 人日）

- **T3.1** `occ.py`：新增精修参数——3a 位 Li 占有率（约束 Li3a + Ni3a = 1，即
  Li/Ni 混排度）、3b 位过渡金属比例（可选，默认由 CIF 固定）、O 的 z 坐标；
  边界保护（occ ∈ [0.5, 1.0] 等）。
- **T3.2** 参数相关性：LM 收敛后算参数协方差矩阵（Jacobian 的伪逆），输出
  Li占位 ↔ 整体B因子 / 择优取向 / 标度 的相关性矩阵。
- **T3.3** 择优取向（00l 板状 March-Dollase）默认开启（复用 RutAna v1.4 机制）。
  - 验收：`tests/test_occ.py` 合成谱 Li 占位 round-trip ≤0.03；相关性矩阵可导出；
    真实 NCA_80_15_05 的 Li 占位与 ICP（如有）一致或在文献区间内。

### 7.5 P4 · RB/GOF + 空间群元数据 + 结果卡（2–3 人日）

- **T4.1** `rbragg.py`：Bragg R 因子，按"实测强度在重叠峰间的分配"标准口径实现
  （写明口径，与 GSAS-II 对齐用于交叉验证）；GOF = Rwp/Rexp。
- **T4.2** `symmetry.py`：230 空间群 → 晶系/点群/劳厄群/布拉维格子静态映射表
  （R-3m → 三方 / -3m(D₃d) / -3m / hR），覆盖内置 CIF 与前驱体空间群。
- **T4.3** 结果卡数据结构：晶胞 a/b/c/α/β/γ、体积、晶系、点群、空间群、劳厄群、
  晶粒尺寸、Li 占位、Rwp/Rp/Rexp/RB/GOF、择优取向 r、参数相关性。
  - 验收：`tests/test_rfactors.py` RB 与合成真值一致；结果卡字段齐全且口径可溯。

### 7.6 P5 · 判定门 + 第二相槽位 + 前驱体分流（3–5 人日）

- **T5.1** `identify.py` 判定门（§4 三重判据 + 前驱体空间群分流）。
- **T5.2** 第二相槽位：精修引擎支持 2 相（主相 NCA + 杂质相 Li₂CO₃/NiO/碳），
  杂质相 CIF 内置常用 3–5 个（Li₂CO₃、NiO、LiOH、石墨碳）。
- **T5.3** 判定门样本库：收集正样本（多型号 NCA）与负样本（LCO/尖晶石/LFP/岩盐/
  前驱体）各 ≥10 条，标定阈值并冻结。
  - 验收：`tests/test_identify.py` 正负样本判定正确率 100%；第二相可精修且 wt% 合理。

### 7.7 P6 · 单窗口 GUI（5–7 人日）

- **T6.1** `main_window.py`：左侧谱图区（matplotlib canvas）+ 工具条（平滑/扣背景/
  去Kα2/寻峰/拟合），右侧 CIF 面板（下拉选内置 + 导入 CIF + 显示晶胞/占位）。
- **T6.2** 流程编排：导入 → 判定（不过弹 `identify_dialog`）→ 精修（后台线程 +
  忙窗进度）→ 结果卡 → 导出。
- **T6.3** 三语 i18n（vendor RutAna 词条改造）+ 关于/版权/联系方式页。
  - 验收：`python -m nca.gui` 手动冒烟，全流程无卡死。

### 7.8 P7 · 导出 + 家族化（2–3 人日）

- **T7.1** 谱图导出 PNG/SVG；结果导出 TXT/CSV/HTML（含公式与判据口径）。
- **T7.2** 家族化：logo/splash（PyInstaller bootloader splash）、版次
  （注册版 Setup + 内部评估版单文件）、注册机制（复用 PolyXRD-Reg 出码）、
  Inno Setup 安装包（`build_exe.bat`/`build_setup.bat`）。
  - 验收：两版 exe 可运行，注册/未注册行为符合家族规范。

### 7.9 P8 · 验证 + GSAS-II 交叉 + 回归（4–6 人日）

- **T8.1** 合成 round-trip（多晶粒尺寸/多混排度/多择优取向）。
- **T8.2** 真实数据：三元层状锚点（Data040，NCM811）+ 真实三元氢氧化物前驱体
  （Data043，附 Maud 独立精修值）+ 四档配比合成谱。
- **T8.3** GSAS-II 交叉比对（§6.2 口径），输出 `docs/VALIDATION_vX.md` 比对表。
- **T8.4** 全量回归（golden + 判定门 + 第二相），打包验证。
  - 验收：VALIDATION 报告齐备，交叉比对判据全过，两版安装包验证通过。

---

## 8. 风险与依赖

| # | 风险 | 等级 | 缓解 |
|---|---|---|---|
| 1 | COD 无任何 Ni-Co-Al 层状 CIF | 中 | 派生四档配比（母本 COD 1520789），标注来源与"名义配比"（§3.1 / CIF_SOURCE） |
| 2 | Li 占位可辨识性差（Li 散射弱 + 与 B 因子/择优强相关） | 高 | 输出相关性矩阵；界面措辞明示精度 ±0.02–0.03 |
| 3 | 判定门误拦 | 高 | 初版只挡明显非层状（稳妥）+ 样本库标定 + 弹窗给理由 |
| 4 | GSAS-II 未安装 | 高（当前阻断） | P0 先装；未装时 CLI 跳过交叉验证并提示 |
| 5 | 尺寸-应变与占位参数耦合 | 中 | 分阶段精修（先晶胞/峰形，后占位），相关性矩阵监控 |
| 6 | 前驱体结构域扩展 | 中 | 独立分流模块，不污染 NCA 主线 |

---

## 9. 下一步（待你确认后立即执行）

1. 建任务卡 `docs/TASKS.md`（按 §7 逐卡、一卡一验收，沿 RutAna 执行规范）。
2. P0 开工：装 GSAS-II、落 CIF（811/111/LiNiO₂/LCO + 622/523 改写 + 前驱体检索）、
   冻结 STD_FREEZE。
3. 把 `_spike/spike_refine.py` 扩成正式 `nca/cli/nca.py` + golden 基线。

> 需要你确认的一点：**前驱体是否纳入 v1.0**？若纳入，P0 要额外花 ~0.5 人日
> 检索氢氧化物/碳酸盐 CIF 并定前驱体判定阈值；若先做成品 NCA，前驱体移到 v1.1。
