# PolyXRD NCA

**三元锂电正极材料（NCA，镍钴铝）及其前驱体 XRD 粉末衍射精修专用工具。**

PolyXRD 家族成员：[PolyXRD](https://github.com/PolyXRD)（数据分析）·
PolyXRD-Converter（格式转换）· PolyXRD-SynthPattern（谱图合成）·
PolyXRD-RutAna（TiO₂ 晶型比）· PolyXRD-DoCryst（结晶度）·
PolyXRD-Graphite（石墨化度）· PolyXRD-NCM（NCM 三元）· **NCA（本项目）**

## 这是什么

打开即一个界面：左侧谱图（读谱 + 平滑/扣背景/去 Kα2/寻峰/单峰拟合），
右侧内置 CIF（NCA 四档配比 80:15 / 82:13 / 85:10 / 90:05、对照层状相、
前驱体、铝系副相、杂质），并支持导入自己的 CIF。

导入 XRD 数据后**先判「是不是层状 NCA」**：不是则弹窗停止并说明理由；
是则以 CIF 为起点做全谱 Rietveld 精修，输出

- 晶胞 a / b / c / α / β / γ、晶系 / 点群 / 空间群 / 劳厄群
- 晶粒尺寸（Scherrer，相干衍射畴）与应变
- Li 占位（3a 位）与 O 的 z 坐标、Li/Ni 混排
- Rwp / Rp / Rexp / RB / GOF / χ²，以及参数相关性矩阵
- 谱图与结果导出（PNG / CSV / JSON）

> 状态：**v0.1.0**（功能完整，431+ 测试）。任务与判据见
> `docs/TASKS.md`、`docs/STD_FREEZE.md`。

## 与 PolyXRD-NCM 的关系（重要）

本项目以 **PolyXRD-NCM v0.1.0 为骨架移植**：两者结构同构（R-3m #166 层状
锂过渡金属氧化物），差别只在 3b 位是 **Al** 还是 **Mn**。移植时**判据阈值
一个没动**，改的是数据面（CIF 库 / 候选池 / 样本库 / 真实锚点）。

因此：

- **机制性结论**（判定门、晶胞回收、占位收敛、round-trip 真值可逆）两版共用；
- **绝对数值**必须在 NCA 口径下重新锚定，不能把 NCM 的数当 NCA 的基准；
- 家族数据根里**没有 NCA 实测谱**，唯一可用的三元层状正极实样是 NCM811
  （`XRData/csuHJW/Data040`），它被当作**同构锚点**使用，代价已量化登记在
  `docs/STD_FREEZE.md` §0。

## 能力边界（界面里也会写）

1. 判定门只保证"是层状 R-3m 锂过渡金属氧化物结构类型"。**NCA 各配比 /
   NCM / LiNiO₂ 的晶胞几乎相同、峰强差 0.03–0.14 %，XRD 分不开** ——
   Ni:Co:Al 比例须由 ICP/EDS 给。
2. XRD 的"晶粒尺寸"是**相干衍射畴尺寸**，不是 SEM 的二次颗粒粒径。
3. Li 占位的绝对值精度约 ±0.02–0.03，适合**同批次相对比较**，不宜当绝对真值。

## 文档

| 文档 | 内容 |
|---|---|
| `docs/PROJECT_PLAN.md` | 方案、可行性实测、阶段路线图 |
| `docs/TASKS.md` | 原子任务卡与验收命令 |
| `docs/STD_FREEZE.md` | 判据与口径冻结（**§0 是 NCA/NCM 差异登记，先读它**） |
| `docs/CIF_SOURCE.md` | 内置 CIF 的 COD 来源与派生规则 |
| `docs/VALIDATION*.md` | 真实数据验证记录（含 NCA 口径重锚） |
| `docs/UI_GUIDE.md` | 界面结构、措辞与交互约定 |
| `docs/GSAS2_STATUS.md` | GSAS-II 交叉验证的现状 |

## 跑起来

```bash
D:/Project/XRD/python310/python.exe -m nca.gui                  # 界面
D:/Project/XRD/python310/python.exe -m nca.cli.nca --help       # 命令行
```

## 跑测试

```bash
cd /d/Project/PolyXRD-Tools/PolyXRD-NCA
D:/Project/XRD/python310/python.exe -m pytest -q                  # 全量
D:/Project/XRD/python310/python.exe -m pytest -q -m "not slow"    # 跳过全谱精修
```

## 打包

```bat
build.bat             :: = build_exe.bat + build_setup.bat，产物收进 dist\
tidy_dist.bat dry     :: 只列清单（Git Bash 里必须写 dry，不能写 /dry）
```

`dist\` 只保留两个交付物：`PolyXRD-NCA-Internal-vX.Y.Z.exe` 与
`PolyXRD-NCA-vX.Y.Z-Setup.exe`。

## 许可

双许可（Academic / Commercial），见 `LICENSE` 与 `LICENSE-CN`。
商用联系：sshztx@outlook.com · https://github.com/PolyXRD
