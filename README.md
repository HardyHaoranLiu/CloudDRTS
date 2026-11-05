# CloudDRTS：四阶段实验复核与补强 v7

本版完成自编完整入口的重新运行（13 个脚本），新增保留电容动态的戴维南等效、独立频域解、三频率/五延迟比较，以及原计划要求的分阶段说明。ANDES 原始运行来自上一版已执行的独立基准；本版重新核对它与本次自编故障 CSV 的比较。

## 四项实验现在怎样评价

| 实验 | 已有证据 | 当前问题与判断 |
|---|---|---|
| 1 机电动态 | 六阶四机、三相故障、三模态、ANDES 同参数潮流/模态/故障比较；两页说明 | 计算和比较已完成，但原论文工况、饱和/补偿对应以及教材阻尼差异未闭合，保留“文献对齐部分完成” |
| 2 EMT 基础 | R/L/C 伴随模型、固定电压源、开关、无损 Bergeron；解析响应、多个步长、反射、能量及事件连续性 | 入门练习检查通过；线路要求整数传播步，事件电容约束必须独立一致，属于已声明适用范围 |
| 3 切分/延迟/等效 | 历史项交换、单向整数步附加延迟、静态近似和动态戴维南、误差归因、三页说明和 Kundur 推广方案 | 入门练习检查通过；动态多端口机组 EMT 与真实网络通信尚属后续研究 |
| 4 实时入门 | 每步分位数、超时比例、调用热点、环境和第一周讨论清单 | 已完成计划允许的 Python 剖析路径；当前代码没有硬实时保证，HLS 是另一选项 |

**不把全部四项写成无保留验收完成。** 主要数值缺口在实验 1。两页模态说明和三页接口说明均已按原计划页数单独输出；阶段四清单也不再只散落于旧版综合报告。

## 本轮新增结果

- 400 Hz、零附加延迟下，原静态戴维南近似的受端电压误差为约 32.987%；保留动态电容后为浮点舍入量级。
- 三频率、五延迟下，独立动态等效与完整离散子电路的最大电压差为 1.85e-13 V，最大电流差为 3.11e-15 A。
- 独立连续/离散频域核对得到约 2.00042 的步长收敛阶；电压、电流分别按各自单位保存。
- 保留原静态变体用于比较；不能把模型零延迟误差算作通信延迟误差，也不能把误差抵消解释为延迟改善了算法。

## 入口和交付物

本 Git 仓库包含源码、依赖锁、必要案例和第三方许可。下列生成结果、PDF 报告、`ACCEPTANCE_STATUS.json` 和 Git bundle 随完整的 `CloudDRTS_Experiments_v7.zip` 交付，不全部纳入 Git；运行下方入口可重新生成实验结果和说明。

- `CloudDRTS/run_all.py`：完整自编实验和验证入口。
- `CloudDRTS/dynamic_thevenin.py`：新增独立标量动态等效；不调用原 Circuit / Bergeron 类。
- `CloudDRTS/results/`：全部本轮自编 CSV、JSON、Nature 风格 PNG/SVG/矢量 PDF；保留上、右边框。
- `output/pdf/Phase1_Note_EN.pdf`、`output/pdf/Phase1_Note_ZH.pdf`：阶段一两页说明，英文、中文各一份，IEEE 会议双栏排版；由 `phase1_note/build_phase1_note.py` 生成，正文数字全部取自结果文件。
- `CloudDRTS/phase1_sensitivity.py`：复算文献 [2] 表 1、表 2 的四个对称工况，对照教材网络与原文无升压变压器网络，并扫描升压电抗和补偿量（只扫描，不拟合）。
- `output/pdf/Phase3_Partitioning_Note.pdf`：阶段三三页说明。
- `Phase2_Validation.md`：阶段二验证表与适用范围。
- `Phase4_Questions.md`：本轮实测性能与第一周讨论清单。
- `BOUNDARY_CONTRACT.json`：电流方向、单位、交换顺序和单向/往返延迟定义。
- `ACCEPTANCE_STATUS.json`：四项状态、检查数值和未完成范围。
- `READING_MAP.md`：13 条指定阅读与实验问题的对应，不声称替代个人阅读和导师讨论。
- `CloudDRTS_Source_v7.bundle`：可克隆 Git 源码仓库，使用方法见 `Git_Reproduction.md`。
- `external_run/`：上一轮实际 ANDES 原始日志、输入、依赖锁与复跑程序；本版已更新相对路径并重算比较统计。
- `external_run/output/pdf/CloudDRTS_External_Benchmark_Report.pdf`：上一版六页外部核对记录，作为历史报告保留；最新四项状态以本 README 和分阶段说明为准。

## 复跑

自编实验使用 Python 3.12，单独环境：

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r CloudDRTS/requirements.txt
.venv/bin/python CloudDRTS/run_all.py
```

ANDES 基准使用另一个环境，避免它的 NumPy/SciPy 锁与自编实验锁混装：

```bash
python3.12 -m venv .venv-andes
.venv-andes/bin/python -m pip install -r external_run/requirements_lock.txt
.venv-andes/bin/python external_run/run_reproduction.py
```

Git 源码交付不带所有大体积输出；首次从源码复跑请先运行自编完整入口，再运行 ANDES 入口。若仅运行 ANDES，追加 `--rerun-local-transient` 可先生成其依赖的自编故障波形。

Windows 将 `.venv/bin/python` 换为 `.venv\Scripts\python.exe`。生成阶段三说明另安装 `requirements_notes.txt`，再执行 `python build_delivery_notes.py`。

阶段一中英文说明需要 XeLaTeX（IEEEtran、ctex/xeCJK、siunitx）以及 TeX Gyre Termes、Noto Serif CJK SC 字体；先完成上面两个入口，再在自编环境中执行：

```bash
.venv/bin/python -m pip install -r requirements_notes.txt
.venv/bin/python phase1_note/build_phase1_note.py
```

构建会检查：页数为 2、LaTeX 无未定义引用或溢出、ANDES 比较所用波形与本次 `detailed_fault.csv` 哈希一致，以及正文结论对应的数值条件。

## 实验 1 剩余问题的准确范围

已对照原文核实（`CloudDRTS/phase1_sensitivity.py`，结果在 `results/phase1_sensitivity.json`）：表 2 的 0.481 Hz、0.023 是手动励磁、1 回联络线、零传输工况，与 v7 的单回联络线一致；原文没有给出本地模态数值。按附录 A 加入快速静止励磁（KA=200、TR=0.01 s、无 TGR）后计算原文表 1、表 2 的四个工况：教材网络（0.15 p.u. 升压变压器）与原文图 1 的无升压变压器网络，把四个已发表的频率和阻尼比全部夹在中间；每区补偿 100–350 Mvar 只使频率在约 0.415–0.441 Hz 内变化，四阶与六阶只差约 0.001 Hz。差距主要来自原文未给出的机组接入阻抗；没有拟合或采用任何电抗值。

已对照 Kundur 教材例 12.6（第 813–816 页）核实：系统、负荷、补偿和运行点数据与代码一致；表 E12.3 的 0.545/1.087/1.117 Hz、阻尼比 0.032/0.072/0.072 及本地模态的区域归属与代码参考值一致。例题计入饱和，按其饱和数据计算，频率差在约 1.1% 以内，表 E12.4(ii) 晶闸管励磁的三个频率差在约 0.7% 以内，且同样得到不稳定的区间模态；本地模态阻尼比高约 1.3–1.4 个百分点的原因已查明：教材采用经典参数换算的等效电路机组模型（不计饱和时与本仓库的 GENROU 模型特征值一致），并在小干扰分析中用增量饱和系数 1/(1+Bsat·ψI) 线性化两轴饱和；GENROU 的割线饱和与之不同。`CloudDRTS/textbook_reproduction.py` 按教材方法复现了表 E12.3 全部 24 个特征值（最大差 0.012 s⁻¹，来自书中虚部 6.82 的印刷舍入）和表 E12.4 的 (ii)、(iii) 两种励磁（差 ≤0.004 s⁻¹）；(iv) 加 PSS 时最多差约 0.10 s⁻¹，原因未明；(i) 直流励磁未建模。ANDES 与本仓库模型没有错误。不计饱和的基本工况有一个 +0.0049 s⁻¹ 的缓慢实特征值（ANDES 同样给出），计入饱和后消失；`results/detailed_checks.json` 和 `external_run/results/modal_comparison_checks.json` 记录全部正实部特征值。


原论文表 2 的零传输、手动励磁区域间频率为 0.481 Hz、阻尼比 0.023；现有对称暂定配置约为 0.425546 Hz、4.475%，频率差约 -11.53%。原文使用恒阻抗负荷，已有对称测试已作该设置。补偿和磁路/饱和的完整对应尚未唯一确定，不能以临时开放基准替代对原论文的解释。

官方 ANDES 网页属于 1.9.3 文档，但打印输出来自 1.7.5 时代；同一输入跨版本会变化。默认深故障流程在清除后失败，只有明示的事件代数一致性适配运行完成。报告和原始失败记录均保留，没有关闭稳定性检查或改故障参数。

## 共享代码

目标 GitHub 仓库为 https://github.com/HardyHaoranLiu/CloudDRTS ，克隆与复跑方法见 `Git_Reproduction.md`。完整程序包包含运行结果及说明；其中的 Git bundle 保留 v7 源码提交和输入来源，可作为离线副本。
