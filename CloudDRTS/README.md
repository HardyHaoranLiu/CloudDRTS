# CloudDRTS 四阶段实验代码 v7

Python 3.12：

```bash
python -m pip install -r requirements.txt
python run_all.py
```

完整入口按依赖执行十五个脚本，所有输出写入 results/。再次运行覆盖同名生成结果。各组图提供 PNG、可编辑 SVG 和矢量 PDF，保留四个边框。

| 模块 | 作用 |
|---|---|
| kundur.py | 四阶简化相量模型，作为原始基线保留 |
| detailed_machine.py / detailed_transient.py | 六阶模型、三模态（标签由振型计算）、饱和变体、小扰动及三相故障 |
| textbook_reproduction.py | 按教材方法（等效电路机组、增量饱和）复现 Kundur 例 12.6 表 E12.3、E12.4 |
| phase1_sensitivity.py | 复算文献 [2] 表 1、表 2 的四个对称工况（含附录 A 快速励磁），两种网络对照及升压电抗、补偿量扫描 |
| emt.py | R/L/C 节点梯形伴随、指定电位源、有限导通电阻开关、无损 Bergeron |
| experiments.py / waveform.py / extended.py | 解析响应、初始分区/静态近似、时间与频率扫描、逐步性能计时 |
| audit_checks.py / event_validation.py | 输入回归、多次反射、RC/RLC 投切、状态连续性、LC 储能 |
| literature_checks.py | 负荷模型比较及原论文宏观工况的暂定配置 |
| grid_cuts.py / working_point_validation.py | 两处派生 EMT 网格切分、独立物理潮流、功率平衡与步长核对 |
| interface_attribution.py | 区分零延迟模型偏差、附加延迟偏差和总偏差 |
| dynamic_thevenin.py | 独立标量动态等效与独立两端口频域校验 |

四项验收状态、原论文剩余差异、两套独立依赖环境、Git bundle、分阶段说明及第一周问题清单，见包根 README.md 和 ACCEPTANCE_STATUS.json。阶段三说明用包根 build_delivery_notes.py 重建；阶段一中英文说明用 phase1_note/build_phase1_note.py 重建。

当前同参数 ANDES 核对支持六阶无饱和、固定励磁的机电相量模型。派生 EMT 网络仍是固定正弦源和平衡单相瞬时模型，未与动态机组和控制统一耦合。论文对称工况的每区 200 Mvar 补偿为暂定输入，原论文频率及教材阻尼差异尚未闭合。
