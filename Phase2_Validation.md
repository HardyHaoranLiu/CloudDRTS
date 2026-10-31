# 阶段二：EMT 基础验收

本阶段按原计划以已知解析响应和至少两个时间步为主要证据，不把尚未做另一套 EMT 软件验证视为原计划必做项。

| 要求 | 实现/核对 | 证据 |
|---|---|---|
| R、L、C 节点伴随模型 | 梯形 Norton 印章与状态历史更新 | CloudDRTS/emt.py |
| 电压源和开关 | 负节点为指定电位；有限导通电阻、断开时去除印章 | emt.py，experiments.py |
| 无损 Bergeron 线 | 双向出射历史项，物理传播时间为整数步 | emt.py，audit_checks.py |
| 解析校验 | RC/RL 阶跃、RLC 状态空间解析响应、线路传播/反射 | results/validation.csv、rlc_convergence.csv、audit_checks.json |
| 至少两个步长 | 包含 50 µs 和 5 µs，并增加 25/10/2.5 µs 事件细化 | results/event_rlc_convergence.csv、audit_switch_convergence.csv |
| 投切状态连续 | 精确事件时刻保持电容电压与电感电流；事件代数量重新求解 | results/event_validation.json |
| 反射与能量 | 两种步长多次反射检查；无源 LC 储能检查 | results/audit_checks.json |
| 代码和图 | 全部源码、Git bundle、逐点数据、PNG/SVG/PDF | 主程序包与 Git_Reproduction.md |

本版完整入口已重新执行。RLC 投切电流/电压误差约二阶收敛；事件储能状态差为零；八项非法输入回归通过。具体数值以本轮 JSON 为准。

已声明的边界：固定步长；物理线路延迟必须是步长的正整数倍；非理想开关闭合电阻须为正；事件电容约束必须独立、一致。程序是教学用线性节点求解器，不覆盖任意理想源回路、任意拓扑退化或三相不平衡的通用工业 EMT 建模。

下一步若用于论文或完整动态电网，可增加独立 EMT 平台的同模型验证；这属于证据增强，不改变本阶段已有解析检查的结论。
