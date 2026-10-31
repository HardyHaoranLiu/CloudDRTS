"""Build the Chinese audit report from the actual CSV/JSON run evidence."""
from pathlib import Path
import csv,json,html,re
import pymupdf as fitz
from matplotlib.font_manager import findfont

ROOT=Path(__file__).resolve().parent
RES=ROOT/'results';FIG=ROOT/'figures';OUT=ROOT/'output/pdf'
OUT.mkdir(parents=True,exist_ok=True)
def rows(name):return list(csv.DictReader((RES/name).open()))
def info(name):return json.loads((RES/name).read_text())

doc=fitz.open()
css='*{font-family:sans-serif}p{margin:0}table{border-collapse:collapse;width:100%;font-size:10.4px}th,td{padding:6px 5px;border-bottom:1px solid #d2dfdd;text-align:left;vertical-align:top}th{background:#e7efed}'
def text(s,size=11.3,gap=7):
    global y
    superscript=str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺','0123456789-+')
    s=re.sub('[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+',lambda m:'<sup>'+m.group().translate(superscript)+'</sup>',s)
    rect=fitz.Rect(43,y,552,752)
    spare,_=page.insert_htmlbox(rect,f'<div style="font-size:{size}px;line-height:1.4;color:#203638">{s}</div>',css=css,scale_low=1)
    if spare<0:raise ValueError(f'Page {len(doc)} overflow at y={y}: {s[:50]}')
    y+=rect.height-spare+gap
def start(title,sub):
    global y,page
    page=doc.new_page(width=595,height=810);y=35
    y=62
    text('<b>'+title+'</b>',23,5);text(sub,12.3,12)

def heading(s):text('<b>'+s+'</b>',14,4)
def table(head,body):
    text('<table><tr>'+''.join('<th>'+s+'</th>' for s in head)+'</tr>'+''.join('<tr>'+''.join('<td>'+str(v)+'</td>' for v in r)+'</tr>' for r in body)+'</table>',10.6,10)
def graph(name,height):
    global y
    with fitz.open(FIG/(name+'.pdf')) as figure:page.show_pdf_page(fitz.Rect(43,y,552,y+height),figure,0)
    y+=height+7
def note(s):text('<span style="color:#4d696b">'+s+'</span>',9.3,7)

summary=info('verification_summary.json')
modes=rows('matched_mode_comparison.csv')
official=rows('official_version_comparison.csv')
metric=rows('fault_waveform_comparison.csv')
coarse=info('fault_consistent_run.json');fine=info('fault_fine_run.json')
def metric_value(run,name):return next(float(r['max_abs_error']) for r in metric if r['run']==run and r['quantity']==name)

start('外部基准已实际运行','本次使用独立 ANDES 引擎复现公开示例，再把同参数十一母线模型与现有六阶程序比较。公开值、当前软件值和本程序值分别记录。')
table(['核对项目','本轮实际结果','结论'],[
    ['公开历史模态','ANDES 1.7.5 三对机电模态的实部、虚部、频率与阻尼均落在公开表的舍入区间内','公开精度复现通过'],
    ['版本敏感性','同一 XLSX 文件在 1.9.3 中的区域间频率为 0.646897 Hz；1.7.5 为 0.672389 Hz','不能把新版输出直接替代历史值'],
    ['十一母线同参数模态','三模态最大频差 6.55×10⁻⁷ Hz；最大复特征值差 4.13×10⁻⁶ s⁻¹','现有六阶实现获独立引擎支持'],
    ['母线 8 三相故障','默认 ANDES 流程失败；明示事件一致性适配后，1 ms 和 0.5 ms 两次均完成 12 s','适配后波形可核对，默认流程未通过']])
heading('先明确比较的是哪两个案例')
table(['配置','用途'],[
    ['官方 kundur_full：10 母线，GENROU + EXDC2 + TGOV1','复现 ANDES 网页示例及检查版本变化'],
    ['本次匹配案例：11 母线，六阶无饱和、固定励磁和机械转矩','核对现有自编模型；不与上行模态强行逐项等同']])
text('本轮补足了实验 1 的独立相量软件验证，但尚未把 1991 年论文的原始补偿、饱和与控制配置全部还原。因此不能据此宣布四项实验及原论文所有结果均已精确复现。')
note('执行环境：Python 3.12.14；NumPy 1.26.4；SciPy 1.13.1；ANDES 1.9.3 安装版及 1.7.5 官方发布源码。完整依赖锁和原始日志随包提供。')

start('公开历史输出：已复现','网页属于 ANDES 1.9.3 文档，但其中打印的示例报告来自 2022-07-05 的 1.7.5.post2.dev0+ga53a31fa。此处使用正式发布版 1.7.5 核对公开精度。')
graph('official_versions',190)
note('图 1　两套 ANDES 版本使用同一官方输入文件。圆点为本地历史版本实跑，方点为 1.9.3 实跑，黑叉为网页公开值。')
table(['模态','公开频率 / Hz','1.7.5 实跑 / Hz','1.9.3 实跑 / Hz'],[
    [r['mode'],r['published_frequency_Hz'],f'{float(r["andes_175_frequency_Hz"]):.9f}',f'{float(r["andes_193_frequency_Hz"]):.9f}'] for r in official])
text('历史版本的实部、虚部、阻尼和频率共 12 个数值，全部通过各自最后一位的半单位舍入检查。该结论仅指这三对机电模态的公开精度，不等同于证明整个开发版与发布版逐状态完全相同。')
heading('发现的问题：相同输入不保证跨版本同结果')
text('两个版本自带的 kundur_full.xlsx 的 SHA-256 完全相同。区域间频率变化约 −3.79%，阻尼由 4.53457% 变为 3.43092%。差异与软件版本及其默认实现行为相关；本轮尚未把变化定位到某一个具体提交或模型公式。')
note('XLSX SHA-256：f725e03ba12d8207616f68acdd606bbd35e7c4a68f13e66d7db43925adac2ed8。正式复跑应固定版本，不能只固定案例文件。')
note('公开来源：[1] <a href="https://docs.andes.app/en/v1.9.3/_examples/ex4.html">ANDES Eigenvalue Analysis</a>；历史源码：[2] <a href="https://pypi.org/project/andes/1.7.5/">PyPI ANDES 1.7.5</a>。全部实跑数值来自 results/official_version_comparison.csv。')

start('现有六阶模型：同参数核对通过','重新构造 100 MVA、60 Hz 的十一母线网络。四台 900 MVA 机组保留定子电阻和次暂态状态，动态负荷为有功恒电流、无功恒阻抗。')
graph('matched_modes',179)
note('图 2　左：两套独立引擎的机电频率重合。右：绝对频差，单位为微赫兹。')
table(['模态','自编 / ANDES 频率 Hz','自编 / ANDES 阻尼 %'],[
    [r['mode'],f'{float(r["local_frequency_Hz"]):.9f} / {float(r["andes_frequency_Hz"]):.9f}',f'{float(r["local_damping_pct"]):.6f} / {float(r["andes_damping_pct"]):.6f}'] for r in modes])
heading('一致性与数值差异')
text('三对特征值最大差 4.13×10⁻⁶ s⁻¹，满足本次数值核对阈值 10⁻⁴ s⁻¹；最大频差 6.55×10⁻⁷ Hz，满足 10⁻⁵ Hz。11 个母线的最大电压幅值差为 6.53×10⁻⁷ p.u.。潮流、特征值和匹配规则均保存在 CSV / JSON 中。')
text('ANDES 1.9.3 在线路导纳公式中给 r、x 各加 10⁻⁸ p.u. 数值正则项；本程序沿用原始线路参数，因此两者并非逐比特相同的有效网络。未通过调参拟合输出。自动匹配使用复特征值距离的匈牙利算法。')
heading('本次匹配条件的边界')
text('无饱和、无励磁调节器、无调速器；励磁电压和机械转矩固定。惯量 H 为 6.5 / 6.175 s，定子电阻 0.0025 p.u.，d/q 轴次暂态电抗均为 0.25 p.u.。这些选择用于核对现有程序，不代表已验证带 PSS 或原论文的饱和模型。')
note('近零参考特征值出现约 ±1.92×10⁻⁶ s⁻¹ 的数值分裂。它们不作为机电模态误差或失稳结论。参数背景：[3] <a href="https://colib.net/testCases/kundur_two_area_system/">CoLib Kundur 案例</a>。')

start('三相故障：波形完成独立比较','母线 8 在 1.00 s 投入、1.08 s 清除电阻 10⁻⁴ p.u. 的对称故障，仿真至 12 s。自编 DOP853 结果与 ANDES 隐式梯形结果比较。')
graph('fault_overlay',280)
note('图 3　上：自编连续积分与 ANDES 0.5 ms 的转角及联络功率叠加。下：ANDES 减自编程序的误差。两次事件的右极限均参与统计。')
table(['最大绝对差','ANDES 1 ms','ANDES 0.5 ms'],[
    [label,f'{metric_value("consistent",key):.4g} {unit}',f'{metric_value("fine",key):.4g} {unit}']
    for key,label,unit in [('COI_angle','COI 相对转角','°'),('speed_deviation','各机频差','Hz'),('bus8_voltage','母线 8 电压','p.u.'),('tie_power','7–8 联络有功','MW')]])
text('按自编结果的 2401 个输出时刻比较，未剔除采样点。ANDES 原始时间序列在事件时刻记录回调前的值；两处精确事件行改用其独立求解的右极限，其余行线性插值。原始序列与右极限值均保留。联络有功统一按输入 r/x 与两端复电压计算；定义及事件逐项值见 README 和 event_side_alignment.csv。')
note('步长减半降低了转角、电压和功率的最大差异，但频差的最大误差未单调下降。不能据两次结果宣称严格二阶跨软件收敛；模型正则化、求解与插值误差仍混合存在。自编波形为 v5 已执行的 DOP853 数据，输出间隔 5 ms，rtol=10⁻⁸、atol=10⁻¹⁰。')

start('故障流程的失败与修正：完整披露','默认 ANDES 深故障流程在清除后保持近零电压支路，并在约 4.097 s 触发转角稳定性停止。该失败记录保留，未关闭稳定性判据。')
graph('fault_event_audit',180)
note('图 4　左：默认流程与一致性事件处理的电压差异；右：适配后 0.5 ms 波形相对自编程序的电压差。默认失败波形不能作为正确物理解。')
heading('事件适配具体做了什么')
text('仍调用 ANDES 原始投入 / 清除回调。在事件瞬间固定全部微分状态，使用 ANDES 自身的网络、机器状态和负荷系数建立直角坐标电流平衡方程，求解新的代数工作点，再用未修改的 ANDES 全部代数残差 g 复核。未读取自编波形来构造该解。')
table(['0.5 ms 运行事件','电流残差','ANDES 全 g 残差','状态跳变'],[
    [f'{r["t"]:.2f} s',f'{r["current_balance_residual"]:.3e}',f'{r["andes_algebraic_residual"]:.3e}',str(r['max_state_change'])] for r in fine['events']])
text('适配器不改变连续模型方程、梯形积分器、故障电阻、事件时刻或保护阈值；它改变了事件后的代数初始化。适用范围被限制为本次固定励磁、无饱和、等次暂态电抗、有功 I / 无功 Z 负荷案例，不应直接套用于任意 ANDES 模型。')
heading('本轮仍未闭合的证据')
text('原生深故障流程尚未找到无需适配的可靠配置；本轮没有完成上游缺陷定责。原论文对称单回、零传输案例仍存在参数与模型映射差异；此前的频率偏差和阻尼偏差不能由这次同参数核对自动消除。STEPSS 的公开 PSS 结果仍仅作补充资料，本轮没有运行 STEPSS。')
note('为受限运行环境中缺失的 /proc/cpuinfo，启动脚本将 psutil.cpu_count 固定返回 1，仅用于串行调度。ANDES 供应方方程代码未修改。')

start('程序包、复跑与四项实验状态','本包将实际日志、输入、参考精度、完整波形和失败尝试一起保留。每个结论都可回溯到对应的 CSV、JSON 或原始 ANDES 报告。')
table(['文件 / 入口','用途'],[
    ['external_run/run_reproduction.py','按顺序运行历史 EIG、新版 EIG、同参数模态、两步长故障和结果统计'],
    ['external_run/requirements_lock.txt','本次实际安装依赖的完整冻结版本；Python 3.12'],
    ['external_run/results/','原始报告、全部波形、事件右极限、数值比较与核对结论'],
    ['external_run/historical/andes-1.7.5.tar.gz','官方历史源码；附 SHA-256 与 GPL 许可证'],
    ['CloudDRTS/','自编程序源码与本次使用的 v5 原始故障结果']])
heading('复跑命令')
text('<span style="font-family:monospace;font-size:10px">python3.12 -m venv .venv<br/>.venv/bin/python -m pip install -r external_run/requirements_lock.txt<br/>.venv/bin/python external_run/run_reproduction.py<br/># 需要同时重算自编故障波形时追加：<br/># --rerun-local-transient</span>',10.8,9)
note('Windows 使用 .venv\\Scripts\\python.exe。历史版旧打包脚本不兼容 Python 3.12，因此直接使用经哈希验证的官方源码执行；模型代码未修改。首次安装需要可访问 Python 包源。')
heading('本次结论在四项任务中的位置')
table(['任务','本轮之后可作出的结论'],[
    ['1　机电相量实验','临时外部基准与同参数六阶实现核对完成；精确原论文复现仍未完成'],
    ['2　节点 EMT 元件','沿用已有解析及步长验证；本轮没有增加外部 EMT 软件对照'],
    ['3　切分与延迟','沿用已有接口与延迟实验；本轮 ANDES 结果不验证 EMT 分区误差'],
    ['4　剖析 / HLS 入门','已有 Python 性能剖析属于计划中的可选路径；本轮没有新增 HLS 或硬件实测']])
note('来源 [1] ANDES 官方示例；[2] ANDES 1.7.5 发布源码；[3] CoLib 参数说明。完整链接、文件来源、许可证和哈希见 README.md 与 THIRD_PARTY_NOTICES.md。')

header_font=findfont('DejaVu Sans')
for number,p in enumerate(doc,1):
    font_name=f'Header{number}'
    p.insert_font(fontname=font_name,fontfile=header_font)
    p.insert_text((532,785),str(number),fontsize=9,fontname=font_name)

doc.set_metadata({'title':'CloudDRTS 临时外部基准实跑与核对报告','author':'CloudDRTS experiment workspace'})
# Do not re-subset embedded Matplotlib font subsets: it corrupts some plot glyphs.
target=OUT/'CloudDRTS_External_Benchmark_Report.pdf'
doc.save(target,garbage=4,deflate=True)
qa=ROOT/'qa';qa.mkdir(exist_ok=True)
doc=fitz.open(target)
for i,p in enumerate(doc):p.get_pixmap(matrix=fitz.Matrix(1.2,1.2)).save(qa/f'page_{i+1}.png')
print(target)
print('Pages:',len(doc))
