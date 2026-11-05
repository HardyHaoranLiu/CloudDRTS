"""Create phase-specific notes and an evidence-based acceptance matrix."""
from pathlib import Path
import csv,json,re,hashlib
import pymupdf as fitz
from matplotlib.font_manager import findfont

ROOT=Path(__file__).resolve().parent;RES=ROOT/'CloudDRTS/results';EXT=ROOT/'external_run'
OUT=ROOT/'output/pdf';OUT.mkdir(parents=True,exist_ok=True)
def rows(path):return list(csv.DictReader(path.open()))
def info(path):return json.loads(path.read_text())
css='*{font-family:sans-serif}table{border-collapse:collapse;font-size:10.3px;width:100%}th,td{padding:5px;border-bottom:1px solid #ccdcd9;text-align:left;vertical-align:top}th{background:#e7efed}'

class Note:
    def __init__(self):self.doc=fitz.open()
    def text(self,s,size=11.1,gap=7):
        mapping=str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺','0123456789-+')
        s=re.sub('[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+',lambda m:'<sup>'+m.group().translate(mapping)+'</sup>',s)
        r=fitz.Rect(43,self.y,552,749)
        spare,_=self.page.insert_htmlbox(r,f'<div style="font-size:{size}px;line-height:1.4;color:#203638">{s}</div>',css=css,scale_low=1)
        if spare<0:raise ValueError(f'Page {len(self.doc)} overflow: {s[:60]}')
        self.y+=r.height-spare+gap
    def start(self,title,subtitle):
        self.page=self.doc.new_page(width=595,height=810);self.y=43
        self.text('<b>'+title+'</b>',22,6);self.text(subtitle,12,12)
    def heading(self,s):self.text('<b>'+s+'</b>',13.5,4)
    def table(self,head,data):
        self.text('<table><tr>'+''.join('<th>'+s+'</th>' for s in head)+'</tr>'+''.join('<tr>'+''.join('<td>'+str(v)+'</td>' for v in row)+'</tr>' for row in data)+'</table>',10.4,9)
    def graph(self,path,height):
        with fitz.open(path) as source:self.page.show_pdf_page(fitz.Rect(43,self.y,552,self.y+height),source,0)
        self.y+=height+6
    def note(self,s):self.text(s,9.2,7)
    def save(self,name,count):
        assert len(self.doc)==count
        for i,p in enumerate(self.doc,1):
            p.insert_font(fontname=f'Page{i}',fontfile=str(findfont('DejaVu Sans')))
            p.insert_text((531,785),str(i),fontname=f'Page{i}',fontsize=9)
        target=OUT/name;self.doc.save(target,garbage=4,deflate=True)
        qa=ROOT/'qa'/target.stem;qa.mkdir(parents=True,exist_ok=True)
        with fitz.open(target) as reopened:
            for i,p in enumerate(reopened,1):p.get_pixmap(matrix=fitz.Matrix(1.2,1.2)).save(qa/f'{i}.png')
        print(target)

def main():
    dynamic=info(RES/'dynamic_thevenin_checks.json')
    event=info(RES/'event_validation.json');audit=info(RES/'audit_checks.json')
    timing=rows(RES/'timing.csv');grid=rows(RES/'grid_cut_metrics.csv')
    reduced=rows(RES/'dynamic_thevenin_metrics.csv')

    # Phase 1 notes (Chinese and English) are built by phase1_note/build_phase1_note.py.
    n=Note();n.start('阶段三｜接口、延迟与等效','三页实验记录：同一个 Bergeron 小电路、相同初态和投切事件，分别比较历史项交换、静态戴维南近似和保留电容状态的动态等效。')
    n.heading('边界合同与两种戴维南模型')
    n.text('两端电流均定义为流入线路：i = v/Z + h；每端发布的出射历史为 −2v/Z − h。只交换线路历史，不传送内部状态。Z=50 Ω，物理传播时间 250 µs；每个方向另加 kΔt 的人工延迟，双向往返附加量为 2kΔt。')
    n.table(['左侧子系统','等效定义','保留的物理状态'],[
        ['25 Ω 源阻抗，50 Ω 并联负荷，10 µF 电容','G=1/25+1/50=0.06 S','电容电压'],
        ['静态近似','E=2Vs/3，R=1/G=16.667 Ω','省略电容，只匹配直流'],
        ['动态戴维南','Z(s)=1/(G+sC)，E(s)=(Vs/Rs)/(G+sC)','保留受端加载下的电容历史']])
    n.graph(RES/'dynamic_thevenin_delay.pdf',177)
    n.note('图 3　末十个周期的电压 NRMSE；动态等效与完整历史项模型重合。静态近似可能因延迟相移出现误差抵消。')
    n.heading('新增结论：零延迟模型误差不可算作通信误差')
    n.table(['频率','静态零延迟电压误差','动态零延迟电压误差'],[[f'{f} Hz']+[f'{float(next(r["total_v_pct"] for r in reduced if r["frequency_Hz"]==str(f) and r["added_delay_us"]=="0" and r["variant"]==variant and r["window"]=="steady_last_10_cycles")):.4g}%' for variant in ['static_thevenin','dynamic_thevenin']] for f in [60,400,1200]])
    n.note('误差分别记录：对未切分完整模型的总偏差、各变体自身零延迟模型的附加延迟偏差、静态模型本身的零延迟偏差。三者是不同的 NRMSE 比较，不能直接按百分数相加。')

    n.start('阶段三｜动态等效的独立核验','新增求解器使用独立标量公式与独立线路历史队列，不调用原 Circuit 或 Bergeron 类；另用频域两端口方程进行解析交叉检查。')
    n.heading('梯形伴随等效保留加载状态')
    n.text('每步令 Gc=2C/Δt，Rcomp=1/(G+Gc)，Ecomp=(Vs/Rs−Hc)/(G+Gc)。边界满足 v=Ecomp−Rcomp·i。求出实际受载电压后，iC=Gc·v+Hc，并更新 Hc,new=−(Gc·v+iC)。此处的历史源来自受载电容状态，不是把空载电压滤波后直接串联一个静态电阻。')
    n.graph(RES/'dynamic_thevenin_validation.pdf',182)
    n.note('图 4　左：固定物理系统的独立频域解步长误差；右：400 Hz 零延迟下静态近似与动态等效的模型误差。')
    n.table(['验证','实际结果'],[
        ['3 频率 × 5 延迟，与完整离散子电路比较',f'最大电压差 {dynamic["exact_discrete_equivalent_max_voltage_difference_V"]:.3e} V；电流差 {dynamic["exact_discrete_equivalent_max_current_difference_A"]:.3e} A'],
        ['每步独立 KCL 检查',f'最大残差 {dynamic["max_kcl_residual_A"]:.3e} A'],
        ['稳态波形拟合 vs 独立离散频域解',f'最大电压差 {dynamic["max_steady_fit_vs_independent_discrete_voltage_V"]:.3e} V'],
        ['50、25、10、5、2.5 µs，连续频域目标',f'电压相量误差收敛阶 {dynamic["continuous_phasor_observed_order"]:.6f}']])
    n.text('频域核对使用线路两端口导纳和 RC 连续阻抗；离散目标仅把电容中的角频率替换为梯形积分的 2tan(ωΔt/2)/Δt。物理传播时间及附加 100 µs 延迟在各步长中保持不变，避免改变模型后再比较收敛。电压与电流误差在 CSV 中按各自单位分别保存。')
    n.note('复跑：python CloudDRTS/dynamic_thevenin.py；完整入口 run_all.py 共十三个实验/验证脚本。该解析式适用于本次线性 RC 单端口，不宣称已完成一般动态电网的高阶降阶。')

    n.start('阶段三｜推广到 Kundur 的方案','已有十一母线派生 EMT 网络提供两处线路切分与延迟扫描。本轮新增动态 RC 等效作可验证起点；完整机组 EMT 和多端口等效属于后续扩展。')
    n.heading('已有网格试验与适用边界')
    n.table(['100 µs 单向附加延迟、全窗','母线 9 电压 NRMSE','受端线路电流 NRMSE'],[[cut]+[f'{float(next(r[key] for r in grid if r["cut"]==cut and r["delay_us"]=="100" and r["window"]=="all")):.3f}%' for key in ['v9_rms_pct','i9_rms_pct']] for cut in ['7-8','8-9']])
    n.text('当前派生网络使用固定正弦内电势与平衡单相瞬时模型，长线无损且传播时间固定为 390 µs。独立潮流已重新求解并通过功率平衡；它与机电相量模型具有物理线路差异。网格每个延迟工况从各自稳态初始化，总误差含稳态偏移；该口径不能混同于从同一历史状态突然开启通信延迟。')
    n.heading('下一步试验合同')
    n.table(['步骤','固定条件与输出'],[
        ['1　统一模型和初始化','先固定同一物理 EMT 网格、元件、负荷、源、故障时刻；区分各自稳态与共同历史两种初始化。'],
        ['2　构造边界等效','先做线性子网的动态多端口消元，再研究机组状态耦合；核对稳定性、被动性和初态映射。'],
        ['3　切分与误差扫描','分别切 7–8、8–9；保留零延迟、相同步长和相同事件基线，分故障前、故障中、清除后统计电压、电流及边界功率。'],
        ['4　接入实际运行环境','再加入非整数延迟、抖动和丢包，统计时间戳、每步尾部耗时及超时；与确定性整数延迟结果分开报告。']])
    n.heading('入门任务与后续研究的分界')
    n.text('小电路切分、历史项交换、整数步延迟、第二种等效及误差曲线已完成；本轮补上精确保留电容动态的等效对照。Kundur 推广方案和接口符号在本说明及 BOUNDARY_CONTRACT.json 中给出。完整动态 EMT、真实通信、FPGA 和 HIL 尚未由这些结果验证，亦非阶段三入门练习的全部必做项。')
    n.note('阅读对应：原计划 [6]、[7] 提供地理分布式仿真背景，[8]、[9] 对应接口问题，[10] 为频率相关等效背景。此处没有把尚未运行的抖动、丢包或硬件试验写成已完成。')
    n.save('Phase3_Partitioning_Note.pdf',3)

    # Machine-readable scope and the requested first-week discussion list.
    perf={}
    for budget in [5,50]:
        rr=[r for r in timing if float(r['budget_us'])==budget]
        perf[str(budget)]=dict(median_us_range=[min(float(r['median_us']) for r in rr),max(float(r['median_us']) for r in rr)],
            miss_pct_range=[min(float(r['miss_pct']) for r in rr),max(float(r['miss_pct']) for r in rr)],
            max_us=max(float(r['max_us']) for r in rr),runs=len(rr),steps_each=4000)
    matrix=dict(
        phase1=dict(status='partial_reference_alignment',practice='fault, linearization, computed mode labels, external matched-case check, sensitivity sweep and two-page notes (zh, en) complete',
          notes=['output/pdf/Phase1_Note_EN.pdf','output/pdf/Phase1_Note_ZH.pdf'],
          unresolved='original-paper operating point and textbook damping difference; see CloudDRTS/results/phase1_sensitivity.json'),
        phase2=dict(status='practice_checks_passed',checks=['RC/RL/RLC analytic validation','50 and 5 us plus refinements',
          'Bergeron reflection','switch event state continuity','LC energy','eight invalid input regressions'],
          emt_external_software='not performed; not required by onboarding plan'),
        phase3=dict(status='practice_checks_passed',dynamic_thevenin=dynamic,
          unresolved='dynamic multiport Kundur EMT coupling and real communications are future work'),
        phase4=dict(status='python_profiling_path_complete',performance=perf,
          real_time_guarantee=False,hls='alternative path not selected; not a missing mandatory step'),
        shared_code='Portable Git bundle provided; no remote repository has been configured or pushed.',
        reading_and_reviews='Reading map supplied; no claim that the user completed all readings or mentor reviews.',
        external_scope='Prior ANDES raw runs retained; external result analysis rechecked against this run local fault CSV.')
    (ROOT/'ACCEPTANCE_STATUS.json').write_text(json.dumps(matrix,indent=2,ensure_ascii=False))
    contract=dict(units=dict(voltage='V',current='A',time='s'),direction='current positive into line at each port',
        incoming='h',port_relation='i=v/Z+h',outgoing='-2*v/Z-h',
        exchange_order=['read incoming histories','solve local subsystem','publish outgoing history'],
        payload='one scalar outgoing history per direction per step',
        physical_delay_s=.00025,artificial_delay_s='k*dt each direction',round_trip_added_s='2*k*dt',
        initialization='zero histories and zero capacitor state for RC toy circuit',
        network_implementation='in-process deterministic queue; no actual network transport',
        recommended_future_metadata=['sample time','sequence number','declared validity interval'])
    (ROOT/'BOUNDARY_CONTRACT.json').write_text(json.dumps(contract,indent=2))
    text=f'''# 阶段四：实际剖析与第一周讨论清单

本阶段按原计划选择 Python 每步剖析路径，HLS 入门是另一条可选路径。没有把未上 FPGA 视为本阶段未完成，也没有把离线脚本当作实时原型。

## 本次测量

- 5 µs 预算：5 轮，每轮 4000 步；单步中位数范围 {perf['5']['median_us_range'][0]:.3f}–{perf['5']['median_us_range'][1]:.3f} µs；超时比例 {perf['5']['miss_pct_range'][0]:.3f}–{perf['5']['miss_pct_range'][1]:.3f}%。
- 50 µs 预算：同样 5 轮；中位数范围 {perf['50']['median_us_range'][0]:.3f}–{perf['50']['median_us_range'][1]:.3f} µs；超时比例 {perf['50']['miss_pct_range'][0]:.3f}–{perf['50']['miss_pct_range'][1]:.3f}%。
- 实测最大耗时分别为 {perf['5']['max_us']:.3f} / {perf['50']['max_us']:.3f} µs。它们是有限样本的最大值，不是最坏执行时间上界。
- 原始计时、分位数、运行环境、调用热点分别在 CloudDRTS/results/timing_raw.csv、timing.csv、environment.json、profile.txt。计时包括模型计算，不代表完整网络 I/O、记录和实时调度预算。

## 第一周讨论问题

1. 目标平台、时钟及采样步长是什么？本地求解、线路收发、同步和记录各分配多少预算？
2. 固定拓扑能否缓存分解或标量求解系数？遇到未预计算的拓扑怎样处置并保持明确的时限？
3. Python 的 NumPy 小矩阵调用与状态更新各占多少耗时？应先转为标量、编译内核还是直接 HLS？
4. 用什么计时与负载条件区分求解耗时、主机调度抖动和真实网络抖动？超时一次时采用什么保持或丢弃策略？
5. 边界历史项用何种定点字长、缩放与饱和策略？怎样核对符号、时间戳、初态和能量交换？
6. 电压/电流 NRMSE、峰值误差、故障清除时刻和每步超时比例，各自的验收阈值由谁确定？
7. 动态戴维南单端口推广为 Kundur 多端口时，哪些机组、控制与负荷状态必须保留？怎样验证被动性和稳定性？
8. KAUST 首轮验证用哪套 RTDS、OPAL-RT、Typhoon HIL 或 FPGA 接口？哪些对照算例和波形能够作为共享参考？

## 当前结论

剖析与问题清单已完成。现有解释执行的 Python 求解器没有硬实时保证；后续优化和上板是项目工作，不能由本轮均值或中位数替代验收。
'''
    (ROOT/'Phase4_Questions.md').write_text(text)
    print('Acceptance matrix and phase 4 questions written.')

if __name__=='__main__':main()
