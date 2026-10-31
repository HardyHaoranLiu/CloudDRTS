"""Publication-style scientific figures, with all four plot spines retained."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'CloudDRTS'))
from common import np,plt

RES=ROOT/'results';OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
def data(name):return np.genfromtxt(RES/name,delimiter=',',names=True,encoding='utf8')
def save(fig,name):
    for i,ax in enumerate(fig.axes):
        for side in ['left','right','top','bottom']:ax.spines[side].set_visible(True)
        ax.annotate(chr(97+i),(0,1),xycoords='axes fraction',xytext=(-10,6),
                    textcoords='offset points',fontsize=10,fontweight='bold')
    fig.tight_layout()
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{name}.{ext}',bbox_inches='tight',dpi=320)
    plt.close(fig)

def main():
    d=data('official_version_comparison.csv')
    fig,ax=plt.subplots(1,2,figsize=(7.15,2.6))
    x=np.arange(3)
    for a,col,pub in [(ax[0],'frequency_Hz','published_frequency_Hz'),(ax[1],'damping_pct','published_damping_pct')]:
        a.plot(x,d[pub],'kx',ms=7,label='Published 2022 output')
        a.plot(x,d['andes_175_'+col],'o',mfc='none',ms=7,label='ANDES 1.7.5 run')
        a.plot(x,d['andes_193_'+col],'s--',label='ANDES 1.9.3 run')
        a.set_xticks(x,['Inter-area','Local low','Local high'])
    ax[0].set_ylabel('Damped frequency (Hz)');ax[1].set_ylabel('Damping ratio (%)')
    ax[0].legend(loc='upper left',fontsize=6.4)
    save(fig,'official_versions')

    d=data('matched_mode_comparison.csv')
    fig,ax=plt.subplots(1,2,figsize=(7.15,2.45))
    ax[0].plot(x,d['local_frequency_Hz'],'o',ms=7,mfc='none',label='Custom six-state')
    ax[0].plot(x,d['andes_frequency_Hz'],'x',ms=6,label='ANDES 1.9.3')
    ax[0].set(xticks=x,xticklabels=['Inter-area','Local low','Local high'],ylabel='Damped frequency (Hz)')
    ax[0].legend()
    ax[1].bar(x,d['frequency_abs_difference_Hz']*1e6,width=.55)
    ax[1].set(xticks=x,xticklabels=['Inter-area','Local low','Local high'],ylabel='Absolute difference (microHz)')
    save(fig,'matched_modes')

    mine=np.genfromtxt(ROOT.parent/'CloudDRTS/results/detailed_fault.csv',delimiter=',',names=True)
    ext=data('aligned_fault_fine.csv');err=data('fault_error_fine.csv')
    fig,ax=plt.subplots(2,2,figsize=(7.15,4.3),sharex=True)
    for j in range(4):
        line=ax[0,0].plot(mine['t_s'],mine[f'G{j+1}_coi_deg'],label=f'G{j+1}')[0]
        ax[0,0].plot(ext['t_s'],ext[f'G{j+1}_coi_deg'],'--',color=line.get_color())
        ax[1,0].plot(err['t_s'],err[f'G{j+1}_coi_deg']*1e3,color=line.get_color())
    ax[0,0].legend(ncol=4,fontsize=6);ax[0,0].set_ylabel('COI angle change (deg)')
    ax[0,1].plot(mine['t_s'],mine['P7to8_MW'],label='Custom DOP853')
    ax[0,1].plot(ext['t_s'],ext['P7to8_MW'],'--',label='ANDES + event adapter')
    ax[0,1].set_ylabel('Tie-line power (MW)');ax[0,1].legend(fontsize=6)
    ax[1,0].set_ylabel('Angle difference (millideg)')
    ax[1,1].plot(err['t_s'],err['P7to8_MW']*1e3)
    ax[1,1].set_ylabel('Power difference (kW)')
    for a in ax.ravel():a.axvspan(1,1.08,color='grey',alpha=.12)
    for a in ax[1]:a.set_xlabel('Time (s)')
    save(fig,'fault_overlay')

    raw=data('matched_fault_waveform.csv');fine=data('matched_fault_fine_waveform.csv')
    fig,ax=plt.subplots(1,2,figsize=(7.15,2.5))
    ax[0].plot(raw['t_s'],raw['V8_pu'],label='Default ANDES (failed)')
    ax[0].plot(fine['t_s'],fine['V8_pu'],'--',label='With event adapter')
    ax[0].set(xlim=(.8,1.5),xlabel='Time (s)',ylabel='Bus 8 voltage (p.u.)');ax[0].legend(fontsize=6)
    ax[1].plot(err['t_s'],err['V8_pu']*1e6)
    ax[1].set(xlabel='Time (s)',ylabel='Voltage difference (micro p.u.)')
    save(fig,'fault_event_audit')

if __name__=='__main__':main()
