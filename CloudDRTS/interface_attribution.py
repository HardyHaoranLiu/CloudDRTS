"""Compare each model to its own zero-delay baseline, and retain steady windows."""
from common import np, plt, table, figure
from experiments import benchmark

def pct(y,ref):
    return 100*np.sqrt(np.sum((y-ref)**2,axis=0)/np.sum(ref**2,axis=0))

def main():
    rows=[]
    freqs=[60,200,400,800,1200];delays=[0,50,100,200,400]
    for f in freqs:
        duration=max(.06,25/f)
        t,full,_=benchmark(10e-6,duration=duration,frequency=f)
        masks={'full_window':np.ones(len(t),bool),'steady_last_10_cycles':t>=t[-1]-10/f}
        for mode in ['split','thevenin']:
            _,zero,_=benchmark(10e-6,mode=mode,duration=duration,frequency=f)
            for us in delays:
                _,y,_=benchmark(10e-6,us//10,mode,duration=duration,frequency=f)
                for window,m in masks.items():
                    total=pct(y[m,1:],full[m,1:])
                    communication=pct(y[m,1:],zero[m,1:])
                    reduction=pct(zero[m,1:],full[m,1:])
                    rows.append([f,us,mode,window,*total,*communication,*reduction])
    table('interface_attribution',['source_Hz','delay_us','variant','window',
          'total_v_pct','total_i_pct','delay_only_v_pct','delay_only_i_pct',
          'zero_delay_model_v_pct','zero_delay_model_i_pct'],rows)
    fig,axes=plt.subplots(1,2,figsize=(7.15,2.75))
    for mode,marker in [('split','o'),('thevenin','s')]:
        rr=[r for r in rows if r[0]==400 and r[2]==mode and r[3]=='steady_last_10_cycles']
        axes[0].plot(delays,[r[4] for r in rr],marker=marker,label=mode)
        axes[1].plot(delays,[r[6] for r in rr],marker=marker,label=mode)
    for ax in axes:ax.set(xlabel='Added one-way delay (µs)',ylabel='Relative RMS voltage error (%)')
    axes[0].set_title('Against full zero-delay model')
    axes[1].set_title('Against own zero-delay model')
    axes[0].legend();figure(fig,'interface_attribution')
if __name__=='__main__':main()
