"""Analytic RLC convergence, fault-duration sweep and excitation-frequency sensitivity."""
from common import np, plt, figure, table, OUT
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import kundur as kd
from emt import Circuit
from experiments import benchmark

def rlc():
    R,L,C=6.,.01,10e-6
    alpha=R/(2*L); wd=np.sqrt(1/(L*C)-alpha**2)
    fig,ax=plt.subplots(1,2,figsize=(7.15,2.9))
    rows=[]; curves=[]
    for us in [100,50,25,10,5]:
        dt=us*1e-6; t=np.arange(round(.01/dt)+1)*dt
        cir=Circuit(2,dt)
        cir.resistor(-1,1,R);cir.reactive('L',1,2,L,v0=1)
        cir.reactive('C',2,0,C)
        v=np.r_[0,[cir.step({-1:1.},[0.,0.])[1] for _ in t[1:]]]
        a=1-np.exp(-alpha*t)*(np.cos(wd*t)+alpha/wd*np.sin(wd*t))
        err=max(abs(v-a)); rows.append((us,err,np.sqrt(np.mean((v-a)**2)),v.max()))
        curves.extend(zip(np.full(t.shape,us),t,v,a))
        if us in (100,5):ax[0].plot(t*1000,v,label=f'{us} µs',ls='--' if us==100 else '-')
    ax[0].plot(t[::24]*1000,a[::24],'.',color='black',label='Analytic',ms=2.5)
    ax[0].legend(); ax[0].set(xlabel='Time (ms)',ylabel='Capacitor voltage (V)',title='Series RLC step')
    steps,err=np.array(rows)[:,:2].T
    slope=np.polyfit(np.log(steps),np.log(err),1)[0]
    assert 1.95<slope<2.05
    ax[1].loglog(steps,err,'o-',label=f'Fitted order {slope:.2f}')
    ax[1].legend();ax[1].set(xlabel='Time step (µs)',ylabel='Maximum error (V)',title='Time-step refinement')
    table('rlc_convergence',['dt_us','max_error_V','rms_error_V','v_peak_V'],rows)
    table('rlc_waveforms',['dt_us','time_s','voltage_V','analytic_V'],curves)
    figure(fig,'rlc_convergence')

def clearing():
    fig,ax=plt.subplots(1,2,figsize=(7.15,2.9))
    rows=[];waves=[]
    for ms in [40,80,120,160]:
        tf=ms*.001
        def run(z,x,end):
            ts=np.linspace(0,end,round(end/.005)+1)
            ans=solve_ivp(lambda t,x:kd.rhs(x,z),(0,end),x,t_eval=ts,method='DOP853',rtol=1e-9,atol=1e-11,max_step=.03)
            assert ans.success
            return ans
        fault=run(kd.short,kd.x0,tf)
        rec=run(kd.normal,fault.y[:,-1],8.)
        angles=np.rad2deg(rec.y[:4]-kd.delta[:,None])
        gap=angles[:2].mean(axis=0)-angles[2:].mean(axis=0)
        peak=np.max(abs(gap)); peaks=find_peaks(gap,prominence=.05)[0]
        assert len(peaks)>=3
        freq=1/np.median(np.diff(rec.t[peaks]))
        rows.append([ms,peak,freq]); waves.extend(zip(np.full(len(gap),ms),rec.t,gap))
        ax[0].plot(rec.t,gap,label=f'{ms} ms')
    ax[0].legend(ncol=2);ax[0].set(xlabel='Time after clearing (s)',ylabel='Area angle difference (deg)')
    ax[1].plot(np.array(rows)[:,0],np.array(rows)[:,1],'o-')
    ax[1].set(xlabel='Fault duration (ms)',ylabel='Peak area separation change (deg)')
    figure(fig,'clearing');table('clearing',['duration_ms','peak_delta_deg','peak_frequency_Hz'],rows)
    table('clearing_traces',['duration_ms','time_after_clear_s','area_delta_deg'],waves)

def spectrum():
    rows=[]
    freqs=[60,200,400,800,1200]
    delays=[0,25,50,100,200,400]
    fig,axes=plt.subplots(1,2,figsize=(7.15,2.85),sharey=True)
    for i,f in enumerate(freqs):
        dt=5e-6
        _,ref,_=benchmark(dt,frequency=f)
        for ax,mode in zip(axes,['split','thevenin']):
            values=[]
            for us in delays:
                _,y,_=benchmark(dt,us//5,mode,frequency=f)
                dif=y[:,1:]-ref[:,1:]
                pct=100*np.sqrt(np.sum(dif**2,axis=0)/np.sum(ref[:,1:]**2,axis=0))
                rows.append([f,us,mode,*pct]);values.append(pct[0])
            ax.plot(delays,values,marker=['o','s','^','v','D'][i],ls=['-','--','-.',':','-'][i],label=f'{f} Hz')
    axes[0].set_ylabel('Voltage relative RMS error (%)')
    for ax,title in zip(axes,['Full split','Static Thévenin']):
        ax.set(xlabel='Added one-way delay (µs)',title=title)
    axes[0].legend(ncol=2);figure(fig,'frequency_delay')
    table('frequency_delay',['source_Hz','delay_us','variant','v_rms_pct','i_rms_pct'],rows)

def timing_plot():
    a=np.genfromtxt(OUT/'timing_raw.csv',delimiter=',',names=True)
    fig,axs=plt.subplots(1,2,figsize=(7.15,2.65),sharey=True)
    for ax,us in zip(axs,[5,50]):
        for rep in range(1,6):
            x=np.sort(a['elapsed_ns'][(a['dt_us']==us)&(a['rep']==rep)]/1000)
            ax.semilogx(x,(np.arange(len(x))+1)/len(x)*100,alpha=.7)
        ax.axvline(us,color='black',ls='--',label='Step budget');ax.legend(loc='lower right')
        ax.set(xlabel='Measured step time (µs)',title=f'{us} µs budget')
    axs[0].set_ylabel('Cumulative share (%)');figure(fig,'timing_cdf')
if __name__=='__main__':
    rlc(); clearing(); spectrum(); timing_plot()
