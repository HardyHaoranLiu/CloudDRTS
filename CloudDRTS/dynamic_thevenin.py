"""Independent dynamic Thevenin reduction and analytic frequency-response checks.

The retained C state makes the left subsystem an exact dynamic one-port.
No Circuit or Bergeron implementation is used in the reduced solver.
"""
from collections import deque
import json
from common import np,plt,table,figure,OUT
from experiments import benchmark

RS,RP,C,Z,TAU,LOAD=25.,50.,10e-6,50.,250e-6,80.
G=1/RS+1/RP

def reduced(dt,k=0,duration=.06,frequency=400):
    if dt<=0 or not isinstance(k,(int,np.integer)) or k<0:
        raise ValueError('Use positive dt and an integer nonnegative delay.')
    physical=round(TAU/dt)
    if physical<1 or not np.isclose(physical*dt,TAU,atol=1e-15,rtol=1e-12):
        raise ValueError('The physical travel time must contain an integer number of steps.')
    if any(not np.isclose(x/dt,round(x/dt),atol=1e-9,rtol=0) for x in [.005,.009]):
        raise ValueError('Switch events must lie on the time grid.')
    lr=deque([0.]*(physical+k));rl=deque([0.]*(physical+k))
    t=np.arange(round(duration/dt)+1)*dt
    gc=2*C/dt;req=1/(G+gc);history=0.
    wave=np.empty((len(t),3));audit=np.empty((len(t),3))
    for j,ti in enumerate(t):
        hl,hr=rl[0],lr[0]
        source=100*np.sin(2*np.pi*frequency*ti)
        # v = E_history - R_companion * i_port; keep the loaded capacitor state.
        eth=(source/RS-history)*req
        vl=(eth/req-hl)/(1/req+1/Z)
        closed=round(.005/dt)<=j<round(.009/dt)
        vr=-hr/(1/Z+1/LOAD+(1/18 if closed else 0.))
        ic=gc*vl+history
        port=vl/Z+hl
        audit[j]=ic,port,(source-vl)/RS-vl/RP-ic-port
        history=-(gc*vl+ic)
        lr.popleft();rl.popleft()
        lr.append(-2*vl/Z-hl);rl.append(-2*vr/Z-hr)
        wave[j]=vl,vr,vr/Z+hr
    return t,wave,audit

def phasor(frequency,delay,dt=None):
    """Independent continuous/discrete harmonic solution of the two-port circuit."""
    w=2*np.pi*frequency
    wd=w if dt is None else 2*np.tan(w*dt/2)/dt
    phi=w*(TAU+delay)
    self_y=-1j/(Z*np.tan(phi));mutual_y=1j/(Z*np.sin(phi))
    matrix=np.array([[G+1j*wd*C+self_y,mutual_y],[mutual_y,1/LOAD+self_y]])
    v=np.linalg.solve(matrix,np.array([-100j/RS,0]))
    return np.r_[v,-v[1]/LOAD]

def relative(y,ref):return 100*np.sqrt(np.sum((y-ref)**2,axis=0)/np.sum(ref**2,axis=0))

def main():
    dt=5e-6;rows=[];equiv=[];max_kcl=0.;traces={}
    freqs=[60,400,1200];delays=[0,50,100,200,400]
    fig,axes=plt.subplots(1,3,figsize=(7.15,2.65))
    for ax,f in zip(axes,freqs):
        duration=max(.06,25/f)
        t,ref,_=benchmark(dt,duration=duration,frequency=f)
        masks={'all':np.ones(len(t),bool),'steady_last_10_cycles':t>=t[-1]-10/f}
        _,static_zero,_=benchmark(dt,mode='thevenin',duration=duration,frequency=f)
        _,dynamic_zero,_=reduced(dt,duration=duration,frequency=f)
        curves={name:[] for name in ['history_exchange','static_thevenin','dynamic_thevenin']}
        for us in delays:
            k=round(us/(dt*1e6))
            _,full_delayed,_=benchmark(dt,k,'split',duration=duration,frequency=f)
            _,static,_=benchmark(dt,k,'thevenin',duration=duration,frequency=f)
            _,dynamic,audit=reduced(dt,k,duration,f)
            difference=float(np.max(abs(dynamic[:,:2]-full_delayed[:,:2])))
            current_difference=float(np.max(abs(dynamic[:,2]-full_delayed[:,2])))
            max_kcl=max(max_kcl,float(np.max(abs(audit[:,2]))))
            equiv.append([f,us,difference,current_difference])
            assert difference<1e-8 and current_difference<1e-8
            for name,y,zero in [('history_exchange',full_delayed,ref),
                                ('static_thevenin',static,static_zero),
                                ('dynamic_thevenin',dynamic,dynamic_zero)]:
                for window,mask in masks.items():
                    total=relative(y[mask,1:],ref[mask,1:])
                    model=relative(zero[mask,1:],ref[mask,1:])
                    delay=relative(y[mask,1:],zero[mask,1:])
                    rows.append([f,dt*1e6,us,name,window,*total,*model,*delay])
                    if window=='steady_last_10_cycles':curves[name].append(total[0])
            if f==400 and us in [0,100]:
                traces[f't_s_{us}']=t
                for name,y in [('full_delayed',full_delayed),('static',static),('dynamic',dynamic)]:
                    traces[f'vr_{name}_{us}']=y[:,1]
                    traces[f'ir_{name}_{us}']=y[:,2]
        for name,ls,marker in [('history_exchange','-','o'),('static_thevenin','--','s'),('dynamic_thevenin',':','x')]:
            ax.plot(delays,curves[name],ls,marker=marker,label=name.replace('_',' '))
        ax.set(xlabel='One-way added delay (us)',title=f'{f} Hz')
    axes[0].set_ylabel('Receiving voltage RMS error (%)')
    axes[0].legend(fontsize=5.8)
    figure(fig,'dynamic_thevenin_delay')
    table('dynamic_thevenin_metrics',['frequency_Hz','dt_us','added_delay_us','variant','window',
          'total_v_pct','total_i_pct','zero_delay_model_v_pct','zero_delay_model_i_pct',
          'delay_only_v_pct','delay_only_i_pct'],rows)
    table('dynamic_thevenin_equivalence',['frequency_Hz','added_delay_us','max_voltage_difference_V','max_current_difference_A'],equiv)
    table('dynamic_thevenin_traces',traces.keys(),zip(*traces.values()))

    # A separate harmonic solution checks the scalar transient solver, and
    # compares the fixed physical model across five time steps.
    refinement=[]
    for us in [50,25,10,5,2.5]:
        dt=us*1e-6;f=400.;delay=100e-6
        t,y,audit=reduced(dt,round(delay/dt),.06,f)
        mask=t>=t[-1]-10/f
        basis=np.column_stack([np.cos(2*np.pi*f*t[mask]),np.sin(2*np.pi*f*t[mask])])
        coeff=np.linalg.lstsq(basis,y[mask],rcond=None)[0]
        fitted=coeff[0]-1j*coeff[1]
        discrete=phasor(f,delay,dt);continuous=phasor(f,delay)
        steady_error=float(max(abs(fitted[:2]-discrete[:2])))
        steady_current_error=float(abs(fitted[2]-discrete[2]))
        physical_error=float(max(abs(discrete[:2]-continuous[:2])))
        physical_current_error=float(abs(discrete[2]-continuous[2]))
        assert steady_error<1e-7 and steady_current_error<1e-7
        refinement.append([us,steady_error,physical_error,steady_current_error,physical_current_error])
    order=float(np.polyfit(np.log(np.array(refinement)[:,0]),np.log(np.array(refinement)[:,2]),1)[0])
    assert 1.98<order<2.02 and max_kcl<1e-9
    table('dynamic_thevenin_refinement',['dt_us','fit_vs_discrete_max_voltage_difference_V','discrete_vs_continuous_max_voltage_difference_V','fit_vs_discrete_current_difference_A','discrete_vs_continuous_current_difference_A'],refinement)
    fig,ax=plt.subplots(1,2,figsize=(7.15,2.6))
    ax[0].loglog(np.array(refinement)[:,0],np.array(refinement)[:,2],'o-',label=f'Order {order:.3f}')
    ax[0].set(xlabel='Time step (us)',ylabel='Maximum voltage phasor difference (V)');ax[0].legend()
    static_rows=[r for r in rows if r[0]==400 and r[2]==0 and r[3]=='static_thevenin' and r[4]=='steady_last_10_cycles']
    dynamic_rows=[r for r in rows if r[0]==400 and r[2]==0 and r[3]=='dynamic_thevenin' and r[4]=='steady_last_10_cycles']
    # Error quantities retain their separate units in CSV; this panel is voltage only.
    ax[1].bar(['Static','Dynamic'],[static_rows[0][5],dynamic_rows[0][5]])
    ax[1].set_ylabel('Zero-delay RMS voltage error (%)');ax[1].set_title('400 Hz, same left subsystem')
    figure(fig,'dynamic_thevenin_validation')
    checks=dict(exact_discrete_equivalent_max_voltage_difference_V=max(r[2] for r in equiv),
        exact_discrete_equivalent_max_current_difference_A=max(r[3] for r in equiv),
        max_kcl_residual_A=max_kcl,continuous_phasor_observed_order=order,
        max_steady_fit_vs_independent_discrete_voltage_V=max(r[1] for r in refinement),
        max_steady_fit_vs_independent_discrete_current_A=max(r[3] for r in refinement),
        static_zero_delay_400Hz_voltage_error_pct=static_rows[0][5],
        dynamic_zero_delay_400Hz_voltage_error_pct=dynamic_rows[0][5],
        static_equivalent='DC-only approximation, C omitted',
        dynamic_equivalent='Zth(s)=1/(1/Rs+1/Rp+s*C); loaded capacitor state retained',
        delay_convention='One-way k*dt added to both directional wave histories, in addition to physical TAU',
        voltage_tolerance_V=1e-8,current_tolerance_A=1e-8,passed=True)
    (OUT/'dynamic_thevenin_checks.json').write_text(json.dumps(checks,indent=2))
    print(json.dumps(checks,indent=2),flush=True)

if __name__=='__main__':main()
