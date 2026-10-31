from common import np, plt, table, figure, OUT
from emt import Circuit, Bergeron
import time, platform, sys, scipy, json, cProfile, pstats, io

def validation():
    fig,axes=plt.subplots(1,3,figsize=(9,2.9))
    rows=[]
    for ax,case in zip(axes,['RC','RL','line']):
        last=None
        for dt in [50e-6,5e-6]:
            n=round(.004/dt); t=np.arange(n+1)*dt
            if case=='line':
                z,tau=50.,.0005
                # Source is stepped at t=0; matched source, load=3Z.
                line=Bergeron(z,tau,dt)
                circuit=Circuit(2,dt)
                circuit.resistor(-1,1,z); circuit.resistor(1,0,z)
                circuit.resistor(2,0,3*z); circuit.resistor(2,0,z)
                y=[]
                for _ in t:
                    hl,hr=line.incoming()
                    vl,vr=circuit.step({-1:1.},[-hl,-hr])
                    line.publish(vl,vr,hl,hr); y.append(vr)
                exact=.75*(np.arange(n+1)>=round(tau/dt))
            else:
                circuit=Circuit(1,dt); circuit.resistor(-1,1,10.)
                if case=='RC':
                    circuit.reactive('C',1,0,100e-6,i0=.1)
                else:
                    circuit.reactive('L',1,0,.01,v0=1.)
                y=[0.]
                for _ in t[1:]:
                    v=circuit.step({-1:1.},[0.])[0]
                    y.append(v if case=='RC' else (1-v)/10)
                exact=(1-np.exp(-t/.001))*(1. if case=='RC' else .1)
            err=np.asarray(y)-exact
            maxerr=max(abs(err)); rms=np.sqrt(np.mean(err**2))
            rows.append([case,dt*1e6,maxerr,rms])
            table(f'validate_{case}_{round(dt*1e6)}us',['time_s','computed','analytic','error'],zip(t,y,exact,err))
            ax.plot(t*1000,y,label=f'{dt*1e6:g} µs',linestyle='-' if dt>1e-5 else '--')
            if last is not None and case!='line': assert maxerr < last/90
            if case=='line': assert maxerr < 1e-12
            last=maxerr
        ax.plot(t[::20]*1000,exact[::20],'.',ms=2.5,color='black',label='analytic')
        ax.set(xlabel='Time (ms)',ylabel='A' if case=='RL' else 'V',title=case)
    axes[0].legend(); figure(fig,'validation')
    table('validation',['case','dt_us','max_abs_error','rms_error'],rows)

Z,TAU=50.,250e-6
def benchmark(dt, k=0, mode='unsplit',timed=False,duration=.020,frequency=400):
    if mode not in {"unsplit", "split", "thevenin"}:
        raise ValueError("Unknown boundary model.")
    n=round(duration/dt)+1
    t=np.arange(n)*dt
    line=Bergeron(Z,TAU,dt,k)
    same=mode=='unsplit'
    left=Circuit(2 if same else 1,dt)
    right=left if same else Circuit(1,dt)
    rn=2 if same else 1
    eq=mode=='thevenin'
    if eq:
        left.resistor(-1,1,1/(1/25+1/50))
    else:
        left.resistor(-1,1,25); left.resistor(1,0,50)
        left.reactive('C',1,0,10e-6)
    left.resistor(1,0,Z); right.resistor(rn,0,Z)
    right.resistor(rn,0,80); right.switch(rn,0,18)
    wave=np.empty((n,3)); ns=np.empty(n)
    for j,ti in enumerate(t):
        if timed: start=time.perf_counter_ns()
        vs=100*np.sin(2*np.pi*frequency*ti)
        on,off=round(.005/dt),round(.009/dt)
        closed=on<=j<off
        event=j in (on,off)
        before=not closed if event else closed
        hl,hr=line.incoming()
        if same:
            vl,vr=left.step({-1:vs},[-hl,-hr],(before,))
            if event:
                vl,vr=left.reinitialize({-1:vs},[-hl,-hr],(closed,))
        else:
            vl=left.step({-1:vs*(2/3 if eq else 1.)},[-hl])[0]
            vr=right.step({},[-hr],(before,))[0]
            if event:
                vr=right.reinitialize({},[-hr],(closed,))[0]
        line.publish(vl,vr,hl,hr)
        wave[j]=vl,vr,vr/Z+hr
        if timed: ns[j]=time.perf_counter_ns()-start
    return t,wave,ns

def interfaces():
    dt=5e-6
    t,ref,_=benchmark(dt)
    table('interface_reference',['t_s','vl_V','vr_V','ir_into_line_A'],np.column_stack([t,ref]))
    rows=[]; fig,axs=plt.subplots(1,2,figsize=(8.4,3.2))
    for name in ['split','thevenin']:
        curve=[]
        for k in [0,1,2,5,10,20,40,80]:
            _,wave,_=benchmark(dt,k,name)
            diff=wave[:,1:]-ref[:,1:]
            rms=np.sqrt(np.mean(diff**2,axis=0))
            pct=100*rms/np.sqrt(np.mean(ref[:,1:]**2,axis=0))
            rows.append([name,k,k*5,*rms,*pct,*np.max(abs(diff),axis=0)])
            curve.append([k*5,*pct])
            if name=='split' and k==0: assert np.max(abs(diff))<1e-10
            if k in (0,20):
                table(f'{name}_k{k}',['t_s','vl_V','vr_V','ir_into_line_A'],np.column_stack([t,wave]))
        curve=np.array(curve)
        for col,ax in enumerate(axs,1):
            ax.plot(curve[:,0],curve[:,col],'o-',label=name)
            ax.set(xlabel='Added one-way delay (µs)',ylabel='Relative RMS error (%)',title=['Receiving voltage','Current into line'][col-1])
    axs[0].legend(); figure(fig,'delay')
    table('interface_errors',['variant','k','extra_us','rms_V','rms_A','relative_v_pct','relative_i_pct','max_V','max_A'],rows)

def profile():
    rows=[]; timings=[]
    for us in [5,50]:
        dt=us*1e-6
        benchmark(dt,duration=dt*500)
        for rep in range(5):
            _,_,ns=benchmark(dt,timed=True,duration=dt*3999)
            timings.extend((us,rep+1,j,x) for j,x in enumerate(ns))
            rows.append([us,rep+1,len(ns),*np.percentile(ns/1000,[50,95,99,100]),np.mean(ns>us*1000)*100])
    table('timing',['budget_us','rep','steps','median_us','p95_us','p99_us','max_us','miss_pct'],rows)
    table('timing_raw',['dt_us','rep','step','elapsed_ns'],timings)
    p=cProfile.Profile(); p.runcall(benchmark,5e-6)
    buf=io.StringIO(); pstats.Stats(p,stream=buf).strip_dirs().sort_stats('cumulative').print_stats(20)
    (OUT/'profile.txt').write_text(buf.getvalue())
    (OUT/'environment.json').write_text(json.dumps(dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),architecture=platform.machine()),indent=2))

if __name__=='__main__':
    validation(); interfaces(); profile()
