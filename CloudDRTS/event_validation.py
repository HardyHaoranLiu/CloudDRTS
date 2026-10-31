"""Validate both capacitor voltage and inductor current across two switch events."""
import json
from scipy.linalg import expm
from common import np, plt, figure, table, OUT
from emt import Circuit

def exact_segment(x,t,rsh):
    R,L,C=6.,.01,10e-6
    A=np.array([[-R/L,-1/L],[1/C,-(0 if rsh is None else 1/rsh)/C]])
    b=np.array([1/L,0.])
    steady=-np.linalg.solve(A,b)
    return steady+expm(A*t)@(x-steady)

def main():
    rows=[]; continuity=[]
    for us in [50,25,10,5,2.5]:
        dt=us*1e-6;on,off=round(.002/dt),round(.005/dt)
        c=Circuit(2,dt)
        c.resistor(-1,1,6.);c.reactive('L',1,2,.01,v0=1.)
        c.reactive('C',2,0,10e-6);c.switch(2,0,100.)
        times=np.arange(round(.008/dt)+1)*dt
        at_on=exact_segment(np.zeros(2),.002,None)
        at_off=exact_segment(at_on,.003,100.)
        wave=[];exact=[]
        for j,t in enumerate(times):
            if j==0:
                voltage=0.;current=0.
            else:
                before=on<j<=off
                point=c.step({-1:1.},[0.,0.],(before,))
                if j in (on,off):
                    energy_before=np.array([c.storage_state[0][1],c.storage_state[1][0]])
                    point=c.reinitialize({-1:1.},[0.,0.],(j==on,))
                    energy_after=np.array([c.storage_state[0][1],c.storage_state[1][0]])
                    continuity.append(max(abs(energy_after-energy_before)))
                voltage=point[1];current=c.storage_state[0][1]
            wave.append([current,voltage])
            if j<=on:ref=exact_segment(np.zeros(2),t,None)
            elif j<=off:ref=exact_segment(at_on,t-.002,100.)
            else:ref=exact_segment(at_off,t-.005,None)
            exact.append(ref)
        wave=np.array(wave);exact=np.array(exact)
        errors=np.max(abs(wave-exact),axis=0);rows.append([us,*errors])
        table(f'event_rlc_{us:g}us',['t_s','i_A','v_V','exact_i_A','exact_v_V'],np.column_stack([times,wave,exact]))
    orders=[float(np.polyfit(np.log(np.array(rows)[:,0]),np.log(np.array(rows)[:,i]),1)[0]) for i in [1,2]]
    assert all(1.95<x<2.05 for x in orders)
    assert max(continuity)<1e-12
    table('event_rlc_convergence',['dt_us','max_current_error_A','max_voltage_error_V'],rows)
    rc=np.genfromtxt(OUT/'audit_switch_convergence.csv',delimiter=',',names=True)
    fig,axes=plt.subplots(1,2,figsize=(7.15,2.7))
    axes[0].loglog(rc['dt_us'],rc['max_error_V'],'o-',label='Switched RC')
    axes[0].loglog(np.array(rows)[:,0],np.array(rows)[:,2],'s--',label='Switched RLC')
    axes[0].legend();axes[0].set(xlabel='Time step (µs)',ylabel='Maximum voltage error (V)')
    axes[1].loglog(np.array(rows)[:,0],np.array(rows)[:,1],'o-')
    axes[1].set(xlabel='Time step (µs)',ylabel='Maximum inductor-current error (A)')
    figure(fig,'event_convergence')
    (OUT/'event_validation.json').write_text(json.dumps(dict(rlc_current_order=orders[0],
            rlc_voltage_order=orders[1],max_event_state_discontinuity=max(continuity)),indent=2))
if __name__=='__main__':main()
