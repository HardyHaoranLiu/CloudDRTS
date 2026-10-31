"""Independent analytic checks and regression evidence for the audit."""
import json
import numpy as np
from common import OUT, table
from emt import Circuit, Bergeron
from experiments import benchmark
from grid_cuts import simulate, model, PORTS, ZC, NPHYS, DT
from kundur import V0, W

def require_error(func):
    try:
        func()
    except ValueError:
        return True
    raise AssertionError("Invalid input was accepted.")

def main():
    guards = [
        require_error(lambda: Circuit(1, 0)),
        require_error(lambda: Circuit(1, 1e-3).reactive('wrong',1,0,1)),
        require_error(lambda: Circuit(1, 1e-3).resistor(1,0,-1)),
        require_error(lambda: Bergeron(-50,.001,.001)),
        require_error(lambda: Bergeron(50,.001,.001,k=1.5)),
        require_error(lambda: benchmark(5e-6,mode='typo')),
        require_error(lambda: simulate('typo')),
        require_error(lambda: simulate(dt=7e-6)),
    ]
    # Mismatched source and load cause repeated reflections with a known lattice solution.
    # With Z=50, Rs=100, RL=150: source reflection=1/3, load reflection=1/2.
    errors=[]
    for dt in (50e-6, 5e-6):
        tau=.0005; n=round(tau/dt); line=Bergeron(50,tau,dt)
        c=Circuit(2,dt)
        c.resistor(-1,1,100); c.resistor(1,0,50)
        c.resistor(2,0,150); c.resistor(2,0,50)
        t=np.arange(round(.006/dt)+1)*dt; out=[]
        for _ in t:
            hl,hr=line.incoming()
            vl,vr=c.step({-1:1},[-hl,-hr])
            line.publish(vl,vr,hl,hr);out.append(vr)
        j=np.arange(len(t)); count=np.maximum(0,(j-n)//(2*n)+1)
        exact=.6*(1-(1/6)**count)
        err=max(abs(np.array(out)-exact)); assert err<1e-12
        errors.append([dt*1e6,err])
        table(f'audit_reflections_{round(dt*1e6)}us',['t_s','computed_V','analytic_V'],zip(t,out,exact))
    # Right-hand event limit preserves charge and restarts the correct history.
    switch_rows=[]
    for us in [50,25,10,5,2.5]:
        dt=us*1e-6; c=Circuit(1,dt)
        c.resistor(-1,1,10);c.reactive('C',1,0,100e-6,i0=.1);c.switch(1,0,10)
        t=np.arange(round(.005/dt)+1)*dt;v=[0.]
        for j,ti in enumerate(t[1:],1):
            event=round(.002/dt)
            state=j>event
            point=c.step({-1:1},[0],(state,))[0]
            if j==event:
                point=c.reinitialize({-1:1},[0],(True,))[0]
            v.append(point)
        exact=np.where(t<.002,1-np.exp(-t/.001),
                       .5+(1-np.exp(-2)-.5)*np.exp(-(t-.002)/.0005))
        err=max(abs(np.asarray(v)-exact));switch_rows.append([us,err])
        table(f'audit_switch_{us:g}us',['t_s','computed_V','analytic_V','error_V'],
              zip(t,v,exact,np.array(v)-exact))
    switch_order=float(np.polyfit(np.log(np.array(switch_rows)[:,0]),
                                  np.log(np.array(switch_rows)[:,1]),1)[0])
    assert 1.95 < switch_order < 2.05
    table('audit_switch_convergence',['dt_us','max_error_V'],switch_rows)
    dt=5e-6
    # Undamped LC should conserve its quadratic energy under trapezoidal integration.
    c=Circuit(1,dt);c.reactive('L',1,0,.01,v0=1);c.reactive('C',1,0,10e-6,v0=1)
    peak=0.
    for _ in range(4000):
        vl=c.step({},[0])[0]
        branch=c.storage[0]
        # Stored L history is i+g*v. Recover the just-computed inductor current.
        il=branch[-1]-branch[2]*vl
        energy=.5*.01*il**2+.5*10e-6*vl**2
        peak=max(peak,abs(energy/(.5*10e-6)-1))
    assert peak<1e-10
    # Quantify changed operating point in the derived EMT model.
    B,L,R,C,G,E=model();wd=2*np.tan(W*DT/2)/DT;ad=1/(R+1j*wd*L)
    Y=(B*ad)@B.T+np.diag(G+1j*wd*C)
    for j in (0,2):
        a,b=PORTS[j:j+2];phi=W*DT*NPHYS
        Y[a,a]+=-1j/(ZC*np.tan(phi));Y[b,b]+=-1j/(ZC*np.tan(phi))
        Y[a,b]+=1j/(ZC*np.sin(phi));Y[b,a]+=1j/(ZC*np.sin(phi))
    V=np.linalg.solve(Y,B@(ad*E))
    table('audit_operating_point',['bus','phasor_v_pu','derived_emt_v_pu','delta_pu'],
          zip(range(1,12),abs(V0),abs(V),abs(V)-abs(V0)))
    summary=dict(input_regressions_passed=sum(guards),multiple_reflection_checks=errors,
                 switch_convergence=switch_rows,switch_observed_order=switch_order,lc_max_relative_energy_error=peak,
                 emt_vs_phasor_max_voltage_difference_pu=float(max(abs(abs(V)-abs(V0)))))
    (OUT/'audit_checks.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
