"""Positive-voltage, fixed-state algebraic initialization at ANDES fault events."""
import numpy as np
from scipy.optimize import least_squares

def install(system):
    records=[]
    for timer in (system.Fault.tf,system.Fault.tc):
        original=timer.callback
        def callback(mask,original=original):
            action=original(mask)
            if not np.any(mask):return action
            dae=system.dae
            states=dae.x.copy()
            seed=dae.y.copy()
            lower=np.full(dae.m,-np.inf);upper=np.full(dae.m,np.inf)
            lower[system.Bus.v.a]=1e-8
            seed[system.Bus.v.a]=np.maximum(seed[system.Bus.v.a],1e-5)
            def residual(y):
                dae.y[:]=y
                system.vars_to_models()
                system.TDS.fg_update(system.exist.pflow_tds)
                r=dae.g.copy()
                voltage=y[system.Bus.v.a]
                r[system.Bus.a.a]/=voltage
                r[system.Bus.v.a]/=voltage
                return r
            target=system.Fault.uf.v.copy()
            schedule=([.00001,.00003,.0001,.0003,.001,.003,.01,.03,.1,.3,1.]
                      if np.any(target) else [.3,.1,.03,.01,.003,.001,.0003,.0001,.00003,.00001,0.])
            for factor in schedule:
                system.Fault.uf.v[:]=factor
                answer=least_squares(residual,seed,bounds=(lower,upper),x_scale='jac',
                                     ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=100)
                seed=answer.x.copy()
                mismatch=float(max(abs(residual(seed))))
                if mismatch>1e-7:
                    raise RuntimeError(f'Event continuation failed at {factor}: {mismatch}')
            system.Fault.uf.v[:]=target
            mismatch=float(max(abs(residual(seed))))
            if mismatch>1e-7:raise RuntimeError(f'Event algebraic initialization failed: {mismatch}')
            assert np.array_equal(states,dae.x)
            records.append(dict(t=float(dae.t),max_current_scaled_residual=mismatch,
                                max_state_change=float(max(abs(dae.x-states))),
                                minimum_bus_voltage=float(min(system.Bus.v.v))))
            print('CONSISTENT EVENT',records[-1],flush=True)
            return action
        timer.callback=callback
    return records
