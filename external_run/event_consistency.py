"""Fixed-state event initialization using Cartesian current balance.

This is an explicit adapter around the unmodified ANDES 1.9.3 integrator.
It solves Cartesian current balance at a deep three-phase fault and then
checks the unmodified ANDES algebraic equations. Only the unsaturated,
fixed-field/fixed-torque GENROU case with active-I / reactive-Z loads is supported.
"""
import numpy as np
from scipy.optimize import root

def install(system):
    if not np.allclose(system.GENROU.xd2.v,system.GENROU.xq2.v):
        raise ValueError('This adapter requires equal subtransient reactances')
    if np.any(system.GENROU.S10.v):
        raise ValueError('This adapter supports the unsaturated case only')
    if not (system.PQ.config.p2i == 1 and system.PQ.config.q2z == 1):
        raise ValueError('This adapter supports active-I/reactive-Z loads only')
    if system.Bus.n != 11 or not np.all(system.Line.u.v):
        raise ValueError('This adapter is limited to the supplied interconnected 11-bus case')
    records=[]
    nb=system.Bus.n
    line=system.Line;machine=system.GENROU
    base=np.zeros((nb,nb),complex)
    for i in range(line.n):
        a=system.Bus.idx2uid(line.bus1.v[i]);b=system.Bus.idx2uid(line.bus2.v[i])
        y=line.yhk.v[i];m=line.tap.v[i]*np.exp(1j*line.phi.v[i])
        base[a,a]+=(y+line.yh.v[i])/abs(m)**2
        base[b,b]+=y+line.yk.v[i]
        base[a,b]-=y/np.conj(m);base[b,a]-=y/m
    for i in range(system.Shunt.n):
        b=system.Bus.idx2uid(system.Shunt.bus.v[i])
        base[b,b]+=system.Shunt.u.v[i]*(system.Shunt.g.v[i]+1j*system.Shunt.b.v[i])
    buses=np.array([system.Bus.idx2uid(x) for x in machine.bus.v])
    load_buses=np.array([system.Bus.idx2uid(x) for x in system.PQ.bus.v])
    for timer in (system.Fault.tf,system.Fault.tc):
        original=timer.callback
        def callback(mask,original=original):
            action=original(mask)
            if not np.any(mask):return action
            dae=system.dae
            states=dae.x.copy()
            matrix=base.copy()
            for i in range(system.Fault.n):
                b=system.Bus.idx2uid(system.Fault.bus.v[i])
                matrix[b,b]+=system.Fault.uf.v[i]/complex(system.Fault.rf.v[i],system.Fault.xf.v[i])
            delta=machine.delta.v
            pd=machine.gd1.v*machine.e1q.v+machine.gd2.v*(machine.xd1.v-machine.xl.v)*machine.e2d.v
            pq=machine.gq1.v*machine.e1d.v+(1-machine.gq1.v)*machine.e2q.v
            impedance=machine.ra.v+1j*machine.xd2.v
            internal=np.exp(1j*delta)*(pd-1j*pq)
            source=np.zeros(nb,complex)
            np.add.at(matrix,(buses,buses),1/impedance)
            np.add.at(source,buses,internal/impedance)
            old=system.Bus.v.v*np.exp(1j*system.Bus.a.v)
            def residual(z):
                v=z[:nb]+1j*z[nb:]
                magnitude=abs(v[load_buses])
                power=(system.PQ.Ipeq.v*magnitude+1j*system.PQ.Xeq.v*magnitude**2)
                balance=matrix@v-source
                np.add.at(balance,load_buses,np.conj(power/v[load_buses]))
                return np.r_[balance.real,balance.imag]
            answer=root(residual,np.r_[old.real,old.imag],tol=1e-11)
            mismatch=float(max(abs(residual(answer.x))))
            if mismatch>1e-7:raise RuntimeError(f'Current-balance event initialization failed: {mismatch}')
            v=answer.x[:nb]+1j*answer.x[nb:]
            current=(internal-v[buses])/impedance
            rotated=current*np.exp(-1j*delta)
            vdq=v[buses]*np.exp(-1j*delta)
            iq=rotated.real;id_=-rotated.imag
            assignments={
                system.Bus.v:abs(v),system.Bus.a:np.angle(v),
                machine.Iq:iq,machine.Id:id_,machine.vq:vdq.real,machine.vd:-vdq.imag,
                machine.psi2d:pd,machine.psi2q:pq,machine.psi2:np.hypot(pd,pq),
                machine.psid:vdq.real+machine.ra.v*iq,
                machine.psiq:vdq.imag-machine.ra.v*id_,
                machine.Pe:(v[buses]*current.conj()).real,
                machine.Qe:(v[buses]*current.conj()).imag,
            }
            for variable,values in assignments.items():dae.y[variable.a]=values
            system.vars_to_models()
            system.TDS.fg_update(system.exist.pflow_tds)
            # Algebraic field-current/torque helpers have explicit definitions.
            dae.y[machine.XadIfd.a]=(machine.e1q.v+(machine.xd.v-machine.xd1.v)*
                (machine.gd1.v*id_+machine.gd2.v*(machine.e1q.v-machine.e2d.v)))
            dae.y[machine.XaqI1q.a]=(machine.e1d.v+(machine.xq.v-machine.xq1.v)*
                (-machine.gq1.v*iq+machine.gq2.v*(machine.e1d.v-machine.e2q.v)))
            dae.y[machine.te.a]=(v[buses]*current.conj()).real+machine.ra.v*abs(current)**2
            system.vars_to_models();system.TDS.fg_update(system.exist.pflow_tds)
            g=float(max(abs(dae.g)))
            if g>1e-7:raise RuntimeError(f'ANDES event residual remains {g}')
            assert np.array_equal(states,dae.x)
            records.append(dict(t=float(dae.t),current_balance_residual=mismatch,
                andes_algebraic_residual=g,max_state_change=0.,
                minimum_bus_voltage=float(min(abs(v))),
                right_limit=dict(delta_rad=machine.delta.v.tolist(),
                    omega_pu=machine.omega.v.tolist(),
                    V8_pu=float(abs(v[7])),
                    P7to8_MW=float((v[6]*np.conj((v[6]-v[7])*2/(110*(.0001+.001j)))).real*100))))
            print('CONSISTENT EVENT',records[-1],flush=True)
            return action
        timer.callback=callback
    return records
