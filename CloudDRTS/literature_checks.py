"""Explicit model and reference comparison, with unmatched parameters recorded."""
import json
from scipy.optimize import root
from scipy.linalg import eig
from common import np, plt, table, figure, OUT
import kundur as kd

class Case:
    def __init__(self,Y,load,V,H,active_exponent):
        self.Y=Y; self.load=load;self.V=V;self.H=H;self.exponent=active_exponent
        I=((Y@V)+np.conj(load/V))[:4]/9
        d=np.angle(V[:4]+1j*kd.XQ*I)
        ir=I*np.exp(-1j*d);iq=ir.real;id_=-ir.imag
        vr=V[:4]*np.exp(-1j*d)
        eq=vr.real+kd.XDP*id_;ed=-vr.imag-kd.XQP*iq
        self.x=np.r_[d,np.zeros(4),eq,ed]
        self.efd=eq+(kd.XD-kd.XDP)*id_
        self.pm=(V[:4]*I.conj()).real
        self.guess=np.r_[V.real,V.imag,iq,id_]

    def network(self,x,reactance_d=kd.XDP,reactance_q=kd.XQP,resistance=0.,return_voltage=False):
        d,w,eq,ed=x.reshape(4,4)
        rotation=np.exp(1j*d)
        def residual(z):
            v=z[:11]+1j*z[11:22];iq=z[22:26];id_=z[26:30]
            scale=abs(v)/abs(self.V)
            power=self.load.real*scale**self.exponent+1j*self.load.imag*scale**2
            current=rotation*(iq-1j*id_)
            kcl=self.Y@v+np.conj(power/v)
            kcl[:4]-=9*current
            vr=v[:4]/rotation
            return np.r_[kcl.real,kcl.imag,vr.real+resistance*iq+reactance_d*id_-eq,
                         -vr.imag+resistance*id_-reactance_q*iq-ed]
        answer=root(residual,self.guess,tol=1e-10)
        if np.max(abs(residual(answer.x)))>1e-8:
            raise RuntimeError("Network equations did not converge.")
        z=answer.x;v=z[:11]+1j*z[11:22];iq=z[22:26];id_=z[26:30]
        pe=(v[:4]*np.conj(rotation*(iq-1j*id_))).real
        return (iq,id_,pe,v) if return_voltage else (iq,id_,pe)

    def rhs(self,x):
        d,w,eq,ed=x.reshape(4,4);iq,id_,pe=self.network(x)
        return np.r_[kd.W*w,(self.pm-pe)/(2*self.H),
                     (self.efd-eq-(kd.XD-kd.XDP)*id_)/8,
                     (-ed+(kd.XQ-kd.XQP)*iq)/.4]

    def modes(self):
        h=1e-5
        A=np.column_stack([(self.rhs(self.x+e)-self.rhs(self.x-e))/(2*h)
                           for e in np.eye(16)*h])
        ev,vec=eig(A)
        ix=sorted(np.where((ev.imag>1)&(ev.imag<12))[0],key=lambda i:ev[i].imag)
        return ev[ix],vec[:,ix]

def symmetric(compensation_pu=2.0,tie_circuits=1,transformer_x_pu=.15):
    """Symmetric, zero-transfer case of [2] (Table 2: 1400 MW / 1367 MW per area).

    compensation_pu: shunt capacitor at each load bus (bus 7 and bus 9), on
    the 100 MVA base; 2.0 p.u. = 200 Mvar is the tentative v7 value ([2] only
    says capacitors were added as needed).
    tie_circuits: parallel 220 km circuits between the areas (1 to 3). [2]
    Fig. 1 shows 3 x 0.22; its Table 2 manual-exciter case has 1 in service.
    transformer_x_pu: generator step-up reactance on the 900 MVA machine base.
    0.15 is the textbook network [1]; [2] Fig. 1 draws no step-up transformers
    and gives no transformer data, approximated here by a small value.
    Defaults reproduce the v7 case exactly.
    """
    if tie_circuits not in (1,2,3):
        raise ValueError('tie_circuits must be 1, 2 or 3.')
    if transformer_x_pu<=0:
        raise ValueError('transformer_x_pu must be positive; use 1e-4 to approximate no transformer.')
    Y=kd.Y.copy()
    Y[6,6]+=1j*(compensation_pu-2.0)
    Y[8,8]+=1j*(compensation_pu-3.5)
    extra=tie_circuits-2  # kd.Y holds two circuits on each 110 km section 7-8 and 8-9
    g=1/(110*(.0001+.001j));bc=1j*110*.00175/2
    for a,b in [(6,7),(7,8)]:
        Y[a,a]+=extra*(g+bc);Y[b,b]+=extra*(g+bc);Y[a,b]-=extra*g;Y[b,a]-=extra*g
    change=9/(1j*transformer_x_pu)-9/(1j*.15)
    for a,b in [(0,4),(1,5),(2,10),(3,9)]:
        Y[a,a]+=change;Y[b,b]+=change;Y[a,b]-=change;Y[b,a]-=change
    load=np.zeros(11,complex);load[[6,8]]=13.67+1j
    angle_ix=np.array([1,3,4,5,6,7,8,9,10])
    def volts(z):
        mag=np.r_[1.03,1.01,1.03,1.01,z[9:]]
        angle=np.zeros(11);angle[angle_ix]=z[:9]
        return mag*np.exp(1j*angle)
    def mismatch(z):
        v=volts(z);S=v*np.conj(Y@v)+load
        p=S.real[angle_ix].copy();p[:2]-=7
        return np.r_[p,S.imag[4:]]
    sol=root(mismatch,np.r_[np.zeros(9),np.ones(7)],tol=1e-11)
    if max(abs(sol.fun))>1e-9:
        raise RuntimeError('Symmetric power flow did not converge.')
    v=volts(sol.x)
    return Case(Y,load,v,np.full(4,6.5),2)

def main():
    # Independent all-bus nonlinear formulation verifies the original reduced network.
    z=Case(kd.Y,kd.load,kd.V0,kd.H,2)
    rng=np.random.default_rng(43)
    errors=[]
    for _ in range(12):
        x=z.x+rng.normal(0,1e-4,16)
        errors.append(max(abs(z.rhs(x)-kd.rhs(x,kd.normal))))
    assert max(errors)<1e-8
    rows=[]
    refs=[.545,1.087,1.117] # Kundur textbook manual-excitation example.
    for name,exponent in [('constant_Z',2),('constant_I_P_and_Z_Q',1)]:
        c=Case(kd.Y,kd.load,kd.V0,kd.H,exponent)
        ev,vec=c.modes()
        assert len(ev)==3
        for i,e in enumerate(ev):
            f=e.imag/(2*np.pi); s=vec[4:8,i];s=s/s[np.argmax(abs(s))]
            rows.append([name,i+1,e.real,e.imag,f,-e.real/abs(e)*100,refs[i],
                         (f/refs[i]-1)*100,*np.angle(s,deg=True)])
    table('literature_model_comparison',
          ['load_model','mode_index','sigma_per_s','omega_rad_s','frequency_Hz',
           'damping_pct','textbook_reference_Hz','frequency_difference_pct',
           'G1_phase_deg','G2_phase_deg','G3_phase_deg','G4_phase_deg'],rows)
    c=symmetric();ev,vec=c.modes()
    s=vec[4:8,0];s/=s[np.argmax(abs(s))]
    phase=abs(np.angle(np.mean(s[:2])/np.mean(s[2:]),deg=True))
    assert abs(phase-180)<1e-5
    sg=c.V[:4]*np.conj((c.Y@c.V)[:4])
    summary=dict(independent_network_rhs_max_difference=max(errors),
                 symmetric_one_tie_frequency_Hz=ev[0].imag/(2*np.pi),
                 symmetric_one_tie_damping_pct=-ev[0].real/abs(ev[0])*100,
                 symmetric_inter_area_phase_deg=phase,
                 symmetric_generation_MW=sg.real*100,
                 original_1991_manual_no_transfer_frequency_Hz=.481,
                 original_1991_manual_no_transfer_damping_pct=2.3,
                 unmatched_parameters=['fourth-order machine instead of detailed subtransient machine',
                                       'no stator resistance or saturation',
                                       'assumed 200 Mvar compensation per area, not fully specified in paper',
                                       'outer-unit output determined by power flow'])
    def native(o):
        return o.tolist() if hasattr(o,'tolist') else float(o)
    (OUT/'literature_checks.json').write_text(json.dumps(summary,default=native,indent=2))
    fig,ax=plt.subplots(1,2,figsize=(7.15,2.7))
    for name,marker in [('constant_Z','o'),('constant_I_P_and_Z_Q','s')]:
        rr=[r for r in rows if r[0]==name]
        ax[0].plot([1,2,3],[r[4] for r in rr],marker=marker,label=name.replace('_',' '))
        ax[1].plot([1,2,3],[r[5] for r in rr],marker=marker)
    ax[0].plot([1,2,3],refs,'kx',label='Textbook frequencies')
    ax[0].set(ylabel='Frequency (Hz)',xlabel='Mode index');ax[0].legend(fontsize=6)
    ax[1].set(ylabel='Damping ratio (%)',xlabel='Mode index')
    for a in ax:a.set_xticks([1,2,3])
    figure(fig,'literature_comparison')
    print(json.dumps(summary,default=native,indent=2))
if __name__=='__main__':main()
