"""Six-state round-rotor model using transient and damper flux states.

Machine topology follows standard GENROU equations documented by ANDES.
Two declared saturation laws permit comparison without fitting eigenvalues.
"""
import json
from scipy.linalg import eig
from common import np, plt, table, figure, OUT
from literature_checks import Case, symmetric
import kundur as kd

RA, XL, XD, XQ, XD1, XQ1, X2 = .0025,.2,1.8,1.7,.3,.55,.25
GD1=(X2-XL)/(XD1-XL); GQ1=(X2-XL)/(XQ1-XL)
GD2=(XD1-X2)/(XD1-XL)**2;GQ2=(XQ1-X2)/(XQ1-XL)**2
GQD=(XQ-XL)/(XD-XL)

def saturation(psi,kind):
    psi=np.asarray(psi)
    if kind=='none':return np.zeros_like(psi)
    if kind=='kundur_exponential':
        return np.where(psi>.9,.015*np.exp(9.6*(psi-.9))/psi,0.)
    if kind=='genrou_quadratic':
        s1=.015*np.exp(.96)
        s12=.015*np.exp(2.88)/1.2
        ratio=np.sqrt(1.2*s12/s1)
        a=(1.2-ratio)/(1-ratio)
        b=s1/(1-a)**2
        return b*np.maximum(psi-a,0)**2/psi
    raise ValueError('Unknown saturation model.')

class Detailed:
    def __init__(self,case,sat='none'):
        self.case=case;self.sat=sat
        V=case.V[:4]
        I=(case.Y@case.V+np.conj(case.load/case.V))[:4]/9
        psi_mag=abs(V+(RA+1j*X2)*I)
        se=saturation(psi_mag,sat)
        xq_effective=X2+(XQ-X2)/(1+se*GQD)
        delta=np.angle(V+(RA+1j*xq_effective)*I)
        vr=V*np.exp(-1j*delta);ir=I*np.exp(-1j*delta)
        iq=ir.real;id_=-ir.imag
        psi2d=vr.real+RA*iq+X2*id_
        psi2q=-vr.imag+RA*id_-X2*iq
        e1q=psi2d+(XD1-X2)*id_
        e1d=psi2q-(XQ1-X2)*iq
        e2d=e1q-(XD1-XL)*id_
        e2q=e1d+(XQ1-XL)*iq
        self.x=np.r_[delta,np.zeros(4),e1q,e1d,e2d,e2q]
        self.efd=se*psi2d+e1q+(XD-XD1)*id_
        self.tm=(V*np.conj(I)).real+RA*abs(I)**2
        assert max(abs(self.rhs(self.x)))<1e-8

    def electrical(self,x,return_voltage=False):
        d,w,e1q,e1d,e2d,e2q=x.reshape(6,4)
        psi2d=GD1*e1q+GD2*(XD1-XL)*e2d
        psi2q=GQ1*e1d+(1-GQ1)*e2q
        result=self.case.network(np.r_[d,w,psi2d,psi2q],X2,X2,RA,return_voltage)
        iq,id_,pe=result[:3]
        quantities=(psi2d,psi2q,iq,id_,pe)
        return quantities+(result[3],) if return_voltage else quantities

    def rhs(self,x):
        d,w,e1q,e1d,e2d,e2q=x.reshape(6,4)
        pd,pq,iq,id_,pe=self.electrical(x)
        se=saturation(np.hypot(pd,pq),self.sat)
        field=se*pd+e1q+(XD-XD1)*(GD1*id_+GD2*(e1q-e2d))
        quadrature=se*GQD*pq+e1d+(XQ-XQ1)*(-GQ1*iq+GQ2*(e1d-e2q))
        torque=pe+RA*(iq**2+id_**2)
        return np.r_[kd.W*w,(self.tm-torque)/(2*self.case.H),
                     (self.efd-field)/8,-quadrature/.4,
                     (-id_*(XD1-XL)+e1q-e2d)/.03,
                     (iq*(XQ1-XL)+e1d-e2q)/.05]

    def state_matrix(self,h=1e-5):
        return np.column_stack([(self.rhs(self.x+e)-self.rhs(self.x-e))/(2*h)
                                for e in np.eye(24)*h])

    def spectrum(self,h=1e-5):
        """All 24 eigenvalues, used to check slow aperiodic modes as well."""
        return eig(self.state_matrix(h),right=False)

    def modes(self,h=1e-5):
        vals,vec=eig(self.state_matrix(h))
        ix=sorted(np.where((vals.imag>1)&(vals.imag<12))[0],key=lambda i:vals[i].imag)
        if len(ix)!=3:raise RuntimeError('Expected three electromechanical modes.')
        return vals[ix],vec[:,ix]

def static_exciter_modes(model,KA=200.,TR=.01,h=1e-6):
    """Electromechanical modes with a thyristor (static) exciter on every unit.

    Efd = KA (Vref - Vm), dVm/dt = (Vt - Vm)/TR, no transient gain reduction;
    [1] Table E12.4(ii) and [2] Appendix A use KA = 200, TR = 0.01 s. Limits are
    not active in a small-signal study; Vref keeps the manual equilibrium.
    Returns the three modes sorted by frequency with their speed mode shapes.
    """
    vt0=abs(model.electrical(model.x,True)[-1][:4]);efd0=model.efd.copy();vref=vt0+efd0/KA
    x0=np.r_[model.x,vt0]
    def f(x):
        xs,vm=x[:24],x[24:]
        vt=abs(model.electrical(xs,True)[-1][:4])
        model.efd=KA*(vref-vm)
        try:
            return np.r_[model.rhs(xs),(vt-vm)/TR]
        finally:
            model.efd=efd0
    if max(abs(f(x0)))>1e-7:raise RuntimeError('Exciter equilibrium is inconsistent.')
    A=np.column_stack([(f(x0+e)-f(x0-e))/(2*h) for e in np.eye(28)*h])
    vals,left,vec=eig(A,left=True)
    participation=abs(left.conj()*vec)
    mechanical=participation[:8].sum(0)/participation.sum(0)
    ix=sorted([i for i in range(28) if .5<vals[i].imag<12 and mechanical[i]>.5],key=lambda i:vals[i].imag)
    if len(ix)!=3:raise RuntimeError(f'Expected three electromechanical modes, found {len(ix)}.')
    return vals[ix],vec[4:8][:,ix]

# Reference values keyed by the computed mode label (frequency Hz, damping %),
# all checked against the sources.
# textbook: Kundur [1] Example 12.6, Table E12.3, manual excitation; [1] identifies
# 1.087 Hz as the area-1 (G1 vs G2) mode and 1.117 Hz as the area-2 mode. The
# example includes saturation (Asat 0.015, Bsat 9.6, psiT1 0.9), i.e. the
# 'kundur_exponential' variant here. That variant applies secant (GENROU-type)
# saturation; textbook_reproduction.py uses the incremental treatment of [1] and
# reproduces Tables E12.3 and E12.4 eigenvalue by eigenvalue.
# paper_symmetric: Klein, Rogers and Kundur [2], Table 2, manual excitation, one tie
# circuit, zero transfer. [2] gives no local-mode values.
REFERENCES={'textbook':{'inter-area':(.545,3.2),'area-1':(1.087,7.2),'area-2':(1.117,7.2)},
            'paper_symmetric':{'inter-area':(.481,2.3)}}
# Positive real eigenvalues are reported above 1e-3 s^-1: the double zero eigenvalue
# (no infinite bus, KD = 0) splits to about +-2e-4 under finite differences.
# Eigenvalues of [1]: Table E12.3 (manual) and Table E12.4(ii) (thyristor, KA 200, TR 0.01).
TEXTBOOK_EIGENVALUES={'manual':{'inter-area':-.111+3.43j,'area-1':-.492+6.82j,'area-2':-.506+7.02j},
                      'thyristor':{'inter-area':.031+3.84j,'area-1':-.490+7.15j,'area-2':-.496+7.35j}}

def main():
    data=[];summaries=[];book=[]
    cases={'textbook':Case(kd.Y,kd.load,kd.V0,kd.H,1),
           'paper_symmetric':symmetric()}
    for name,case in cases.items():
        for sat in ['none','kundur_exponential','genrou_quadratic']:
            model=Detailed(case,sat)
            values,vectors=model.modes()
            check=model.modes(5e-6)[0]
            difference=max(abs(values-check))
            assert difference<1e-6
            labels=[kd.mode_label(vectors[4:8,j]) for j in range(3)]
            expected=(['area-1','area-2','inter-area'] if name=='textbook'
                      else ['inter-area','local-both-areas','local-both-areas'])
            if sorted(labels)!=expected or labels[0]!='inter-area':
                raise RuntimeError(f'Unexpected mode classification {labels} for {name}/{sat}.')
            for j,val in enumerate(values):
                speed=vectors[4:8,j].copy();speed/=speed[np.argmax(abs(speed))]
                ref,zref=REFERENCES[name].get(labels[j],(np.nan,np.nan))
                data.append([name,sat,j+1,val.real,val.imag,val.imag/(2*np.pi),
                             -val.real/abs(val)*100,ref,zref,
                             (val.imag/(2*np.pi)/ref-1)*100,
                             *np.angle(speed,deg=True),labels[j],*abs(speed)])
            spectrum=model.spectrum()
            summaries.append(dict(case=name,saturation=sat,
                                  equilibrium_max_residual=max(abs(model.rhs(model.x))),
                                  jacobian_halving_max_eigenvalue_change=float(difference),
                                  max_real_eigenvalue=float(spectrum.real.max()),
                                  positive_real_eigenvalues=[float(z.real) for z in spectrum
                                                             if z.real>1e-3],
                                  near_zero_eigenvalues=int(sum(abs(spectrum)<1e-3))))
            if name=='textbook' and sat in ('none','kundur_exponential'):
                for excitation in ('manual','thyristor'):
                    if excitation=='manual':
                        vals,shapes=values,vectors[4:8]
                    else:
                        vals,shapes=static_exciter_modes(Detailed(case,sat))
                    for j in range(3):
                        label=kd.mode_label(shapes[:,j]);ref=TEXTBOOK_EIGENVALUES[excitation][label]
                        val=vals[j]
                        book.append([sat,excitation,label,val.real,val.imag/(2*np.pi),
                                     -100*val.real/abs(val),ref.real,ref.imag/(2*np.pi),
                                     -100*ref.real/abs(ref),(val.imag/ref.imag-1)*100])
    table('detailed_modes',['case','saturation','mode','sigma_per_s','omega_rad_s',
          'frequency_Hz','damping_pct','reference_Hz','reference_damping_pct',
          'frequency_error_pct','G1_phase_deg','G2_phase_deg','G3_phase_deg','G4_phase_deg',
          'mode_label','G1_amp','G2_amp','G3_amp','G4_amp'],data)
    (OUT/'detailed_checks.json').write_text(json.dumps(summaries,indent=2))
    table('textbook_eigenvalue_comparison',['saturation','excitation','mode_label','sigma_per_s',
          'frequency_Hz','damping_pct','ref1_sigma_per_s','ref1_frequency_Hz','ref1_damping_pct',
          'frequency_difference_pct'],book)
    fig,ax=plt.subplots(1,2,figsize=(7.15,2.8))
    for sat,marker in zip(['none','kundur_exponential','genrou_quadratic'],['o','s','^']):
        rs=[r for r in data if r[0]=='textbook' and r[1]==sat]
        ax[0].plot([1,2,3],[r[5] for r in rs],marker=marker,label=sat.replace('_',' '))
        ax[1].plot([1,2,3],[r[6] for r in rs],marker=marker)
    ax[0].plot([1,2,3],[.545,1.087,1.117],'kx',label='Textbook')
    ax[1].plot([1,2,3],[3.2,7.2,7.2],'kx')
    ax[0].legend(fontsize=6)
    for a in ax:a.set(xlabel='Mode index',xticks=[1,2,3])
    ax[0].set_ylabel('Frequency (Hz)');ax[1].set_ylabel('Damping ratio (%)')
    figure(fig,'detailed_modes')
    for r in data:print(r[:10])

if __name__=='__main__':main()
