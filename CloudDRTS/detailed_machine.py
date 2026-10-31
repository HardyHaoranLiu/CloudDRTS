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

    def modes(self,h=1e-5):
        A=np.column_stack([(self.rhs(self.x+e)-self.rhs(self.x-e))/(2*h)
                           for e in np.eye(24)*h])
        vals,vec=eig(A)
        ix=sorted(np.where((vals.imag>1)&(vals.imag<12))[0],key=lambda i:vals[i].imag)
        if len(ix)!=3:raise RuntimeError('Expected three electromechanical modes.')
        return vals[ix],vec[:,ix]

def main():
    data=[];summaries=[]
    cases={'textbook':Case(kd.Y,kd.load,kd.V0,kd.H,1),
           'paper_symmetric':symmetric()}
    for name,case in cases.items():
        for sat in ['none','kundur_exponential','genrou_quadratic']:
            model=Detailed(case,sat)
            values,vectors=model.modes()
            check=model.modes(5e-6)[0]
            difference=max(abs(values-check))
            assert difference<1e-6
            for j,val in enumerate(values):
                speed=vectors[4:8,j];speed/=speed[np.argmax(abs(speed))]
                ref=[.545,1.087,1.117][j] if name=='textbook' else (.481 if j==0 else np.nan)
                zref=[3.2,7.2,7.2][j] if name=='textbook' else (2.3 if j==0 else np.nan)
                data.append([name,sat,j+1,val.real,val.imag,val.imag/(2*np.pi),
                             -val.real/abs(val)*100,ref,zref,
                             (val.imag/(2*np.pi)/ref-1)*100,
                             *np.angle(speed,deg=True)])
            summaries.append(dict(case=name,saturation=sat,
                                  equilibrium_max_residual=max(abs(model.rhs(model.x))),
                                  jacobian_halving_max_eigenvalue_change=float(difference)))
    table('detailed_modes',['case','saturation','mode','sigma_per_s','omega_rad_s',
          'frequency_Hz','damping_pct','reference_Hz','reference_damping_pct',
          'frequency_error_pct','G1_phase_deg','G2_phase_deg','G3_phase_deg','G4_phase_deg'],data)
    (OUT/'detailed_checks.json').write_text(json.dumps(summaries,indent=2))
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
