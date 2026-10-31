"""Balanced instantaneous-waveform study on an 11-bus Kundur-derived network.
Per-phase RMS p.u. phasors; instantaneous p.u. uses the RMS voltage base.
Generators are fixed sinusoidal sources for this short network transient.
"""
import numpy as np
from functools import lru_cache
from scipy.optimize import root
from common import plt, figure, table
from kundur import V0, I0, W

DT = 10e-6
ZC = np.sqrt(.001/.00175)/2  # two identical parallel circuits
NPHYS = round(110*np.sqrt(.001*.00175)/W/DT)
PORTS = np.array([6,7,7,8])
OTHER = [1,0,3,2]
FULL = [np.arange(11)]
CUTS = {'7-8': [[0,1,4,5,6],[2,3,7,8,9,10]],
        '8-9': [[0,1,4,5,6,7],[2,3,8,9,10]]}
EXTRA = [0,10,50,100,200,400]
FAULT_ON, FAULT_OFF = .020, .060

def legacy_model():
    # B[:,j] defines the voltage drop from the + to the - end of RL branch j.
    ends=[(i,None) for i in range(4)]+[(0,4),(1,5),(2,10),(3,9),(4,5),(5,6),(8,9),(9,10)]
    x=np.r_[[.3/9]*4,[.15/9]*4,np.array([25,10,10,25])*.001]
    r=np.r_[np.zeros(8),np.array([25,10,10,25])*.0001]
    # Reactive loads: shunt inductors sized at the RMS power-flow voltage.
    ends.extend([(6,None),(8,None)])
    x=np.r_[x,abs(V0[[6,8]])**2]; r=np.r_[r,0.,0.]
    B=np.zeros((11,len(x)))
    for j,(a,b) in enumerate(ends):
        B[a,j]=1
        if b is not None:B[b,j]=-1
    # Active load is conductance; net shunt susceptance is capacitive at both load buses.
    gload=np.zeros(11); cap=np.zeros(11)
    gload[[6,8]]=[9.67,17.67]/abs(V0[[6,8]])**2
    cap[[6,8]]=[2.,3.5]
    for (a,b),km in zip(ends[8:12],[25,10,10,25]):
        cap[[a,b]]+=km*.00175/2
    e=np.zeros(len(x),complex); e[:4]=V0[:4]+1j*.3*I0
    return B,x/W,r,cap/W,gload,e

@lru_cache(maxsize=1)
def operating_point():
    """Solve physical power flow for the represented 390-us lossless-line network.

    Fixed RLC load values and generator source phasors are derived from this
    operating point once. They remain fixed during all time-step refinements.
    """
    B,L,R,C,G,E=legacy_model()
    series=1/(R[4:12]+1j*W*L[4:12])
    Y=(B[:,4:12]*series)@B[:,4:12].T+np.diag(1j*W*C)
    phi=W*NPHYS*DT
    for j in (0,2):
        a,b=PORTS[j:j+2]
        Y[a,a]+=-1j/(ZC*np.tan(phi));Y[b,b]+=-1j/(ZC*np.tan(phi))
        Y[a,b]+=1j/(ZC*np.sin(phi));Y[b,a]+=1j/(ZC*np.sin(phi))
    loads=np.zeros(11,complex);loads[[6,8]]=[9.67+1j,17.67+1j]
    target=-loads.real.copy();target[1:4]=[7,7.19,7]
    def volts(x):
        return np.r_[1.03,1.01,1.03,1.01,x[10:]]*np.exp(1j*np.r_[0,x[:10]])
    def residual(x):
        v=volts(x);power=v*np.conj(Y@v)
        return np.r_[power.real[1:]-target[1:],(power.imag+loads.imag)[4:]]
    guess=np.r_[np.angle(V0[1:]/V0[0]),abs(V0[4:])]
    answer=root(residual,guess,tol=1e-11)
    if max(abs(answer.fun))>1e-9:
        raise RuntimeError("Derived-network power flow did not converge.")
    v=volts(answer.x)
    L=L.copy();L[-2:]=abs(v[[6,8]])**2/W
    G=G.copy();G[[6,8]]=loads[[6,8]].real/abs(v[[6,8]])**2
    E=E.copy();E[:4]=v[:4]+1j*W*L[:4]*(Y@v)[:4]
    return B,L,R,C,G,E,v,Y,float(max(abs(answer.fun)))

def model():
    return tuple(item.copy() for item in operating_point()[:6])


def simulate(cut=None,k=0,dt=DT):
    if cut is not None and cut not in CUTS:
        raise ValueError("Unknown partition cut.")
    if not isinstance(k, (int, np.integer)) or k < 0:
        raise ValueError("Added delay must be a nonnegative integer.")
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("Time step must be positive and finite.")
    if any(not np.isclose(v/dt, round(v/dt), rtol=0, atol=1e-8)
           for v in (NPHYS*DT, FAULT_ON, FAULT_OFF, .1)):
        raise ValueError("Time step must divide propagation, fault times and final time.")
    B,L,R,C,G,E=model()
    omega_d=2*np.tan(W*dt/2)/dt
    adm=1/(R+1j*omega_d*L)
    g=1/(2*L/dt+R); alpha=(2*L/dt-R)*g
    gc=2*C/dt
    delays=np.full(4,round(NPHYS*DT/dt))
    if cut=='7-8':delays[:2]+=k
    if cut=='8-9':delays[2:]+=k
    # Each configuration starts at its own exact discrete 60-Hz steady state.
    Y=(B*adm)@B.T+np.diag(G+1j*omega_d*C)
    for j in (0,2):
        a,b=PORTS[j:j+2]; phi=W*dt*delays[j]
        Y[a,a]+=-1j/(ZC*np.tan(phi));Y[b,b]+=-1j/(ZC*np.tan(phi))
        Y[a,b]+=1j/(ZC*np.sin(phi));Y[b,a]+=1j/(ZC*np.sin(phi))
    V=np.linalg.solve(Y,B@(adm*E))
    phi=W*dt*delays
    iline=1j*(V[PORTS[OTHER]]-np.cos(phi)*V[PORTS])/(ZC*np.sin(phi))
    emitted=-(V[PORTS]/ZC+iline)
    t=np.arange(round(.1/dt)+1)*dt
    real=lambda ph,when: np.sqrt(2)*np.real(ph*np.exp(1j*W*when))
    # History for the t=0 solve comes from the previous sample.
    d0=B.T@V-E
    h=real(alpha*adm*d0+g*d0,-dt)
    hc=real(-(1j*omega_d*C+gc)*V,-dt)
    outgoing=np.empty((len(t)+max(delays),4))
    pad=max(delays)
    for j in range(pad):outgoing[j]=real(emitted,(j-pad)*dt)
    A=(B*g)@B.T+np.diag(gc+G)
    np.add.at(A,(PORTS,PORTS),1/ZC)
    groups=[np.array(ii) for ii in (CUTS[cut] if cut else FULL)]
    banks=[]
    for fault in [False,True]:
        M=A.copy()
        if fault:M[7,7]+=1/.01
        banks.append([np.linalg.inv(M[np.ix_(ii,ii)]) for ii in groups])
        restored=np.zeros_like(M)
        for ii in groups:restored[np.ix_(ii,ii)]=M[np.ix_(ii,ii)]
        assert np.array_equal(M,restored)
    v=np.zeros(11); y=np.empty((len(t),4))
    expected=real(V[8],t)
    src=real(E[None,:],t[:,None])
    for n,ti in enumerate(t):
        incoming=outgoing[pad+n-delays,OTHER]
        rhs=B@(g*src[n]-h)-hc
        np.add.at(rhs,PORTS,-incoming)
        state=int(round(FAULT_ON/dt)<=n<round(FAULT_OFF/dt))
        for ii,inv in zip(groups,banks[state]):v[ii]=inv@rhs[ii]
        d=B.T@v-src[n]; h=alpha*(g*d+h)+g*d
        hc=-(2*gc*v+hc)
        outgoing[pad+n]=-2*v[PORTS]/ZC-incoming
        y[n]=v[6],v[7],v[8],v[8]/ZC+incoming[3]
    assert np.max(abs(y[t<FAULT_ON,2]-expected[t<FAULT_ON]))<1e-8
    return t,y

def main():
    t,ref=simulate()
    rows=[]; trace={'t_s':t,'v9_full':ref[:,2],'i9_full':ref[:,3]}
    masks={'all':np.ones(len(t),bool),'pre':t<FAULT_ON,
           'fault':(t>=FAULT_ON)&(t<FAULT_OFF),'post':t>=FAULT_OFF}
    fig,axes=plt.subplots(1,2,figsize=(7.15,2.85))
    for cut,marker in zip(CUTS,['o','s']):
        ve=[]
        for delay in EXTRA:
            _,y=simulate(cut,round(delay/(DT*1e6)))
            if delay==0:assert np.max(abs(y-ref))<1e-10
            if delay==100:
                trace[f'v9_{cut}_100us']=y[:,2]; trace[f'i9_{cut}_100us']=y[:,3]
            for window,m in masks.items():
                dif=y[m,2:]-ref[m,2:]
                errors=100*np.sqrt(np.sum(dif**2,axis=0)/np.sum(ref[m,2:]**2,axis=0))
                rows.append([cut,delay,window,*errors,*np.max(abs(dif),axis=0)])
                if window=='all': ve.append(errors)
        for i,ax in enumerate(axes):
            ax.plot(EXTRA,np.array(ve)[:,i],marker=marker,label='Cut '+cut)
    for ax,quantity in zip(axes,['Bus 9 voltage','9-to-8 line current']):
        ax.set(xlabel='Added one-way delay (µs)',ylabel='Relative RMS error (%)',title=quantity)
    axes[0].legend();figure(fig,'grid_cut_errors')
    table('grid_cut_metrics',['cut','delay_us','window','v9_rms_pct','i9_rms_pct','v9_max_pu','i9_max_pu'],rows)
    table('grid_cut_traces',trace.keys(),zip(*trace.values()))
    fig,axes=plt.subplots(2,1,figsize=(7.15,3.45),sharex=True)
    for ax,key,unit in zip(axes,['v9','i9'],['Bus 9 voltage (p.u.)','Current into line (p.u.)']):
        for suffix,label,ls in [('full','Unsplit','-'),('7-8_100us','Cut 7–8','--'),('8-9_100us','Cut 8–9',':')]:
            ax.plot(t*1000,trace[f'{key}_{suffix}'],ls,label=label)
        ax.axvspan(20,60,alpha=.1,color='grey');ax.set_ylabel(unit)
    axes[0].legend(ncol=3);axes[-1].set_xlabel('Time (ms)');figure(fig,'grid_cut_traces')
def refinement():
    rows=[]
    for dt in [10e-6,5e-6,2.5e-6]:
        t,ref=simulate(dt=dt)
        for cut in CUTS:
            _,y=simulate(cut,round(.0001/dt),dt)
            pct=100*np.sqrt(np.sum((y[:,2:]-ref[:,2:])**2,axis=0)/np.sum(ref[:,2:]**2,axis=0))
            rows.append([dt*1e6,cut,*pct])
    table('grid_refinement',['dt_us','cut','v_rms_pct','i_rms_pct'],rows)

if __name__=='__main__':main();refinement()

