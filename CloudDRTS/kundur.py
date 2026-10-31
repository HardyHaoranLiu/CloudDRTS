"""Four two-axis transient machines on the 11-bus Kundur network."""
from common import np, plt, table, figure, OUT
from scipy.optimize import root
from scipy.integrate import solve_ivp
from scipy.linalg import eig
import json

# Bus numbers below are one-based. Network base 100 MVA, machines 900 MVA.
H = np.array([6.5, 6.5, 6.175, 6.175])
XD, XQ, XDP, XQP = 1.8, 1.7, .30, .55
W = 2 * np.pi * 60
Y = np.zeros((11, 11), complex)

def branch(a, b, y, charging=0j):
    a, b = a-1, b-1
    Y[a,a] += y + charging/2
    Y[b,b] += y + charging/2
    Y[a,b] -= y
    Y[b,a] -= y

for a,b in [(1,5),(2,6),(3,11),(4,10)]:
    branch(a,b,1/(1j*.15/9))
for a,b,km,count in [(5,6,25,1),(6,7,10,1),(7,8,110,2),(8,9,110,2),(9,10,10,1),(10,11,25,1)]:
    branch(a,b,count/(km*(.0001+.001j)),1j*count*km*.00175)
Y[6,6] += 2j
Y[8,8] += 3.5j
load = np.zeros(11, complex)
load[[6,8]] = [9.67+1j,17.67+1j]
P = -load.real.copy()
P[1:4] = [7,7.19,7]

def volts(x):
    return np.r_[1.03,1.01,1.03,1.01,x[10:]] * np.exp(1j*np.r_[0,x[:10]])

def mismatch(x):
    V = volts(x)
    S = V*np.conj(Y@V)
    return np.r_[S.real[1:]-P[1:],(S.imag+load.imag)[4:]]

sol = root(mismatch, np.r_[np.zeros(10),np.ones(7)], tol=1e-11)
assert sol.success and np.max(abs(sol.fun)) < 1e-8
V0 = volts(sol.x)
I0 = (Y@V0)[:4]/9
delta = np.angle(V0[:4]+1j*XQ*I0)
rot = np.exp(-1j*delta)
iq, id_ = (I0*rot).real, -(I0*rot).imag
eq, ed = (V0[:4]*rot).real+XDP*id_, -(V0[:4]*rot).imag-XQP*iq
efd = eq+(XD-XDP)*id_
pm = (V0[:4]*I0.conj()).real
x0 = np.r_[delta,np.zeros(4),eq,ed]
Yload = Y+np.diag(load.conj()/abs(V0)**2)
def port(fault=False):
    y = Yload.copy()
    if fault:
        y[7,7] += 1/.0001
    return 9*np.linalg.solve(y,np.eye(11)[:,:4])
normal, short = port(), port(True)

def network(x, z):
    d, speed, eq, ed = x.reshape(4,4)
    r = np.exp(1j*d)
    a = r.conj()[:,None]*z[:4]*r
    # Vq-jVd = a (Iq-jId); Vq=Eq'-X'd Id, Vd=Ed'+X'q Iq.
    mat = np.block([[a.real,a.imag+np.eye(4)*XDP],
                    [-a.imag-np.eye(4)*XQP,a.real]])
    iq,id_ = np.linalg.solve(mat,np.r_[eq,ed]).reshape(2,4)
    i = r*(iq-1j*id_)
    v = z@i
    return iq,id_,v,(v[:4]*i.conj()).real

def rhs(x,z):
    d,w,eq,ed = x.reshape(4,4)
    iq,id_,v,pe = network(x,z)
    return np.r_[W*w,(pm-pe)/(2*H),(efd-eq-(XD-XDP)*id_)/8,
                 (-ed+(XQ-XQP)*iq)/.4]

def main():
    assert max(abs(rhs(x0,normal))) < 1e-8
    def jac(h):
        return np.column_stack([(rhs(x0+e,normal)-rhs(x0-e,normal))/(2*h) for e in np.eye(16)*h])
    vals,vec = eig(jac(1e-5))
    other = eig(jac(5e-6),right=False)
    modes = sorted(np.where((vals.imag>2) & (vals.imag<12))[0],key=lambda j:vals[j].imag)
    assert len(modes)==3
    rows=[]
    for j in modes:
        s=vec[4:8,j].copy()
        s/=s[np.argmax(abs(s))]
        a,b = abs(s[0]-s[1]),abs(s[2]-s[3])
        label = 'inter-area' if abs(s[:2].mean()-s[2:].mean())>max(a,b) else ('area-1' if a>b else 'area-2')
        rows.append([label,vals[j].real,vals[j].imag,vals[j].imag/(2*np.pi),-100*vals[j].real/abs(vals[j]),
                     *abs(s),*np.angle(s,deg=True),min(abs(other-vals[j]))])
    table('modes',['mode','sigma_per_s','omega_rad_s','freq_Hz','damping_pct',
                   *[f'G{i}_amp' for i in range(1,5)],*[f'G{i}_phase_deg' for i in range(1,5)],'jacobian_delta'],rows)
    # Restart integration at both topology transitions. Output is right-continuous.
    tlog=[]; xlog=[]; vlog=[]
    x=x0
    for a,b,z in [(0.,1.,normal),(1.,1.08,short),(1.08,16.,normal)]:
        t=np.linspace(a,b,round((b-a)/.005)+1)
        ans=solve_ivp(lambda t,x:rhs(x,z),(a,b),x,t_eval=t,method='DOP853',rtol=1e-9,atol=1e-11,max_step=.03)
        assert ans.success
        x=ans.y[:,-1]
        count=len(t) if b==16. else len(t)-1
        for ti,xi in zip(t[:count],ans.y.T[:count]):
            tlog.append(ti); xlog.append(xi); vlog.append(network(xi,z)[2])
    t=np.array(tlog); x=np.array(xlog); v=np.array(vlog)
    angle=np.rad2deg(x[:,:4]-x0[:4])
    angle-=np.average(angle,axis=1,weights=H)[:,None]
    tie=(v[:,6]*np.conj((v[:,6]-v[:,7])*2/(110*(.0001+.001j)))).real*100
    table('kundur_fault',['t_s',*[f'G{i}_rel_deg' for i in range(1,5)],*[f'G{i}_df_Hz' for i in range(1,5)],'V8_pu','P7to8_MW'],np.column_stack([t,angle,x[:,4:8]*60,abs(v[:,7]),tie]))
    table('powerflow',['bus','v_pu','angle_deg','net_P_MW','net_Q_Mvar'],zip(range(1,12),abs(V0),np.angle(V0,deg=True),100*(V0*np.conj(Y@V0)).real,100*(V0*np.conj(Y@V0)).imag))
    stats=dict(tie_at_bus7_MW=tie[0],fault_v8_min_pu=min(abs(v[:,7])),pf_residual=max(abs(sol.fun)),initial_derivative=max(abs(rhs(x0,normal))))
    (OUT/'kundur_checks.json').write_text(json.dumps(stats,indent=2))
    fig,axs=plt.subplots(3,1,figsize=(7.5,5.7),sharex=True)
    for i in range(4): axs[0].plot(t,angle[:,i],label=f'G{i+1}')
    axs[0].legend(ncol=4); axs[0].set_ylabel('Rotor angle\nchange (deg, COI)')
    axs[1].plot(t,tie,color='#965527'); axs[1].set_ylabel('7 to 8 power\n(MW)')
    axs[2].plot(t,abs(v[:,7]),color='#28645e'); axs[2].set_ylabel('Bus 8 voltage\n(p.u.)')
    for ax in axs: ax.axvspan(1,1.08,color='red',alpha=.16)
    axs[-1].set_xlabel('Time (s)'); figure(fig,'kundur_fault')
    print(rows,stats)

if __name__=='__main__': main()
