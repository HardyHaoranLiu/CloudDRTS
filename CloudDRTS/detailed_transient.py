"""Time-domain checks of the detailed machine against its modal linearization."""
import json
from scipy.integrate import solve_ivp
from common import np, plt, table, figure, OUT
from detailed_machine import Detailed
from literature_checks import Case
import kundur as kd

def main():
    case=Case(kd.Y,kd.load,kd.V0,kd.H,1)
    model=Detailed(case,'none')
    vals,vecs=model.modes()
    fig,axs=plt.subplots(1,3,figsize=(7.15,2.4))
    errors=[];waves=[]
    for j,(value,vec) in enumerate(zip(vals,vecs.T)):
        vec=vec/max(abs(vec));amp=1e-4
        t=np.linspace(0,5,501)
        ans=solve_ivp(lambda t,x:model.rhs(x),(0,5),model.x+amp*vec.real,
                      method='DOP853',t_eval=t,rtol=1e-10,atol=1e-12,max_step=.025)
        assert ans.success
        exact=model.x[:,None]+amp*np.real(vec[:,None]*np.exp(value*t))
        error=max(np.max(abs(ans.y-exact),axis=0))/amp
        assert error<.005
        errors.append([j+1,value.imag/(2*np.pi),error])
        idx=4+np.argmax(abs(vec[4:8]))
        axs[j].plot(t,(ans.y[idx]-model.x[idx])/amp,label='Nonlinear')
        axs[j].plot(t,(exact[idx]-model.x[idx])/amp,'--',label='Modal prediction')
        axs[j].set(xlabel='Time (s)',title=f'Mode {j+1}: {value.imag/(2*np.pi):.3f} Hz')
        waves.extend(zip(np.full(len(t),j+1),t,ans.y[idx],exact[idx]))
    axs[0].set_ylabel('Normalized speed perturbation');axs[0].legend(fontsize=6)
    figure(fig,'detailed_modal_transient')
    table('detailed_modal_time_check',['mode','freq_Hz','max_state_error_relative_to_initial_scale'],errors)
    table('detailed_modal_waveforms',['mode','t_s','nonlinear_speed_pu','modal_speed_pu'],waves)
    # Match the original 1.00--1.08 s fault timing, with the detailed machine equations.
    logs=[];states=[];volts=[]
    x=model.x.copy()
    for a,b,fault in [(0,1,False),(1,1.08,True),(1.08,12,False)]:
        case.Y=kd.Y.copy()
        if fault:case.Y[7,7]+=1/.0001
        times=np.linspace(a,b,round((b-a)/.005)+1)
        sol=solve_ivp(lambda t,x:model.rhs(x),(a,b),x,method='DOP853',
                       t_eval=times,rtol=1e-8,atol=1e-10,max_step=.025)
        if not sol.success:raise RuntimeError(sol.message)
        x=sol.y[:,-1]
        count=len(times) if b==12 else len(times)-1
        for ti,xi in zip(times[:count],sol.y.T[:count]):
            logs.append(ti);states.append(xi);volts.append(model.electrical(xi,True)[-1])
    case.Y=kd.Y.copy()
    t=np.array(logs);x=np.array(states);v=np.array(volts)
    angle=np.rad2deg(x[:,:4]-model.x[:4])
    angle-=np.average(angle,axis=1,weights=kd.H)[:,None]
    tie=(v[:,6]*np.conj((v[:,6]-v[:,7])*2/(110*(.0001+.001j)))).real*100
    # Inertia-weighted (COI) frequency deviation, summarised in the JSON only:
    # external_run/analyze_results.py compares every column of this CSV.
    coi_df=np.average(x[:,4:8]*60,axis=1,weights=kd.H)
    table('detailed_fault',['t_s',*[f'G{k}_coi_deg' for k in range(1,5)],
          *[f'G{k}_df_Hz' for k in range(1,5)],'V8_pu','P7to8_MW'],
          np.column_stack([t,angle,x[:,4:8]*60,abs(v[:,7]),tie]))
    fig,ax=plt.subplots(2,1,figsize=(7.15,3.3),sharex=True)
    for j in range(4):ax[0].plot(t,angle[:,j],label=f'G{j+1}')
    ax[0].legend(ncol=4);ax[0].set_ylabel('COI angle change (deg)')
    ax[1].plot(t,tie);ax[1].set(xlabel='Time (s)',ylabel='Tie power (MW)')
    for a in ax:a.axvspan(1,1.08,color='grey',alpha=.12)
    figure(fig,'detailed_fault')
    (OUT/'detailed_transient_checks.json').write_text(json.dumps(dict(
        max_modal_state_error_ratio=max(r[2] for r in errors),
        minimum_fault_V8_pu=float(min(abs(v[:,7]))),
        prefault_max_speed_deviation_Hz=float(np.max(abs(x[t<1,4:8]))*60),
        final_COI_df_Hz=float(coi_df[-1]),
        COI_df_note='No governor and no damping term: post-fault COI frequency does not return to 60 Hz.'),indent=2))
if __name__=='__main__':main()
