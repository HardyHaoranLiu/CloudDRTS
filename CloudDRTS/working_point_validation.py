"""Power balance and discrete-sinusoid convergence of the derived EMT network."""
import json
from common import np, plt, table, figure, OUT
from grid_cuts import model, legacy_model, operating_point, PORTS,ZC,NPHYS,DT,W,simulate
from kundur import V0

def discrete_solution(dt,parameters):
    B,L,R,C,G,E=parameters
    wd=2*np.tan(W*dt/2)/dt
    adm=1/(R+1j*wd*L)
    Y=(B*adm)@B.T+np.diag(G+1j*wd*C)
    phi=W*NPHYS*DT
    for j in (0,2):
        a,b=PORTS[j:j+2]
        Y[a,a]+=-1j/(ZC*np.tan(phi));Y[b,b]+=-1j/(ZC*np.tan(phi))
        Y[a,b]+=1j/(ZC*np.sin(phi));Y[b,a]+=1j/(ZC*np.sin(phi))
    v=np.linalg.solve(Y,B@(adm*E))
    return v

def main():
    B,L,R,C,G,E,V,Y,pf_error=operating_point()
    total=(V*np.conj(Y@V)).sum()
    line_currents=(B[:,4:12].T@V)/(R[4:12]+1j*W*L[4:12])
    losses=float(np.sum(R[4:12]*abs(line_currents)**2))
    assert abs(total.real-losses)<1e-10
    original=discrete_solution(DT,legacy_model())
    data=[]; trace=[]
    for dt in [10e-6,5e-6,2.5e-6]:
        v=discrete_solution(dt,model())
        absolute=max(abs(v-V));mag=max(abs(abs(v)-abs(V)))
        data.append([dt*1e6,absolute,mag])
        trace.extend(zip(np.full(11,dt*1e6),range(1,12),abs(V),abs(v),abs(v-V)))
        assert absolute<1e-5
    order=float(np.polyfit(np.log(np.array(data)[:,0]),np.log(np.array(data)[:,1]),1)[0])
    assert 1.99<order<2.01
    table('working_point_refinement',['dt_us','max_phasor_error_pu','max_magnitude_error_pu'],data)
    table('working_point_buses',['dt_us','bus','physical_v_pu','discrete_v_pu','phasor_error_pu'],trace)
    table('working_point_change',['bus','original_phasor_v_pu','legacy_emt_v_pu','new_physical_target_v_pu'],
          zip(range(1,12),abs(V0),abs(original),abs(V)))
    summary=dict(power_flow_max_residual=pf_error,
                 generated_MW=float(np.sum((V*np.conj(Y@V)).real[:4])*100),
                 active_load_MW=2734,branch_loss_MW=losses*100,
                 power_balance_error_MW=abs(total.real-losses)*100,
                 discrete_phasor_error_order=order,
                 note='Operating point belongs to represented 390-us lossless long-line network; physical load powers retained.')
    (OUT/'working_point_checks.json').write_text(json.dumps(summary,indent=2))
    fig,ax=plt.subplots(1,2,figsize=(7.15,2.7))
    ax[0].plot(range(1,12),abs(V0),'o-',label='Original phasor model')
    ax[0].plot(range(1,12),abs(V),'s--',label='Derived EMT target')
    ax[0].legend(fontsize=6);ax[0].set(xlabel='Bus',ylabel='RMS voltage (p.u.)')
    ax[1].loglog(np.array(data)[:,0],np.array(data)[:,1],'o-',label=f'Order {order:.2f}')
    ax[1].legend();ax[1].set(xlabel='Time step (µs)',ylabel='Max discrete phasor error (p.u.)')
    figure(fig,'working_point_validation')
    print(summary)
if __name__=='__main__':main()
