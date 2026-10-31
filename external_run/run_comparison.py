"""Run official and matched ANDES cases and export numerical comparison evidence."""
from pathlib import Path
import psutil
psutil.cpu_count=lambda logical=True:1
import andes,numpy as np
import json,csv,hashlib,sys
from scipy.optimize import linear_sum_assignment
from matched_case import build,ROOT

OUT=ROOT/'results'
def save_csv(path,head,values):
    with path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(head);w.writerows(values)

def mode_rows(values):
    return [[z.real,z.imag,z.imag/(2*np.pi),-z.real/abs(z)*100] for z in values if z.imag>2 and z.imag<10]

def main():
    andes.config_logger(stream_level=20)
    ss=andes.load(str(ROOT/'case/kundur_full.xlsx'),setup=False,no_undill=True,
          default_config=True,pycode_path=str(ROOT/'pycode'),
          output_path=str(OUT/'official_193'))
    ss.prepare(quick=True,models=['Bus','PQ','PV','Slack','Line','Area','GENROU','TGOV1','EXDC2','Toggler'],
               nomp=True,ncpu=1)
    assert ss.setup() and ss.PFlow.run() and ss.EIG.run()
    official=mode_rows(ss.EIG.mu)
    save_csv(OUT/'official_193_modes.csv',['sigma_per_s','omega_rad_s','frequency_Hz','damping_pct'],official)
    np.save(OUT/'official_193_all_eigenvalues.npy',ss.EIG.mu)
    model=build('matched_modal')
    assert model.PFlow.run() and model.EIG.run()
    external=mode_rows(model.EIG.mu)
    np.save(OUT/'matched_all_eigenvalues.npy',model.EIG.mu)
    sys.path.insert(0,str(ROOT.parent/'CloudDRTS'))
    import kundur as kd
    from detailed_machine import Detailed
    from literature_checks import Case
    mine=Detailed(Case(kd.Y,kd.load,kd.V0,kd.H,1),'none')
    mine_vals,_=mine.modes()
    ext_vals=np.array([complex(r[0],r[1]) for r in external])
    ix,jx=linear_sum_assignment(abs(mine_vals[:,None]-ext_vals))
    comparisons=[]
    for i,j in zip(ix,jx):
        z=mine_vals[i];x=ext_vals[j]
        comparisons.append([i+1,z.real,x.real,z.imag/(2*np.pi),x.imag/(2*np.pi),
            -z.real/abs(z)*100,-x.real/abs(x)*100,abs(z-x),
            abs(z.imag-x.imag)/(2*np.pi)])
    save_csv(OUT/'matched_mode_comparison.csv',
      ['mode','local_sigma','andes_sigma','local_frequency_Hz','andes_frequency_Hz',
       'local_damping_pct','andes_damping_pct','eigenvalue_abs_difference_per_s','frequency_abs_difference_Hz'],comparisons)
    # Absolute phase reference is identical, with G1 at zero.
    va=model.Bus.v.v
    angle=model.Bus.a.v
    save_csv(OUT/'matched_powerflow.csv',['bus','andes_v_pu','local_v_pu','andes_angle_deg','local_angle_deg'],
       zip(range(1,12),va,abs(kd.V0),np.rad2deg(angle),np.angle(kd.V0,deg=True)))
    eigdiff=max(x[-2] for x in comparisons); fdiff=max(x[-1] for x in comparisons)
    assert fdiff<1e-5 and eigdiff<1e-4
    summary=dict(andes_version=andes.__version__,
      official_case_buses=ss.Bus.n, matched_case_buses=model.Bus.n,
      official_case_sha256=hashlib.sha256((ROOT/'case/kundur_full.xlsx').read_bytes()).hexdigest(),
      matched_case_sha256=hashlib.sha256((ROOT/'case/matched_modal.json').read_bytes()).hexdigest(),
      matched_frequency_max_abs_difference_Hz=float(fdiff),
      matched_eigenvalue_max_abs_difference_per_s=float(eigdiff),
      matched_voltage_max_abs_difference_pu=float(max(abs(va-abs(kd.V0)))),
      numerical_line_regularization='ANDES 1.9.3 adds 1e-8 to each branch r and x',
      small_eigenvalues=[dict(real=float(z.real),imag=float(z.imag)) for z in model.EIG.mu if abs(z)<1e-4],
      acceptance=dict(freq_abs_tol_Hz=1e-5,eigen_abs_tol_per_s=1e-4,passed=True),
      environment_serial_cpu_override='psutil.cpu_count -> 1 because /proc/cpuinfo is absent; model equations unmodified')
    (OUT/'modal_comparison_checks.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
