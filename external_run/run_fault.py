from pathlib import Path
import psutil
psutil.cpu_count=lambda logical=True:1
import andes,numpy as np,csv,json,argparse
parser=argparse.ArgumentParser()
parser.add_argument("--step",type=float,default=.001)
parser.add_argument("--tag",default="consistent")
args=parser.parse_args()
assert 0 < args.step <= .001
assert args.tag.replace("_", "").isalnum()
from matched_case import build,ROOT
andes.config_logger(stream_level=20)
ss=build('matched_fault_'+args.tag,True,fault_resistance=.0001)
assert ss.PFlow.run()
from event_consistency import install
event_records=install(ss)
ss.Fault.config.scale=1.0
ss.TDS.config.tf=12.
ss.TDS.config.tstep=args.step
ss.TDS.config.fixt=1
ss.TDS.config.shrinkt=0
ss.TDS.config.tol=1e-9
ss.TDS.config.max_iter=25
success=ss.TDS.run()
ts=ss.dae.ts;ts.unpack()
t=ts.t
angles=ts.x[:,ss.GENROU.delta.a]
angles=np.rad2deg(angles-angles[0])
angles-=np.average(angles,axis=1,weights=[6.5,6.5,6.175,6.175])[:,None]
speed=(ts.x[:,ss.GENROU.omega.a]-1)*60
v=ts.y[:,ss.Bus.v.a]*np.exp(1j*ts.y[:,ss.Bus.a.a])
tie=(v[:,6]*np.conj((v[:,6]-v[:,7])*2/(110*(.0001+.001j)))).real*100
data=np.column_stack([t,angles,speed,abs(v[:,7]),tie])
with (ROOT/f'results/matched_fault_{args.tag}_waveform.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['t_s',*[f'G{i}_coi_deg' for i in range(1,5)],
      *[f'G{i}_df_Hz' for i in range(1,5)],'V8_pu','P7to8_MW']);w.writerows(data)
(ROOT/f'results/fault_{args.tag}_run.json').write_text(json.dumps(dict(andes=andes.__version__,
    steps=len(t),initial_delta_rad=ts.x[0,ss.GENROU.delta.a].tolist(),end_time=float(t[-1]),busted=ss.TDS.busted,method=ss.TDS.config.method,
    events=event_records,fault_resistance_pu=.0001,requested_step=args.step,tolerance=ss.TDS.config.tol),indent=2))
print('FAULT FINISHED',len(t),t[-1],np.min(abs(v[:,7])),success)
assert success,ss.TDS.err_msg
