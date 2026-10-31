from pathlib import Path
import psutil
psutil.cpu_count=lambda logical=True:1
import andes,numpy as np,json
ROOT=Path(__file__).resolve().parent
andes.config_logger(stream_level=20)
out=ROOT/'results/historical';out.mkdir(parents=True,exist_ok=True)
case=andes.get_case('kundur/kundur_full.xlsx')
ss=andes.load(case,setup=False,no_undill=True,default_config=True,
     pycode_path=str(ROOT/'historical/pycode'),output_path=str(out))
ss.prepare(quick=True,models=['Bus','PQ','PV','Slack','Line','Area','GENROU','TGOV1','EXDC2','Toggler'],
           nomp=True,ncpu=1)
assert ss.setup()
assert ss.PFlow.run()
assert ss.EIG.run()
vals=[v for v in ss.EIG.mu if v.imag>2 and v.imag<10]
record=dict(andes=andes.__version__,case=case,buses=ss.Bus.n,
  modes=[dict(sigma=v.real,omega=v.imag,freq=v.imag/(2*np.pi),damping=-v.real/abs(v)*100) for v in vals])
(out/'observed.json').write_text(json.dumps(record,indent=2))
print(record)
