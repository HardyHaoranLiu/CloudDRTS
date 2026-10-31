from pathlib import Path
import psutil
# The container lacks /proc/cpuinfo. Use explicit serial execution for scheduling.
psutil.cpu_count=lambda logical=True:1
import andes,numpy as np
ROOT=Path(__file__).resolve().parent
andes.config_logger(stream_level=20)
case=str(ROOT/'case/kundur_full.xlsx')
ss=andes.load(case,setup=False,no_undill=True,default_config=True,
              pycode_path=str(ROOT/'pycode'),output_path=str(ROOT/'results'),no_output=False)
ss.prepare(quick=True,models=['Bus','PQ','PV','Slack','Line','Area','GENROU','TGOV1','EXDC2','Toggler','Fault','BusFreq'],
           nomp=True,ncpu=1)
assert ss.setup()
print('PF',ss.PFlow.run())
print('EIG',ss.EIG.run())
print('VERSION',andes.__version__,'BUS_COUNT',ss.Bus.n)
print('EIGENVALUES')
for v in ss.EIG.mu:
 if v.imag>2 and v.imag<10: print(v, v.imag/(2*np.pi),-v.real/abs(v)*100)
print('MODEL COUNTS',[(n,m.n) for n,m in ss.models.items() if m.n])
print('EXPORTS',list((ROOT/'results').iterdir()))
