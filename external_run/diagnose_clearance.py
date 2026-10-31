import psutil
psutil.cpu_count=lambda logical=True:1
import andes,numpy as np
from matched_case import build
s=build('clearance_diagnostic',True)
andes.config_logger(stream_level=20)
assert s.PFlow.run()
s.TDS.config.tf=1.12;s.TDS.config.tstep=.001
s.TDS.config.tol=1e-9;s.TDS.config.shrinkt=0
s.TDS.config.honest=1
old=s.Fault.tc.callback
def instrument(mask):
    if not np.any(mask):return old(mask)
    print('CLEAR_BEFORE',s.dae.t,'uf',s.Fault.uf.v,'V8',s.dae.y[s.Bus.v.a[7]],'restore',s.Fault.config.restore,'mode',s.Fault.config.mode,'store8',s.Fault._vstore[7])
    v=old(mask)
    print('CLEAR_AFTER','uf',s.Fault.uf.v,'V8',s.dae.y[s.Bus.v.a[7]])
    return v
s.Fault.tc.callback=instrument
print('TDS',s.TDS.run())
s.dae.ts.unpack()
t=s.dae.ts.t;v=s.dae.ts.y[:,s.Bus.v.a[7]]
for target in [1.079,1.08,1.0801,1.081,1.082,1.09,1.1]:
 i=np.argmin(abs(t-target));print(t[i],v[i])
