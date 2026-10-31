"""Build the 11-bus manually excited case with parameters matching our six-state model."""
from pathlib import Path
import psutil
psutil.cpu_count=lambda logical=True:1
import andes
ROOT=Path(__file__).resolve().parent

def build(output='matched',fault=False,fault_resistance=.0001):
    folder=ROOT/'results'/output;folder.mkdir(parents=True,exist_ok=True)
    ss=andes.System(name=output,default_config=True,no_undill=True,
        pycode_path=str(ROOT/'pycode'),output_path=str(folder),no_output=False)
    ss.config.mva=100;ss.config.freq=60
    ss.PFlow.config.tol=1e-11
    ss.PQ.config.pq2z=0
    ss.PQ.config.p2p=0;ss.PQ.config.p2i=1;ss.PQ.config.p2z=0
    ss.PQ.config.q2q=0;ss.PQ.config.q2i=0;ss.PQ.config.q2z=1
    for b in range(1,12):
        ss.add('Bus',idx=b,name=f'B{b}',Vn=20 if b<=4 else 230,
               v0=[1.03,1.01,1.03,1.01][b-1] if b<=4 else 1.,a0=0.)
    for i,(p,v) in enumerate(zip([7.,7.,7.19,7.],[1.03,1.01,1.03,1.01]),1):
        m='Slack' if i==1 else 'PV'
        ss.add(m,idx=f'G{i}',bus=i,Sn=900,Vn=20,p0=p,v0=v,
               qmin=-20,qmax=20,pmin=0,pmax=20,**({'a0':0.} if i==1 else {}))
        ss.add('GENROU',idx=i,name=f'G{i}',bus=i,gen=f'G{i}',Sn=900,Vn=20,fn=60,
               M=13 if i<=2 else 12.35,D=0,ra=.0025,xl=.2,xd=1.8,xq=1.7,
               xd1=.3,xq1=.55,xd2=.25,xq2=.25,Td10=8,Td20=.03,Tq10=.4,Tq20=.05,
               S10=0.,S12=1.)
    for j,(a,b) in enumerate([(1,5),(2,6),(3,11),(4,10)],1):
        ss.add('Line',idx=f'T{j}',bus1=a,bus2=b,Sn=100,Vn1=20,Vn2=230,
               r=0.,x=.15/9,b=0.,trans=1,tap=1.)
    for j,(a,b,km,count) in enumerate([(5,6,25,1),(6,7,10,1),(7,8,110,2),
                                       (8,9,110,2),(9,10,10,1),(10,11,25,1)],1):
        ss.add('Line',idx=f'L{j}',bus1=a,bus2=b,Sn=100,Vn1=230,Vn2=230,
               r=km*.0001/count,x=km*.001/count,b=km*.00175*count)
    for bus,p,qc in [(7,9.67,2.),(9,17.67,3.5)]:
        ss.add('PQ',idx=f'Load{bus}',bus=bus,Vn=230,p0=p,q0=1.)
        ss.add('Shunt',idx=f'C{bus}',bus=bus,Sn=100,Vn=230,b=qc,g=0.)
    if fault:ss.add('Fault',idx='F8',bus=8,tf=1.,tc=1.08,rf=fault_resistance,xf=0.)
    names=['Bus','PQ','PV','Slack','Line','Shunt','GENROU','Fault']
    ss.prepare(quick=True,models=names,nomp=True,ncpu=1)
    assert ss.setup()
    andes.io.json.write(ss,str(ROOT/'case'/(output+'.json')),overwrite=True)
    return ss

if __name__=='__main__':
    import numpy as np
    andes.config_logger(stream_level=20)
    s=build()
    assert s.PFlow.run()
    assert s.EIG.run()
    print('BUS VOLTAGES',s.Bus.v.v)
    for z in s.EIG.mu:
        if z.imag>2 and z.imag<10:print(z,z.imag/(2*np.pi),-z.real/abs(z)*100)
