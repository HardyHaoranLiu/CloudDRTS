"""Linear nodal EMT, fixed-potential ports and trapezoidal Norton companions."""
from collections import deque
from numbers import Integral
import numpy as np

class Circuit:
    """0 is ground. Negative node IDs are ideal prescribed-voltage source nodes."""
    def __init__(self, n, dt):
        if not isinstance(n, Integral) or n < 1 or not np.isfinite(dt) or dt <= 0:
            raise ValueError("Node count and time step must be positive.")
        self.n, self.dt = n, dt
        self.G = np.zeros((n,n))
        self.fixed = {}
        self.storage = []
        self.storage_state = []
        self.switches = []

    def branch(self, a, b, g):
        inc = np.zeros(self.n)
        known = {}
        for node, sign in [(a,1.),(b,-1.)]:
            if node > 0: inc[node-1] += sign
            elif node < 0: known[node] = known.get(node, 0.) + sign
        self.G += g*np.outer(inc,inc)
        for node, sign in known.items():
            self.fixed[node] = self.fixed.get(node,0.)-inc*g*sign
        return inc, known

    def resistor(self,a,b,r):
        if not np.isfinite(r) or r <= 0:
            raise ValueError("Resistance must be positive and finite.")
        self.branch(a,b,1/r)

    def reactive(self,kind,a,b,value,v0=0.,i0=0.):
        if kind not in {"L", "C"} or not np.isfinite(value) or value <= 0:
            raise ValueError("Use kind L or C and a positive finite value.")
        # i = g*v + history, with v and i oriented a -> b.
        g = self.dt/(2*value) if kind=='L' else 2*value/self.dt
        inc,known = self.branch(a,b,g)
        sign = 1 if kind=='L' else -1
        self.storage.append([inc,known,g,sign,sign*(i0+g*v0)])
        self.storage_state.append((float(v0),float(i0)))

    def switch(self,a,b,r_on):
        if a < 0 or b < 0 or not np.isfinite(r_on) or r_on <= 0:
            raise ValueError("Switch requires unknown or ground nodes and positive resistance.")
        inc,known=self.branch(a,b,0)
        self.switches.append(np.outer(inc,inc)/r_on)

    def step(self, sources, injected, states=()):
        rhs=np.array(injected,dtype=float,copy=True)
        for node,col in self.fixed.items(): rhs+=col*sources[node]
        for inc,known,g,sign,h in self.storage: rhs-=inc*h
        matrix=self.G
        for closed,stamp in zip(states,self.switches,strict=True):
            if closed: matrix=matrix+stamp
        v=np.linalg.solve(matrix,rhs)
        for j,item in enumerate(self.storage):
            inc,known,g,sign,h=item
            drop=inc@v+sum(s*sources[n] for n,s in known.items())
            self.storage_state[j]=(drop,g*drop+h)
            item[-1]=sign*(2*g*drop+h)
        return v

    def reinitialize(self, sources, injected, states=()):
        """Resolve a right-hand event limit without advancing time.

        Call step with the old topology up to the exact event timestamp, then
        call this method with the new topology at that same timestamp.
        Capacitor voltage and inductor current remain continuous. Capacitor
        voltage constraints must be independent and consistent.
        """
        matrix=self.G.copy()
        rhs=np.asarray(injected,dtype=float).copy()
        for node,col in self.fixed.items():
            rhs+=col*sources[node]
        for closed,stamp in zip(states,self.switches,strict=True):
            if closed:matrix+=stamp
        constraints=[];values=[]
        for item,(drop,current) in zip(self.storage,self.storage_state):
            inc,known,g,sign,h=item
            source_drop=sum(factor*sources[node] for node,factor in known.items())
            matrix-=g*np.outer(inc,inc)
            rhs+=inc*g*source_drop
            if sign==1:
                rhs-=inc*current
            else:
                constraints.append(inc)
                values.append(drop-source_drop)
        cap=np.column_stack(constraints) if constraints else np.empty((self.n,0))
        m=cap.shape[1]
        system=np.block([[matrix,cap],[cap.T,np.zeros((m,m))]])
        try:
            answer=np.linalg.solve(system,np.r_[rhs,values])
        except np.linalg.LinAlgError as exc:
            raise ValueError("Event constraints must be independent and consistent.") from exc
        v=answer[:self.n]
        cap_index=0
        for j,item in enumerate(self.storage):
            inc,known,g,sign,h=item
            drop=inc@v+sum(factor*sources[node] for node,factor in known.items())
            current=self.storage_state[j][1]
            if sign==-1:
                current=answer[self.n+cap_index]
                cap_index+=1
            self.storage_state[j]=(drop,current)
            item[-1]=sign*(current+g*drop)
        return v

class Bergeron:
    """Current into either end is v/Z + history. Exchange outgoing waves only."""
    def __init__(self, z, tau, dt, k=0):
        if (not all(np.isfinite(v) and v > 0 for v in (z, tau, dt))
                or not isinstance(k, Integral) or k < 0):
            raise ValueError("Use positive finite line parameters and integer nonnegative k.")
        n=round(tau/dt)
        if n<1 or not np.isclose(n*dt,tau,rtol=1e-12,atol=1e-15) or k<0:
            raise ValueError('Use positive integer physical delay and nonnegative added steps.')
        self.z=z
        self.lr=deque([0.]*(n+k)); self.rl=self.lr.copy()

    def incoming(self):
        return self.rl[0],self.lr[0]

    def publish(self, vl, vr, hl, hr):
        self.lr.popleft(); self.rl.popleft()
        self.lr.append(-2*vl/self.z-hl); self.rl.append(-2*vr/self.z-hr)
