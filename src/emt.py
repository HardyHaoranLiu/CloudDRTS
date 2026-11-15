"""Nodal EMT solver with trapezoidal companion models (Dommel, 1969).

Node 0 is ground. Negative node numbers are nodes with a prescribed voltage
(ideal voltage sources); their values are passed to step() in a dict.
An inductor or capacitor is replaced by a conductance g in parallel with a
history current h, i = g v + h, with v and i taken from node a to node b.
"""
from collections import deque
from numbers import Integral

import numpy as np


class Circuit:
    def __init__(self, n, dt):
        if not isinstance(n, Integral) or n < 1 or not (np.isfinite(dt) and dt > 0):
            raise ValueError("Use at least one node and a positive time step.")
        self.n, self.dt = n, dt
        self.G = np.zeros((n, n))
        self.fixed = {}       # source node -> right-hand-side column
        self.storage = []     # [incidence, source incidence, g, sign, history]
        self.state = []       # (voltage, current) of each storage element
        self.switches = []    # conductance stamp of each switch when closed

    def _branch(self, a, b, g):
        inc = np.zeros(self.n)
        known = {}
        for node, sign in ((a, 1.0), (b, -1.0)):
            if node > 0:
                inc[node - 1] += sign
            elif node < 0:
                known[node] = known.get(node, 0.0) + sign
        self.G += g * np.outer(inc, inc)
        for node, sign in known.items():
            self.fixed[node] = self.fixed.get(node, 0.0) - inc * g * sign
        return inc, known

    def resistor(self, a, b, r):
        if not (np.isfinite(r) and r > 0):
            raise ValueError("Resistance must be positive.")
        self._branch(a, b, 1 / r)

    def _storage(self, a, b, g, sign, v0, i0):
        inc, known = self._branch(a, b, g)
        self.storage.append([inc, known, g, sign, sign * (i0 + g * v0)])
        self.state.append((float(v0), float(i0)))

    def inductor(self, a, b, L, v0=0.0, i0=0.0):
        """i(t) = dt/(2L) v(t) + [i(t-dt) + dt/(2L) v(t-dt)].
        v0 and i0 are the voltage and current one step before the first step."""
        if not (np.isfinite(L) and L > 0):
            raise ValueError("Inductance must be positive.")
        self._storage(a, b, self.dt / (2 * L), 1, v0, i0)

    def capacitor(self, a, b, C, v0=0.0, i0=0.0):
        """i(t) = 2C/dt v(t) - [2C/dt v(t-dt) + i(t-dt)]."""
        if not (np.isfinite(C) and C > 0):
            raise ValueError("Capacitance must be positive.")
        self._storage(a, b, 2 * C / self.dt, -1, v0, i0)

    def switch(self, a, b, r_on):
        """Ideal open switch, r_on when closed. Its state is given to step()."""
        if a < 0 or b < 0 or not (np.isfinite(r_on) and r_on > 0):
            raise ValueError("Connect switches between unknown or ground nodes; r_on > 0.")
        inc, _ = self._branch(a, b, 0.0)
        self.switches.append(np.outer(inc, inc) / r_on)

    def _matrix(self, closed):
        matrix = self.G
        for on, stamp in zip(closed, self.switches, strict=True):
            if on:
                matrix = matrix + stamp
        return matrix

    def _rhs(self, sources, injected):
        rhs = np.array(injected, dtype=float, copy=True)
        for node, col in self.fixed.items():
            rhs += col * sources[node]
        return rhs

    def step(self, sources, injected, closed=()):
        """Advance one step. injected: external currents into the unknown nodes."""
        rhs = self._rhs(sources, injected)
        for inc, _, _, _, h in self.storage:
            rhs -= inc * h
        v = np.linalg.solve(self._matrix(closed), rhs)
        for j, item in enumerate(self.storage):
            inc, known, g, sign, h = item
            drop = inc @ v + sum(s * sources[node] for node, s in known.items())
            self.state[j] = (drop, g * drop + h)
            item[-1] = sign * (2 * g * drop + h)
        return v

    def reinitialize(self, sources, injected, closed=()):
        """Re-solve the present instant after a switching event.

        Call step() with the old switch states up to the event time, then this
        method with the new states at the same time. Capacitor voltages and
        inductor currents are held; the other quantities jump, and the
        histories restart from the new values. Capacitor constraints must be
        independent (no loop of capacitors and sources).
        """
        matrix = self._matrix(closed).copy()
        rhs = self._rhs(sources, injected)
        constraints, values = [], []
        for item, (drop, current) in zip(self.storage, self.state):
            inc, known, g, sign, _ = item
            source_drop = sum(f * sources[node] for node, f in known.items())
            matrix -= g * np.outer(inc, inc)
            rhs += inc * g * source_drop
            if sign == 1:
                rhs -= inc * current
            else:
                constraints.append(inc)
                values.append(drop - source_drop)
        cap = np.column_stack(constraints) if constraints else np.empty((self.n, 0))
        m = cap.shape[1]
        system = np.block([[matrix, cap], [cap.T, np.zeros((m, m))]])
        try:
            answer = np.linalg.solve(system, np.r_[rhs, values])
        except np.linalg.LinAlgError as exc:
            raise ValueError("Capacitor constraints are not independent.") from exc
        v = answer[:self.n]
        k = self.n
        for j, item in enumerate(self.storage):
            inc, known, g, sign, _ = item
            drop = inc @ v + sum(f * sources[node] for node, f in known.items())
            current = self.state[j][1]
            if sign == -1:
                current = answer[k]
                k += 1
            self.state[j] = (drop, current)
            item[-1] = sign * (current + g * drop)
        return v


class Bergeron:
    """Lossless line (Bergeron model).

    With currents positive into the line at both ends,
        i_k(t) = v_k(t)/Z + h_k(t),   h_k(t) = -2 v_m(t-T)/Z - h_m(t-T),
    where T = travel time + extra*dt. Each end passes one number per step to
    the other end: its outgoing history -2 v/Z - h. extra > 0 adds an
    interface delay of extra steps in each direction.
    """

    def __init__(self, z, tau, dt, extra=0):
        if not all(np.isfinite(x) and x > 0 for x in (z, tau, dt)):
            raise ValueError("Use positive line parameters and time step.")
        if not isinstance(extra, Integral) or extra < 0:
            raise ValueError("The added delay is a nonnegative number of steps.")
        n = round(tau / dt)
        if n < 1 or not np.isclose(n * dt, tau, rtol=1e-12, atol=1e-15):
            raise ValueError("The travel time must be a whole number of steps.")
        self.z = z
        self.delay = n + extra
        self.lr = deque([0.0] * self.delay)   # left end -> right end
        self.rl = deque([0.0] * self.delay)   # right end -> left end

    def preload(self, from_left, from_right):
        """Past outgoing histories of each end, oldest first."""
        if not len(from_left) == len(from_right) == self.delay:
            raise ValueError("Give one value per step of delay.")
        self.lr = deque(float(x) for x in from_left)
        self.rl = deque(float(x) for x in from_right)

    def incoming(self):
        """Histories (h_left, h_right) for the present step."""
        return self.rl[0], self.lr[0]

    def publish(self, vl, vr, hl, hr):
        self.lr.popleft()
        self.rl.popleft()
        self.lr.append(-2 * vl / self.z - hl)
        self.rl.append(-2 * vr / self.z - hr)
