"""Phase 4: per-step cost of the Python solver against real-time budgets.

The Phase 3 circuit is run as two subsystems coupled by the Bergeron line
(60 Hz source, 5 us step, 10 000 steps including the two switching events).
Four implementations of the same computation are timed step by step:
  generic     emt.Circuit for each subsystem, np.linalg.solve every step
  cached      the same objects with the inverse matrix of each switch
              state computed once
  scalar      the two one-node subsystems written out with Python floats
  compiled    the scalar loop compiled with numba (if installed); only the
              mean step time can be measured
All four give the same waveform.
"""
import cProfile
import io
import platform
import pstats
import sys
import time

import numpy as np

from common import plt, save_csv, save_json, save_figure, RESULTS
from emt import Circuit, Bergeron
from phase3_partition import RS, RP, C, Z, TAU, RLOAD, RSW, VS, DT, T_ON, T_OFF

F = 60.0
STEPS = 10_000
REPEATS = 5
BUDGETS_US = (5, 50)


class CachedCircuit(Circuit):
    """Circuit that reuses the inverse matrix of each switch state."""

    def step(self, sources, injected, closed=()):
        key = tuple(closed)
        if not hasattr(self, "_inv"):
            self._inv = {}
        if key not in self._inv:
            self._inv[key] = np.linalg.inv(self._matrix(closed))
        rhs = self._rhs(sources, injected)
        for inc, _, _, _, h in self.storage:
            rhs -= inc * h
        v = self._inv[key] @ rhs
        for j, item in enumerate(self.storage):
            inc, known, g, sign, h = item
            drop = inc @ v + sum(s * sources[node] for node, s in known.items())
            self.state[j] = (drop, g * drop + h)
            item[-1] = sign * (2 * g * drop + h)
        return v


def run_objects(cls, timed=True):
    a, b = cls(1, DT), cls(1, DT)
    a.resistor(-1, 1, RS)
    a.resistor(1, 0, RP)
    a.capacitor(1, 0, C)
    a.resistor(1, 0, Z)
    b.resistor(1, 0, Z)
    b.resistor(1, 0, RLOAD)
    b.switch(1, 0, RSW)
    line = Bergeron(Z, TAU, DT)
    n_on, n_off = round(T_ON / DT), round(T_OFF / DT)
    w = 2 * np.pi * F
    out, ns = np.empty(STEPS), np.empty(STEPS)
    clock = time.perf_counter_ns
    for j in range(STEPS):
        start = clock()
        closed = n_on <= j < n_off
        event = j == n_on or j == n_off
        before = (not closed) if event else closed
        ha, hb = line.incoming()
        va = a.step({-1: VS * np.sin(w * j * DT)}, [-ha])[0]
        vb = b.step({}, [-hb], (before,))[0]
        if event:
            vb = b.reinitialize({}, [-hb], (closed,))[0]
        line.publish(va, vb, ha, hb)
        out[j] = vb
        ns[j] = clock() - start
    return out, ns


def run_scalar():
    import math
    ga_c = 2 * C / DT
    ga = 1 / RS + 1 / RP + ga_c + 1 / Z
    gb_open, gb_closed = 1 / Z + 1 / RLOAD, 1 / Z + 1 / RLOAD + 1 / RSW
    n = round(TAU / DT)
    lr, rl = [0.0] * n, [0.0] * n
    ptr, hc = 0, 0.0
    n_on, n_off = round(T_ON / DT), round(T_OFF / DT)
    w = 2 * math.pi * F
    out, ns = np.empty(STEPS), np.empty(STEPS)
    clock = time.perf_counter_ns
    for j in range(STEPS):
        start = clock()
        closed = n_on <= j < n_off
        ha, hb = rl[ptr], lr[ptr]
        va = (VS * math.sin(w * j * DT) / RS - hc - ha) / ga
        hc = -(2 * ga_c * va + hc)
        vb = -hb / (gb_closed if closed else gb_open)
        lr[ptr] = -2 * va / Z - ha
        rl[ptr] = -2 * vb / Z - hb
        ptr = (ptr + 1) % n
        out[j] = vb
        ns[j] = clock() - start
    return out, ns


def compiled_kernel():
    try:
        from numba import njit
    except ImportError:
        return None

    @njit(cache=False)
    def kernel(steps, dt, n, n_on, n_off, w):
        ga_c = 2 * C / dt
        ga = 1 / RS + 1 / RP + ga_c + 1 / Z
        gb_open, gb_closed = 1 / Z + 1 / RLOAD, 1 / Z + 1 / RLOAD + 1 / RSW
        lr, rl = np.zeros(n), np.zeros(n)
        out = np.empty(steps)
        ptr, hc = 0, 0.0
        for j in range(steps):
            closed = n_on <= j < n_off
            ha, hb = rl[ptr], lr[ptr]
            va = (VS * np.sin(w * j * dt) / RS - hc - ha) / ga
            hc = -(2 * ga_c * va + hc)
            vb = -hb / (gb_closed if closed else gb_open)
            lr[ptr] = -2 * va / Z - ha
            rl[ptr] = -2 * vb / Z - hb
            ptr = (ptr + 1) % n
            out[j] = vb
        return out

    args = (STEPS, DT, round(TAU / DT), round(T_ON / DT), round(T_OFF / DT), 2 * np.pi * F)
    kernel(*args)                                   # compile
    return lambda: kernel(*args)


def main():
    runners = {"generic": lambda: run_objects(Circuit), "cached": lambda: run_objects(CachedCircuit),
               "scalar": run_scalar}
    reference = runners["generic"]()[0]
    rows, samples, check = [], {}, {}
    for name, fn in runners.items():
        fn()                                        # warm-up
        times = []
        for rep in range(REPEATS):
            out, ns = fn()
            times.append(ns / 1000)
            check[name] = float(np.max(abs(out - reference)))
        us = np.concatenate(times)
        samples[name] = us
        rows.append([name, *np.percentile(us, [50, 99]), us.max(),
                     *[100 * np.mean(us > b) for b in BUDGETS_US]])
    kernel = compiled_kernel()
    if kernel is not None:
        check["compiled"] = float(np.max(abs(kernel() - reference)))
        mean = []
        for _ in range(REPEATS):
            start = time.perf_counter_ns()
            kernel()
            mean.append((time.perf_counter_ns() - start) / 1000 / STEPS)
        rows.append(["compiled (mean only)", min(mean), np.nan, np.nan, np.nan, np.nan])
    if max(check.values()) > 1e-9:
        raise RuntimeError(f"Implementations disagree: {check}")
    save_csv("phase4_timing", ["implementation", "median_us", "p99_us", "max_us",
                               "over_5us_pct", "over_50us_pct"], rows)

    profiler = cProfile.Profile()
    profiler.runcall(run_objects, Circuit)
    buf = io.StringIO()
    pstats.Stats(profiler, stream=buf).strip_dirs().sort_stats("tottime").print_stats(12)
    (RESULTS / "phase4_profile.txt").write_text(buf.getvalue())

    import scipy
    save_json("phase4_environment", dict(python=sys.version.split()[0], numpy=np.__version__,
                                         scipy=scipy.__version__, machine=platform.machine(),
                                         max_waveform_difference=check))

    fig, ax = plt.subplots(figsize=(3.45, 2.4))
    for name, c in zip(samples, ["C0", "C1", "C2"]):
        x = np.sort(samples[name])
        share = 1 - np.arange(len(x)) / len(x)      # share of steps taking at least x
        ax.loglog(x, share, color=c, label=name)
    if kernel is not None:
        ax.axvline(rows[-1][1], color="C3", lw=1, label="compiled (mean)")
    for b in BUDGETS_US:
        ax.axvline(b, color="k", ls="--", lw=0.8)
        ax.text(b * 1.08, 2e-5, f"{b} µs", fontsize=6.5)
    ax.set(xlabel="Time per step (µs)", ylabel="Share of steps taking longer",
           ylim=(1e-5, 1.5))
    ax.legend(loc="lower left", fontsize=6)
    save_figure(fig, "phase4_timing")
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
