"""Phase 2: validation of the nodal EMT solver against analytic solutions.

Cases (unit step source of 1 V applied at t = 0 unless stated otherwise):
  RC        R = 10 ohm, C = 100 uF
  RL        R = 10 ohm, L = 10 mH
  RLC       series R = 6 ohm, L = 10 mH, C = 10 uF (underdamped)
  line      lossless line Z = 50 ohm, travel time 0.5 ms, source resistance
            100 ohm, load 150 ohm (multiple reflections, lattice solution)
  switched  the RLC circuit with a 100 ohm shunt across C closed at 2 ms and
            opened at 5 ms (piecewise exact solution by matrix exponential)
  LC        undamped L = 10 mH, C = 10 uF, energy over 4000 steps
"""
import numpy as np
from scipy.linalg import expm

from common import plt, save_csv, save_json, save_figure, panel_labels
from emt import Circuit, Bergeron

STEPS_US = [100, 50, 25, 10, 5]


def rc(dt, t_end=0.005):
    c = Circuit(1, dt)
    c.resistor(-1, 1, 10.0)
    c.capacitor(1, 0, 100e-6, i0=0.1)          # state just after the step
    t = np.arange(round(t_end / dt) + 1) * dt
    v = np.r_[0.0, [c.step({-1: 1.0}, [0.0])[0] for _ in t[1:]]]
    return t, v, 1 - np.exp(-t / 1e-3)


def rl(dt, t_end=0.005):
    c = Circuit(1, dt)
    c.resistor(-1, 1, 10.0)
    c.inductor(1, 0, 0.01, v0=1.0)
    t = np.arange(round(t_end / dt) + 1) * dt
    i = np.r_[0.0, [(1 - c.step({-1: 1.0}, [0.0])[0]) / 10 for _ in t[1:]]]
    return t, i, 0.1 * (1 - np.exp(-t / 1e-3))


def rlc(dt, t_end=0.01):
    R, L, C = 6.0, 0.01, 10e-6
    alpha = R / (2 * L)
    wd = np.sqrt(1 / (L * C) - alpha ** 2)
    c = Circuit(2, dt)
    c.resistor(-1, 1, R)
    c.inductor(1, 2, L, v0=1.0)
    c.capacitor(2, 0, C)
    t = np.arange(round(t_end / dt) + 1) * dt
    v = np.r_[0.0, [c.step({-1: 1.0}, [0.0, 0.0])[1] for _ in t[1:]]]
    exact = 1 - np.exp(-alpha * t) * (np.cos(wd * t) + alpha / wd * np.sin(wd * t))
    return t, v, exact


def line(dt, t_end=0.006):
    z, tau = 50.0, 5e-4
    ln = Bergeron(z, tau, dt)
    c = Circuit(2, dt)
    c.resistor(-1, 1, 100.0)
    c.resistor(1, 0, z)            # line conductance 1/Z at the sending end
    c.resistor(2, 0, 150.0)
    c.resistor(2, 0, z)            # and at the receiving end
    t = np.arange(round(t_end / dt) + 1) * dt
    out = []
    for _ in t:
        hl, hr = ln.incoming()
        vl, vr = c.step({-1: 1.0}, [-hl, -hr])
        ln.publish(vl, vr, hl, hr)
        out.append(vr)
    # Lattice: reflection 1/3 at the source and 1/2 at the load; the receiving
    # voltage rises by 0.6*(5/6)*(1/6)**m at t = (2m+1) tau.
    n = round(tau / dt)
    arrivals = np.maximum(0, (np.arange(len(t)) - n) // (2 * n) + 1)
    return t, np.array(out), 0.6 * (1 - (1 / 6) ** arrivals)


def switched_exact(t, on=0.002, off=0.005, r_sh=100.0):
    R, L, C = 6.0, 0.01, 10e-6

    def segment(x, tt, shunt):
        A = np.array([[-R / L, -1 / L], [1 / C, -(1 / shunt if shunt else 0.0) / C]])
        steady = -np.linalg.solve(A, np.array([1 / L, 0.0]))
        return steady + expm(A * tt) @ (x - steady)

    x_on = segment(np.zeros(2), on, None)
    x_off = segment(x_on, off - on, r_sh)
    out = []
    for ti in t:
        if ti <= on:
            out.append(segment(np.zeros(2), ti, None))
        elif ti <= off:
            out.append(segment(x_on, ti - on, r_sh))
        else:
            out.append(segment(x_off, ti - off, None))
    return np.array(out)          # columns: inductor current, capacitor voltage


def switched(dt, t_end=0.008, on=0.002, off=0.005, resolve=True):
    """resolve=False skips the second solve at the switching instants."""
    c = Circuit(2, dt)
    c.resistor(-1, 1, 6.0)
    c.inductor(1, 2, 0.01, v0=1.0)
    c.capacitor(2, 0, 10e-6)
    c.switch(2, 0, 100.0)
    j_on, j_off = round(on / dt), round(off / dt)
    t = np.arange(round(t_end / dt) + 1) * dt
    wave, jump = [[0.0, 0.0]], 0.0
    for j in range(1, len(t)):
        if not resolve:
            v = c.step({-1: 1.0}, [0.0, 0.0], (j_on <= j < j_off,))
            wave.append([c.state[0][1], v[1]])
            continue
        v = c.step({-1: 1.0}, [0.0, 0.0], (j_on < j <= j_off,))
        if j in (j_on, j_off):
            before = np.array([c.state[0][1], c.state[1][0]])
            v = c.reinitialize({-1: 1.0}, [0.0, 0.0], (j == j_on,))
            jump = max(jump, np.max(abs(np.array([c.state[0][1], c.state[1][0]]) - before)))
        wave.append([c.state[0][1], v[1]])
    return t, np.array(wave), switched_exact(t), jump


def lc_energy(dt=5e-6, steps=4000):
    L, C = 0.01, 10e-6
    c = Circuit(1, dt)
    c.inductor(1, 0, L, v0=1.0)
    c.capacitor(1, 0, C, v0=1.0)
    worst = 0.0
    for _ in range(steps):
        v = c.step({}, [0.0])[0]
        i = c.state[0][1]
        worst = max(worst, abs((0.5 * L * i ** 2 + 0.5 * C * v ** 2) / (0.5 * C) - 1))
    return worst


def order(steps, errors):
    return float(np.polyfit(np.log(steps), np.log(errors), 1)[0])


def main():
    rows, errs = [], {}
    for name, fn, unit in [("RC", rc, "V"), ("RL", rl, "A"), ("RLC", rlc, "V")]:
        errs[name] = []
        for us in STEPS_US:
            _, y, exact = fn(us * 1e-6)
            e = float(np.max(abs(y - exact)))
            errs[name].append(e)
            rows.append([name, us, e, unit])
    errs["switched"] = []
    jumps = []
    for us in STEPS_US:
        _, y, exact, jump = switched(us * 1e-6)
        e = np.max(abs(y - exact), axis=0)
        errs["switched"].append(float(e[1]))
        jumps.append(jump)
        rows.append(["switched RLC, current", us, float(e[0]), "A"])
        rows.append(["switched RLC, voltage", us, float(e[1]), "V"])
    no_resolve = []
    for us in STEPS_US:
        _, y, exact, _ = switched(us * 1e-6, resolve=False)
        no_resolve.append(float(np.max(abs(y[:, 1] - exact[:, 1]))))
        rows.append(["switched RLC without re-solve, voltage", us, no_resolve[-1], "V"])
    line_err = []
    for us in (50, 5):
        _, y, exact = line(us * 1e-6)
        line_err.append(float(np.max(abs(y - exact))))
        rows.append(["line", us, line_err[-1], "V"])
    save_csv("phase2_validation", ["case", "dt_us", "max_abs_error", "unit"], rows)
    orders = {k: order(STEPS_US, v) for k, v in errs.items()}
    orders["switched_without_resolve"] = order(STEPS_US, no_resolve)
    energy = lc_energy()
    save_json("phase2_summary", dict(observed_order=orders, line_max_error_V=max(line_err),
                                     lc_max_relative_energy_change=energy,
                                     max_state_jump_at_switching=max(jumps)))
    print(orders, max(line_err), energy, max(jumps))

    fig, axes = plt.subplots(2, 2, figsize=(3.45, 3.2))
    axes = axes.ravel()
    for ax, (name, fn, ylabel) in zip(axes[:3], [("RC", rc, "RC: $v_C$ (V)"),
                                                  ("RL", rl, "RL: $i_L$ (A)"),
                                                  ("line", line, "Line: $v_2$ (V)")]):
        for us, ls in [(50, "-"), (5, "--")]:
            t, y, exact = fn(us * 1e-6)
            ax.plot(t * 1e3, y, ls, label=f"{us} µs")
        stride = max(1, len(t) // 60)
        ax.plot(t[::stride] * 1e3, exact[::stride], ".", color="k", ms=2.5, label="exact")
        ax.set(xlabel="Time (ms)", ylabel=ylabel)
    t, y, exact, _ = switched(5e-6)
    t50, y50, _, _ = switched(50e-6)
    axes[3].plot(t50 * 1e3, y50[:, 1], label="50 µs")
    axes[3].plot(t * 1e3, y[:, 1], "--", label="5 µs")
    axes[3].plot(t[::30] * 1e3, exact[::30, 1], ".", color="k", ms=2.5, label="exact")
    axes[3].axvspan(2, 5, color="0.9", lw=0)
    axes[3].set(xlabel="Time (ms)", ylabel="Switched RLC: $v_C$ (V)")
    axes[0].legend(fontsize=6, handlelength=1.5)
    panel_labels(axes)
    fig.tight_layout()
    save_figure(fig, "phase2_waveforms")

    fig, ax = plt.subplots(figsize=(3.45, 2.4))
    for (name, label), marker in zip([("RC", "RC step"), ("RL", "RL step"), ("RLC", "RLC step"),
                                      ("switched", "RLC, switched")], "osv^"):
        ax.loglog(STEPS_US, errs[name], marker + "-", label=f"{label} ({orders[name]:.2f})")
    ref = np.array(STEPS_US, float)
    ax.loglog(ref, errs["RLC"][-1] * (ref / ref[-1]) ** 2, ":", color="0.4", label="slope 2")
    ax.set(xlabel="Time step (µs)", ylabel="Maximum error")
    ax.set_ylim(top=2.0)
    ax.legend(fontsize=6, loc="upper left")
    save_figure(fig, "phase2_convergence")


if __name__ == "__main__":
    main()
