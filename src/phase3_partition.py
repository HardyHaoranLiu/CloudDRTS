"""Phase 3: split the test circuit at its Bergeron line, delay the interface,
and replace the sending side with Thevenin equivalents.

Test circuit
  A (sending side):   source vs(t) = 100 sin(2 pi f t) V behind Rs = 25 ohm,
                      shunt Rp = 50 ohm and C = 10 uF at the line end
  line:               lossless, Z = 50 ohm, travel time 250 us
  B (receiving side): load 80 ohm and an 18 ohm shunt that closes at 20 ms
                      and opens at 35 ms

Variants
  reference  one solver for the whole circuit
  split      solver A and solver B exchange only the line histories, each
             direction delayed by k extra steps
  static     A replaced by its resistive Thevenin equivalent (C ignored)
  phasor     A replaced by its Thevenin equivalent at the source frequency,
             realized as a source behind a series R-C
Every run starts from its own discrete-time sinusoidal steady state, so the
error before the switching event is the steady-state error only.
"""
import numpy as np

from common import plt, save_csv, save_json, save_figure, panel_labels, rel_rms
from emt import Circuit, Bergeron

RS, RP, C = 25.0, 50.0, 10e-6
Z, TAU = 50.0, 250e-6
RLOAD, RSW = 80.0, 18.0
VS = 100.0
DT = 5e-6
T_ON, T_OFF, T_END = 0.020, 0.035, 0.050
DELAYS = [0, 1, 2, 4, 10, 20, 40, 80]          # added steps, 0 ... 400 us
FREQS = [60.0, 400.0, 1200.0]
VARIANTS = ["split", "static", "phasor"]
NAMES = {"split": "Full split (histories only)", "static": "Resistive Thevenin",
         "phasor": "Phasor Thevenin"}


def equivalents(f):
    """Thevenin equivalents of side A seen from the line end."""
    w = 2 * np.pi * f
    zth = 1 / (1 / RS + 1 / RP + 1j * w * C)
    return dict(static=dict(gain=RP / (RS + RP), r=RS * RP / (RS + RP)),
                phasor=dict(gain=zth / RS, r=zth.real, c=-1 / (w * zth.imag)))


def steady_state(variant, k, f, dt=None):
    """Phasors (peak values) with the switch open. With dt, the discrete-time
    steady state of the solver (trapezoidal frequency warping); without dt,
    the continuous-time circuit with the delay k given in seconds."""
    w = 2 * np.pi * f
    if dt is None:
        wd, phi = w, w * (TAU + k)
    else:
        wd, phi = 2 * np.tan(w * dt / 2) / dt, w * dt * (round(TAU / dt) + k)
    ys, ym = -1j / (Z * np.tan(phi)), 1j / (Z * np.sin(phi))
    vs = -1j * VS
    eq = equivalents(f)
    if variant in ("reference", "split"):
        Y = [[1 / RS + 1 / RP + 1j * wd * C + ys, ym], [ym, 1 / RLOAD + ys]]
        I = [vs / RS, 0]
    elif variant == "static":
        s = eq["static"]
        Y = [[1 / s["r"] + ys, ym], [ym, 1 / RLOAD + ys]]
        I = [s["gain"] * vs / s["r"], 0]
    else:
        p = eq["phasor"]
        yc = 1j * wd * p["c"]
        Y = [[1 / p["r"] + yc, -yc, 0], [-yc, yc + ys, ym], [0, ym, 1 / RLOAD + ys]]
        I = [p["gain"] * vs / p["r"], 0, 0]
    V = np.linalg.solve(np.array(Y, complex), np.array(I, complex))
    va, vb = V[-2], V[-1]
    ia, ib = ys * va + ym * vb, ys * vb + ym * va        # currents into the line
    out = dict(va=va, vb=vb, ia=ia, ib=ib, emit_a=-(va / Z + ia), emit_b=-(vb / Z + ib))
    if variant in ("reference", "split"):
        out["cap"] = (va, 1j * wd * C * va)
    elif variant == "phasor":
        out["cap"] = (V[0] - va, 1j * wd * eq["phasor"]["c"] * (V[0] - va))
    return out


def wave(phasor, w, t):
    return np.real(phasor * np.exp(1j * w * np.asarray(t)))


def simulate(variant, k=0, f=60.0, dt=DT):
    """Returns t, v_B and i_B (receiving end, current into the line)."""
    w = 2 * np.pi * f
    ss = steady_state(variant, k, f, dt)
    line = Bergeron(Z, TAU, dt, k)
    past = np.arange(-line.delay, 0) * dt
    line.preload(wave(ss["emit_a"], w, past), wave(ss["emit_b"], w, past))
    cap = ss.get("cap")
    v0 = (wave(cap[0], w, -dt), wave(cap[1], w, -dt)) if cap else None
    eq = equivalents(f)

    if variant == "reference":
        a = b = Circuit(2, dt)
        a.resistor(-1, 1, RS)
        a.resistor(1, 0, RP)
        a.capacitor(1, 0, C, *v0)
        a.resistor(1, 0, Z)
        a.resistor(2, 0, Z)
        a.resistor(2, 0, RLOAD)
        a.switch(2, 0, RSW)
    else:
        if variant == "split":
            a = Circuit(1, dt)
            a.resistor(-1, 1, RS)
            a.resistor(1, 0, RP)
            a.capacitor(1, 0, C, *v0)
            a.resistor(1, 0, Z)
        elif variant == "static":
            a = Circuit(1, dt)
            a.resistor(-1, 1, eq["static"]["r"])
            a.resistor(1, 0, Z)
        else:
            a = Circuit(2, dt)                       # node 1 internal, node 2 line end
            a.resistor(-1, 1, eq["phasor"]["r"])
            a.capacitor(1, 2, eq["phasor"]["c"], *v0)
            a.resistor(2, 0, Z)
        b = Circuit(1, dt)
        b.resistor(1, 0, Z)
        b.resistor(1, 0, RLOAD)
        b.switch(1, 0, RSW)

    if variant == "static":
        source = lambda ti: eq["static"]["gain"] * VS * np.sin(w * ti)
    elif variant == "phasor":
        g = eq["phasor"]["gain"]
        source = lambda ti: abs(g) * VS * np.sin(w * ti + np.angle(g))
    else:
        source = lambda ti: VS * np.sin(w * ti)

    n_on, n_off = round(T_ON / dt), round(T_OFF / dt)
    t = np.arange(round(T_END / dt) + 1) * dt
    vb, ib = np.empty(len(t)), np.empty(len(t))
    for j, tj in enumerate(t):
        closed = n_on <= j < n_off
        event = j in (n_on, n_off)
        before = (not closed) if event else closed
        ha, hb = line.incoming()
        src = {-1: source(tj)}
        if variant == "reference":
            va, vbj = a.step(src, [-ha, -hb], (before,))
            if event:
                va, vbj = a.reinitialize(src, [-ha, -hb], (closed,))
        else:
            va = a.step(src, [-ha] if a.n == 1 else [0.0, -ha])[-1]
            vbj = b.step({}, [-hb], (before,))[0]
            if event:
                vbj = b.reinitialize({}, [-hb], (closed,))[0]
        line.publish(va, vbj, ha, hb)
        vb[j], ib[j] = vbj, vbj / Z + hb
    return t, vb, ib


def fit_phasor(t, y, f):
    basis = np.column_stack([np.cos(2 * np.pi * f * t), -np.sin(2 * np.pi * f * t)])
    re, im = np.linalg.lstsq(basis, y, rcond=None)[0]
    return re + 1j * im


def steady_error(variant, f, delays_s):
    """Continuous-time steady-state error of v_B (percent) against the reference."""
    ref = steady_state("reference", 0.0, f)["vb"]
    return np.array([100 * abs(steady_state(variant, d, f)["vb"] - ref) / abs(ref) for d in delays_s])


def resistor_path(ax, x0, y0, x1, y1, n=6, amp=0.11):
    s = np.linspace(0, 1, 2 * n + 1)
    off = np.r_[0, np.tile([amp, -amp], n)[:2 * n - 1], 0]
    dx, dy = x1 - x0, y1 - y0
    length = np.hypot(dx, dy)
    px, py = -dy / length, dx / length
    ax.plot(x0 + s * dx + off * px, y0 + s * dy + off * py, color="k", lw=1)


def circuit_figure():
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    k = dict(color="k", lw=1)
    ax.plot([0, 0], [0, 0.68], **k)
    ax.plot([0, 0], [1.32, 2], **k)
    ax.add_patch(plt.Circle((0, 1), 0.32, fill=False, lw=1))
    ax.text(0, 1, "~", ha="center", va="center", fontsize=13)
    ax.text(-0.42, 1, "$v_s(t)$", ha="right", va="center")
    ax.plot([0, 0.35], [2, 2], **k)
    resistor_path(ax, 0.35, 2, 1.15, 2)
    ax.text(0.75, 2.2, "$R_s$ 25 Ω", ha="center")
    ax.plot([1.15, 3.0], [2, 2], **k)
    ax.plot([1.6, 1.6], [2, 1.45], **k)
    resistor_path(ax, 1.6, 1.45, 1.6, 0.55)
    ax.plot([1.6, 1.6], [0.55, 0], **k)
    ax.text(1.78, 1.0, "$R_p$\n50 Ω", va="center")
    ax.plot([2.45, 2.45], [2, 1.08], **k)
    ax.plot([2.3, 2.6], [1.08, 1.08], **k)
    ax.plot([2.3, 2.6], [0.92, 0.92], **k)
    ax.plot([2.45, 2.45], [0.92, 0], **k)
    ax.text(2.65, 1.0, "$C$\n10 µF", va="center")
    ax.plot([0, 2.45], [0, 0], **k)
    ax.add_patch(plt.Rectangle((3.0, 1.62), 3.0, 0.76, fill=False, lw=1, zorder=3))
    ax.text(4.5, 2.0, "lossless line\nZ = 50 Ω, τ = 250 µs", ha="center", va="center", fontsize=7.5,
            linespacing=1.4, zorder=4, bbox=dict(fc="white", ec="none", pad=1.5))
    ax.plot([6.0, 8.1], [2, 2], **k)
    ax.plot([7.0, 7.0], [2, 1.45], **k)
    resistor_path(ax, 7.0, 1.45, 7.0, 0.55)
    ax.plot([7.0, 7.0], [0.55, 0], **k)
    ax.text(7.18, 1.0, "$R_L$\n80 Ω", va="center")
    ax.plot([8.1, 8.1], [2, 1.75], **k)
    ax.plot([8.1, 8.33], [1.75, 1.42], **k)
    ax.plot([8.1, 8.1], [1.38, 1.2], **k)
    resistor_path(ax, 8.1, 1.2, 8.1, 0.45)
    ax.plot([8.1, 8.1], [0.45, 0], **k)
    ax.text(8.32, 0.85, "$R_{sw}$ 18 Ω\nclosed\n20–35 ms", va="center", fontsize=7)
    ax.plot([7.0, 8.1], [0, 0], **k)
    for x in (1.2, 7.55):
        for i, half in enumerate([0.22, 0.14, 0.06]):
            ax.plot([x - half, x + half], [-0.12 - 0.08 * i] * 2, **k)
        ax.plot([x, x], [0, -0.12], **k)
    ax.plot([4.5, 4.5], [-0.45, 3.35], color="0.4", lw=1, ls="--", zorder=1)
    ax.text(1.2, 3.05, "Solver A", ha="center", fontsize=8.5, weight="bold")
    ax.text(7.55, 3.05, "Solver B", ha="center", fontsize=8.5, weight="bold")
    ax.annotate("", xy=(5.9, 2.6), xytext=(3.1, 2.6), arrowprops=dict(arrowstyle="-|>", lw=0.9))
    ax.annotate("", xy=(3.1, 2.95), xytext=(5.9, 2.95), arrowprops=dict(arrowstyle="-|>", lw=0.9))
    ax.text(3.15, 2.63, "$-2v_A/Z-h_A$", fontsize=7, va="bottom")
    ax.text(5.85, 2.98, "$-2v_B/Z-h_B$", fontsize=7, va="bottom", ha="right")
    ax.text(4.5, -0.42, "each direction delayed by τ + kΔt", ha="center", va="top", fontsize=7)
    ax.set_xlim(-0.9, 9.3)
    ax.set_ylim(-0.75, 3.45)
    ax.set_aspect("equal")
    ax.axis("off")
    save_figure(fig, "phase3_circuit")


def main():
    circuit_figure()
    rows, steady_rows, runs = [], [], {}
    checks = {}
    for f in FREQS:
        t, vb_ref, ib_ref = simulate("reference", 0, f)
        runs[(f, "reference", 0)] = (vb_ref, ib_ref)
        pre, post = t < T_ON, t >= T_ON
        ref_fit = fit_phasor(t[pre], vb_ref[pre], f)
        for variant in VARIANTS:
            for k in DELAYS:
                _, vb, ib = simulate(variant, k, f)
                runs[(f, variant, k)] = (vb, ib)
                if variant == "split" and k == 0:
                    checks[f"split_vs_unsplit_max_difference_V_{f:g}Hz"] = float(np.max(abs(vb - vb_ref)))
                vb0, ib0 = runs[(f, variant, 0)]
                for window, m in [("full", slice(None)), ("pre", pre), ("post", post)]:
                    total = rel_rms(np.c_[vb, ib][m], np.c_[vb_ref, ib_ref][m])
                    delay = rel_rms(np.c_[vb, ib][m], np.c_[vb0, ib0][m])
                    model = rel_rms(np.c_[vb0, ib0][m], np.c_[vb_ref, ib_ref][m])
                    rows.append([f, variant, k, k * DT * 1e6, window, *total, *delay, *model,
                                 float(np.max(abs(vb - vb_ref)[m])), float(np.max(abs(ib - ib_ref)[m]))])
                fitted = fit_phasor(t[pre], vb[pre], f)
                predicted = steady_state(variant, k, f, DT)["vb"]
                steady_rows.append([f, variant, k, k * DT * 1e6, abs(fitted), np.angle(fitted, deg=True),
                                    abs(fitted - predicted),
                                    100 * abs(fitted - ref_fit) / abs(ref_fit),
                                    steady_error(variant, f, [k * DT])[0]])
    save_csv("phase3_errors", ["f_Hz", "variant", "k", "delay_us", "window",
                               "v_total_pct", "i_total_pct", "v_delay_only_pct", "i_delay_only_pct",
                               "v_model_pct", "i_model_pct", "v_max_abs_V", "i_max_abs_A"], rows)
    save_csv("phase3_steady_state", ["f_Hz", "variant", "k", "delay_us", "vB_amplitude_V", "vB_angle_deg",
                                     "fit_minus_discrete_prediction_V", "steady_error_pct",
                                     "continuous_time_steady_error_pct"], steady_rows)
    checks["max_fit_minus_discrete_prediction_V"] = max(r[6] for r in steady_rows)
    eq = equivalents(60.0)
    checks["equivalents_60Hz"] = dict(static=eq["static"], phasor=dict(
        gain_abs=abs(eq["phasor"]["gain"]), gain_angle_deg=np.angle(eq["phasor"]["gain"], deg=True),
        r=eq["phasor"]["r"], c_uF=eq["phasor"]["c"] * 1e6))
    save_json("phase3_summary", checks)
    print(checks)

    # Waveforms at 60 Hz.
    f = 60.0
    t = np.arange(round(T_END / DT) + 1) * DT
    vb_ref = runs[(f, "reference", 0)][0]
    show = [("split", 20, "Full split, +100 µs"), ("static", 0, "Resistive Thevenin"),
            ("phasor", 0, "Phasor Thevenin")]
    fig, axes = plt.subplots(3, 1, figsize=(6.6, 5.6))
    axes[0].plot(t * 1e3, vb_ref, color="k", lw=1.6, label="Reference (unsplit)")
    for (v, k, label), c in zip(show, ["C0", "C1", "C2"]):
        axes[0].plot(t * 1e3, runs[(f, v, k)][0], color=c, lw=0.9, label=label)
    zoom = (t >= 0.0195) & (t <= 0.0235)
    axes[1].plot(t[zoom] * 1e3, vb_ref[zoom], color="k", lw=1.6)
    for (v, k, label), c in zip(show, ["C0", "C1", "C2"]):
        axes[1].plot(t[zoom] * 1e3, runs[(f, v, k)][0][zoom], color=c, lw=0.9)
        axes[2].plot(t * 1e3, runs[(f, v, k)][0] - vb_ref, color=c, lw=0.9)
    axes[0].legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.02), fontsize=6.5)
    axes[0].set_ylabel("$v_B$ (V)")
    axes[1].set_ylabel("$v_B$ (V)")
    axes[2].set_ylabel("$v_B$ error (V)")
    axes[1].set_xlabel("Time (ms), around the closing of the shunt")
    axes[2].set_xlabel("Time (ms)")
    for ax in (axes[0], axes[2]):
        ax.axvspan(T_ON * 1e3, T_OFF * 1e3, color="0.92", lw=0)
        ax.set_xlim(0, T_END * 1e3)
    axes[1].axvline(T_ON * 1e3, color="0.5", lw=0.6, ls="--")
    panel_labels(axes)
    fig.tight_layout()
    save_figure(fig, "phase3_waveforms")

    # Error against added delay at 60 Hz.
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.6))
    for variant, marker in zip(VARIANTS, "osd"):
        ss = [r[7] for r in steady_rows if r[0] == f and r[1] == variant]
        post = [r[5] for r in rows if r[0] == f and r[1] == variant and r[4] == "post"]
        axes[0].plot([k * 5 for k in DELAYS], ss, marker + "-", label=NAMES[variant])
        axes[1].plot([k * 5 for k in DELAYS], post, marker + "-", label=NAMES[variant])
    axes[0].set(xlabel="Added one-way delay (µs)", ylabel="Steady-state $v_B$ error (%)",
                title="Before the switching event")
    axes[1].set(xlabel="Added one-way delay (µs)", ylabel="RMS $v_B$ error (%)",
                title="From the switching event on")
    axes[0].legend(fontsize=6.5)
    panel_labels(axes)
    fig.tight_layout()
    save_figure(fig, "phase3_error_vs_delay")

    # Steady-state error against delay for three source frequencies.
    d = np.linspace(0, 400e-6, 161)
    fig, ax = plt.subplots(figsize=(4.4, 2.9))
    for f, c in zip(FREQS, ["C0", "C1", "C2"]):
        ax.plot(d * 1e6, steady_error("split", f, d), color=c, label=f"{f:g} Hz, full split")
        ax.plot(d * 1e6, steady_error("static", f, d), color=c, ls="--", lw=0.9,
                label=f"{f:g} Hz, resistive Thevenin")
        sim = [r[7] for r in steady_rows if r[0] == f and r[1] == "split"]
        ax.plot([k * 5 for k in DELAYS], sim, "o", color=c, ms=3.5, mfc="white")
    ax.plot([], [], "o", color="0.3", ms=3.5, mfc="white", label="simulation, full split")
    ax.set(xlabel="Added one-way delay (µs)", ylabel="Steady-state $v_B$ error (%)")
    ax.legend(fontsize=6, ncol=1, loc="upper left")
    save_figure(fig, "phase3_frequency")


if __name__ == "__main__":
    main()
