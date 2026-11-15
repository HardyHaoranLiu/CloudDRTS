"""Phase 1: modal analysis and three-phase fault of the Kundur two-area system.

Six-state machines, manual excitation (constant Efd), constant mechanical
torque, no saturation. Active load constant current, reactive load constant
impedance. The state matrix is formed by central differences at the
equilibrium; the three electromechanical modes are labelled from their speed
mode shapes (kundur_system.mode_label).
"""
import numpy as np
from scipy.integrate import solve_ivp

from common import plt, save_csv, save_json, save_figure, panel_labels
from kundur_system import textbook_case, mode_label, tie_power, H_BOOK
from machines import SixState

FAULT_BUS, FAULT_ON, FAULT_OFF, T_END = 8, 1.0, 1.08, 12.0
R_FAULT = 1e-4


def modal_check(model, values, vectors, amp=1e-4, t_end=5.0):
    """Start the nonlinear model on the real part of each mode shape and compare
    with the linear prediction Re(v exp(lambda t)). Returns the largest state
    error relative to the initial displacement, per mode."""
    errors = []
    t = np.linspace(0, t_end, 501)
    for value, vec in zip(values, vectors.T):
        vec = vec / np.max(abs(vec))
        sol = solve_ivp(lambda _, x: model.rhs(x), (0, t_end), model.x + amp * vec.real,
                        method="DOP853", t_eval=t, rtol=1e-10, atol=1e-12, max_step=0.025)
        if not sol.success:
            raise RuntimeError(sol.message)
        linear = model.x[:, None] + amp * np.real(vec[:, None] * np.exp(value * t))
        errors.append(float(np.max(abs(sol.y - linear)) / amp))
    return errors


def fault_run(model):
    case = model.case
    Y0 = case.Y.copy()
    log_t, log_x, log_v = [], [], []
    x = model.x.copy()
    for a, b, faulted in [(0, FAULT_ON, False), (FAULT_ON, FAULT_OFF, True), (FAULT_OFF, T_END, False)]:
        case.Y = Y0.copy()
        if faulted:
            case.Y[FAULT_BUS - 1, FAULT_BUS - 1] += 1 / R_FAULT
        t = np.linspace(a, b, round((b - a) / 0.005) + 1)
        sol = solve_ivp(lambda _, z: model.rhs(z), (a, b), x, method="DOP853",
                        t_eval=t, rtol=1e-8, atol=1e-10, max_step=0.025)
        if not sol.success:
            raise RuntimeError(sol.message)
        x = sol.y[:, -1]
        keep = len(t) if b == T_END else len(t) - 1   # the next segment starts at b
        for ti, xi in zip(t[:keep], sol.y.T[:keep]):
            log_t.append(ti)
            log_x.append(xi)
            log_v.append(model.electrical(xi, True)[-1])
    case.Y = Y0
    return np.array(log_t), np.array(log_x), np.array(log_v)


def mode_shape_figure(rows):
    fig, axes = plt.subplots(1, 3, figsize=(3.45, 1.6), subplot_kw=dict(projection="polar"))
    titles = {"inter-area": "Inter-area", "area-1": "Local, area 1", "area-2": "Local, area 2"}
    for ax, r in zip(axes, rows):
        for g in range(4):
            ax.annotate("", xy=(np.deg2rad(r[f"G{g+1}_phase"]), r[f"G{g+1}_amp"]), xytext=(0, 0),
                        arrowprops=dict(arrowstyle="-|>", lw=1.2, color=f"C{g}",
                                        shrinkA=0, shrinkB=0, mutation_scale=8))
        ax.set_ylim(0, 1.05)
        ax.set_yticks([0.5, 1])
        ax.set_yticklabels([])
        ax.set_xticks(np.deg2rad([0, 90, 180, 270]))
        ax.set_xticklabels(["0°", "", "180°", ""], fontsize=6.5)
        ax.tick_params(pad=-1)
        ax.grid(lw=0.3)
        ax.set_title(f"{titles[r['mode']]}\n{r['frequency_Hz']:.3f} Hz", pad=6, fontsize=7.5)
    handles = [plt.Line2D([], [], color=f"C{g}", lw=1.5) for g in range(4)]
    fig.legend(handles, ["G1", "G2", "G3", "G4"], ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.06),
               handlelength=1.4, columnspacing=1.2)
    fig.subplots_adjust(wspace=0.45, bottom=0.12, top=0.72, left=0.04, right=0.96)
    save_figure(fig, "phase1_mode_shapes")


def fault_figure(t, angle, tie, coi_df):
    fig, axes = plt.subplots(3, 1, figsize=(3.45, 2.5), sharex=True)
    for g in range(4):
        axes[0].plot(t, angle[:, g], label=f"G{g+1}")
    axes[0].legend(ncol=4, loc="lower right", bbox_to_anchor=(1.0, 0.98), handlelength=1.2,
                   columnspacing=0.8, fontsize=6, borderaxespad=0.1)
    axes[1].plot(t, tie, color="k")
    axes[2].plot(t, coi_df, color="k")
    labels = ["$\\delta - \\delta_\\mathrm{COI}$ (deg)", "$P_{7\\to 8}$ (MW)", "$\\Delta f_\\mathrm{COI}$ (Hz)"]
    for ax, label in zip(axes, labels):
        ax.axvspan(FAULT_ON, FAULT_OFF, color="0.85", lw=0)
        ax.set_ylabel(label)
        ax.yaxis.set_label_coords(-0.11, 0.5)
    axes[2].set_xlabel("Time (s)")
    axes[2].set_xlim(0, T_END)
    panel_labels(axes)
    fig.align_ylabels(axes)
    fig.subplots_adjust(hspace=0.32)
    save_figure(fig, "phase1_fault")


def main():
    case = textbook_case()
    model = SixState(case)
    values, vectors = model.modes(1e-5)
    jac_change = float(np.max(abs(values - model.modes(5e-6)[0])))
    spectrum = model.eigen(1e-5)[0]

    rows = []
    for value, vec in zip(values, vectors.T):
        shape = vec[4:8] / vec[4:8][np.argmax(abs(vec[4:8]))]
        r = dict(mode=mode_label(shape), sigma=value.real, omega=value.imag,
                 frequency_Hz=value.imag / (2 * np.pi), damping_pct=-100 * value.real / abs(value))
        for g in range(4):
            r[f"G{g+1}_amp"] = abs(shape[g])
            r[f"G{g+1}_phase"] = np.angle(shape[g], deg=True)
        rows.append(r)
    if sorted(r["mode"] for r in rows) != ["area-1", "area-2", "inter-area"]:
        raise RuntimeError("Unexpected mode labels.")
    save_csv("phase1_modes", list(rows[0]), [list(r.values()) for r in rows])
    mode_shape_figure(rows)

    modal_errors = modal_check(model, values, vectors)

    t, x, v = fault_run(model)
    angle = np.rad2deg(x[:, :4] - model.x[:4])
    angle -= np.average(angle, axis=1, weights=H_BOOK)[:, None]
    tie = np.array([tie_power(vi) for vi in v])
    coi_df = np.average(x[:, 4:8] * 60, axis=1, weights=H_BOOK)
    fault_figure(t, angle, tie, coi_df)

    save_json("phase1_summary", dict(
        power_flow_residual_pu=case.pf_residual,
        tie_flow_7_to_8_MW=float(tie[0]),
        equilibrium_residual=model.residual,
        jacobian_step_halving_max_change=jac_change,
        positive_real_eigenvalues=[float(z.real) for z in spectrum if z.real > 1e-3],
        modal_check_max_error_ratio=max(modal_errors),
        prefault_max_speed_deviation_Hz=float(np.max(abs(x[t < FAULT_ON, 4:8])) * 60),
        min_bus8_voltage_during_fault_pu=float(np.min(abs(v[(t >= FAULT_ON) & (t < FAULT_OFF), 7]))),
        final_COI_frequency_deviation_Hz=float(coi_df[-1])))
    for r in rows:
        print(f"{r['mode']:>10}: {r['frequency_Hz']:.4f} Hz, {r['damping_pct']:.2f} %")


if __name__ == "__main__":
    main()
