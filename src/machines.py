"""Synchronous machine models for the two-area system.

SixState: round-rotor model with transient and subtransient (damper) states,
           delta, dw, e'q, e'd, psi1d, psi2q; no saturation.
EquivalentCircuit: the machine model used in Kundur's Example 12.6 (field and
           one d-axis damper, two q-axis circuits; stator transients neglected),
           with fundamental parameters from the classical conversion and,
           optionally, the incremental saturation used in the book's
           small-signal analysis.
Both are coupled to the network through kundur_system.Case.network.
"""
import numpy as np
from scipy.linalg import eig

from kundur_system import W0, MBASE, mode_label

# Standard parameters of Example 12.6 (900 MVA base).
RA, XL, XD, XQ, XD1, XQ1, X2 = 0.0025, 0.2, 1.8, 1.7, 0.3, 0.55, 0.25
TD10, TQ10, TD20, TQ20 = 8.0, 0.4, 0.03, 0.05
ASAT, BSAT, PSIT1 = 0.015, 9.6, 0.9

GD1 = (X2 - XL) / (XD1 - XL)
GQ1 = (X2 - XL) / (XQ1 - XL)
GD2 = (XD1 - X2) / (XD1 - XL) ** 2
GQ2 = (XQ1 - X2) / (XQ1 - XL) ** 2


def jacobian(f, x0, h):
    return np.column_stack([(f(x0 + e) - f(x0 - e)) / (2 * h) for e in np.eye(len(x0)) * h])


def electromechanical(values, vectors):
    """Pick the three rotor modes (1 < Im < 12 rad/s), sorted by frequency."""
    ix = sorted(np.where((values.imag > 1) & (values.imag < 12))[0], key=lambda i: values[i].imag)
    if len(ix) != 3:
        raise RuntimeError(f"Expected three electromechanical modes, found {len(ix)}.")
    return values[ix], vectors[:, ix]


class SixState:
    def __init__(self, case):
        self.case = case
        V = case.V[:4]
        I = (case.Y @ case.V)[:4] / MBASE
        delta = np.angle(V + (RA + 1j * XQ) * I)
        vr, ir = V * np.exp(-1j * delta), I * np.exp(-1j * delta)
        iq, id_ = ir.real, -ir.imag
        psi2d = vr.real + RA * iq + X2 * id_
        psi2q = -vr.imag + RA * id_ - X2 * iq
        e1q = psi2d + (XD1 - X2) * id_
        e1d = psi2q - (XQ1 - X2) * iq
        e2d = e1q - (XD1 - XL) * id_
        e2q = e1d + (XQ1 - XL) * iq
        self.x = np.r_[delta, np.zeros(4), e1q, e1d, e2d, e2q]
        self.efd = e1q + (XD - XD1) * id_
        self.tm = (V * np.conj(I)).real + RA * abs(I) ** 2
        self.residual = float(np.max(abs(self.rhs(self.x))))
        if self.residual > 1e-8:
            raise RuntimeError("Inconsistent initial state.")

    def electrical(self, x, return_voltage=False):
        d, w, e1q, e1d, e2d, e2q = x.reshape(6, 4)
        psi2d = GD1 * e1q + GD2 * (XD1 - XL) * e2d
        psi2q = GQ1 * e1d + (1 - GQ1) * e2q
        out = self.case.network(np.r_[d, w, psi2d, psi2q], X2, X2, RA, return_voltage)
        return (psi2d, psi2q) + tuple(out)

    def rhs(self, x, efd=None):
        efd = self.efd if efd is None else efd
        d, w, e1q, e1d, e2d, e2q = x.reshape(6, 4)
        _, _, iq, id_, pe = self.electrical(x)
        field = e1q + (XD - XD1) * (GD1 * id_ + GD2 * (e1q - e2d))
        quad = e1d + (XQ - XQ1) * (-GQ1 * iq + GQ2 * (e1d - e2q))
        torque = pe + RA * (iq ** 2 + id_ ** 2)
        return np.r_[W0 * w, (self.tm - torque) / (2 * self.case.H),
                     (efd - field) / TD10, -quad / TQ10,
                     (-id_ * (XD1 - XL) + e1q - e2d) / TD20,
                     (iq * (XQ1 - XL) + e1d - e2q) / TQ20]

    def terminal_voltage(self, x):
        return abs(self.electrical(x, True)[-1][:4])

    def eigen(self, h=1e-5):
        return eig(jacobian(self.rhs, self.x, h))

    def modes(self, h=1e-5):
        return electromechanical(*self.eigen(h))

    def exciter_modes(self, KA=200.0, TR=0.01, h=1e-6):
        """Rotor modes with a static exciter on every unit,
        Efd = KA (Vref - Vm), dVm/dt = (Vt - Vm)/TR, no transient gain reduction."""
        vt0 = self.terminal_voltage(self.x)
        vref = vt0 + self.efd / KA
        x0 = np.r_[self.x, vt0]

        def f(x):
            xs, vm = x[:24], x[24:]
            return np.r_[self.rhs(xs, KA * (vref - vm)), (self.terminal_voltage(xs) - vm) / TR]

        if np.max(abs(f(x0))) > 1e-7:
            raise RuntimeError("Inconsistent exciter equilibrium.")
        values, left, right = eig(jacobian(f, x0, h), left=True)
        participation = abs(left.conj() * right)
        rotor = participation[:8].sum(0) / participation.sum(0)
        ix = sorted([i for i in range(28) if 0.5 < values[i].imag < 12 and rotor[i] > 0.5],
                    key=lambda i: values[i].imag)
        if len(ix) != 3:
            raise RuntimeError(f"Expected three rotor modes, found {len(ix)}.")
        return values[ix], right[4:8][:, ix]


def fundamental_parameters():
    """Classical conversion of the standard parameters (Kundur, Sec. 4.2)."""
    lad, laq = XD - XL, XQ - XL
    lfd = lad * (XD1 - XL) / (lad - (XD1 - XL))
    l1d = 1 / (1 / (X2 - XL) - 1 / lad - 1 / lfd)
    l1q = laq * (XQ1 - XL) / (laq - (XQ1 - XL))
    l2q = 1 / (1 / (X2 - XL) - 1 / laq - 1 / l1q)
    return dict(Lad=lad, Laq=laq, Lfd=lfd, L1d=l1d, L1q=l1q, L2q=l2q,
                Rfd=(lad + lfd) / (W0 * TD10),
                R1d=(l1d + lad * lfd / (lad + lfd)) / (W0 * TD20),
                R1q=(laq + l1q) / (W0 * TQ10),
                R2q=(l2q + laq * l1q / (laq + l1q)) / (W0 * TQ20))


class EquivalentCircuit:
    """Example 12.6 machine model; saturation 'none' or 'incremental'.

    With 'incremental', the operating point uses the total saturation factor
    Ksd = psi_at/(psi_at + psi_I) and the linear model the incremental factor
    1/(1 + Bsat psi_I) on both axes, as in the book's small-signal analysis.
    """

    def __init__(self, case, saturation="none"):
        self.case, self.p = case, fundamental_parameters()
        p = self.p
        V = case.V[:4]
        I = (case.Y @ case.V)[:4] / MBASE
        if saturation == "incremental":
            psi_at = abs(V + (RA + 1j * XL) * I)
            psi_i = np.where(psi_at > PSIT1, ASAT * np.exp(BSAT * (psi_at - PSIT1)), 0.0)
            ksd, kinc = psi_at / (psi_at + psi_i), 1 / (1 + BSAT * psi_i)
        else:
            ksd = kinc = np.ones(4)
        lads, laqs = ksd * p["Lad"], ksd * p["Laq"]
        delta = np.angle(V + (RA + 1j * (XL + laqs)) * I)
        r = np.exp(-1j * delta)
        iq, id_ = (I * r).real, -(I * r).imag
        vq, vd = (V * r).real, -(V * r).imag
        psid, psiq = vq + RA * iq, -(vd + RA * id_)
        self.psiad0, self.psiaq0 = psid + XL * id_, psiq + XL * iq
        self.iad0, self.iaq0 = self.psiad0 / lads, self.psiaq0 / laqs
        ifd = self.iad0 + id_
        self.lkd, self.lkq = kinc * p["Lad"], kinc * p["Laq"]
        self.efd = p["Lad"] * ifd
        self.tm = psid * iq - psiq * id_
        self.x = np.r_[delta, np.zeros(4), p["Lfd"] * ifd + self.psiad0,
                       self.psiad0, self.psiaq0, self.psiaq0]
        if np.max(abs(self.rhs(self.x))) > 1e-8:
            raise RuntimeError("Inconsistent initial state.")

    def electrical(self, x):
        p = self.p
        d, w, pfd, p1d, p1q, p2q = x.reshape(6, 4)
        ldk = 1 / (1 / self.lkd + 1 / p["Lfd"] + 1 / p["L1d"])
        lqk = 1 / (1 / self.lkq + 1 / p["L1q"] + 1 / p["L2q"])
        psi2d = ldk * (self.psiad0 / self.lkd - self.iad0 + pfd / p["Lfd"] + p1d / p["L1d"])
        psi2q = lqk * (self.psiaq0 / self.lkq - self.iaq0 + p1q / p["L1q"] + p2q / p["L2q"])
        iq, id_, _, v = self.case.network(np.r_[d, w, psi2d, -psi2q], XL + ldk, XL + lqk, RA, True)
        psiad, psiaq = psi2d - ldk * id_, psi2q - lqk * iq
        currents = ((pfd - psiad) / p["Lfd"], (p1d - psiad) / p["L1d"],
                    (p1q - psiaq) / p["L1q"], (p2q - psiaq) / p["L2q"])
        torque = (psiad - XL * id_) * iq - (psiaq - XL * iq) * id_
        return currents, torque, abs(v[:4])

    def rhs(self, x, efd=None):
        p = self.p
        efd = self.efd if efd is None else efd
        (ifd, i1d, i1q, i2q), torque, _ = self.electrical(x)
        return np.r_[W0 * x[4:8], (self.tm - torque) / (2 * self.case.H),
                     W0 * p["Rfd"] * (efd / p["Lad"] - ifd),
                     -W0 * p["R1d"] * i1d, -W0 * p["R1q"] * i1q, -W0 * p["R2q"] * i2q]

    def with_excitation(self, kind="manual"):
        """State function and equilibrium for the excitation cases of Table E12.4:
        manual, thyristor (KA 200, TR 0.01), tgr (TA 1, TB 10) or pss
        (KSTAB 20, TW 10, T1 0.05, T2 0.02, T3 3, T4 5.4)."""
        if kind == "manual":
            return self.rhs, self.x
        KA, TR, TA, TB = 200.0, 0.01, 1.0, 10.0
        KS, TW, T1, T2, T3, T4 = 20.0, 10.0, 0.05, 0.02, 3.0, 5.4
        vt0 = self.electrical(self.x)[2]
        vref = vt0 + self.efd / KA
        tgr, pss = kind == "tgr", kind == "pss"
        x0 = np.r_[self.x, vt0, self.efd / KA if tgr else [], np.zeros(12 if pss else 0)]

        def f(x):
            xs, vm = x[:24], x[24:28]
            u, extra = vref - vm, []
            if tgr:
                z = x[28:32]
                efd = KA * (z + TA / TB * (u - z))
                extra.append((u - z) / TB)
            elif pss:
                xw, x1, x2 = x[28:32], x[32:36], x[36:40]
                y2 = KS * xs[4:8] - xw
                y3 = x1 + T1 / T2 * (y2 - x1)
                y4 = x2 + T3 / T4 * (y3 - x2)
                efd = KA * (u + y4)
                extra += [y2 / TW, (y2 - x1) / T2, (y3 - x2) / T4]
            else:
                efd = KA * u
            return np.r_[self.rhs(xs, efd), (self.electrical(xs)[2] - vm) / TR, *extra]

        return f, x0

    def rotor_modes(self, kind="manual", h=1e-6):
        f, x0 = self.with_excitation(kind)
        if np.max(abs(f(x0))) > 1e-8:
            raise RuntimeError("Inconsistent equilibrium.")
        values, left, right = eig(jacobian(f, x0, h), left=True)
        participation = abs(left.conj() * right)
        rotor = participation[:8].sum(0) / participation.sum(0)
        ix = sorted([i for i in range(len(x0)) if 1 < values[i].imag < 12 and rotor[i] > 0.3],
                    key=lambda i: values[i].imag)
        if len(ix) != 3:
            raise RuntimeError(f"Expected three rotor modes, found {len(ix)}.")
        return {mode_label(right[4:8, i]): values[i] for i in ix}

    def spectrum(self, kind="manual", h=1e-6):
        f, x0 = self.with_excitation(kind)
        return eig(jacobian(f, x0, h), right=False)
