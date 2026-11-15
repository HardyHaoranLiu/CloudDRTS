"""Kundur two-area, four-machine system: network, power flow and network solution.

Data of Kundur, Power System Stability and Control (1994), Example 12.6.
Buses are numbered 1-11 as in the book; array index = bus - 1.
Network base 100 MVA, machine base 900 MVA, 60 Hz.
"""
import numpy as np
from scipy.optimize import root

W0 = 2 * np.pi * 60.0
MBASE = 9.0                              # 900 MVA machine base on the 100 MVA system base
H_BOOK = np.array([6.5, 6.5, 6.175, 6.175])
Z_KM = 1e-4 + 1e-3j                      # line series impedance, p.u./km
B_KM = 1.75e-3                           # line charging, p.u./km
TRANSFORMERS = [(1, 5), (2, 6), (3, 11), (4, 10)]
LINES = [(5, 6, 25), (6, 7, 10), (7, 8, 110), (8, 9, 110), (9, 10, 10), (10, 11, 25)]


def admittance(transformer_x=0.15, tie_circuits=2, shunt=(2.0, 3.5)):
    """Bus admittance matrix.

    transformer_x: step-up transformer reactance on the 900 MVA base.
    tie_circuits: parallel circuits on the two 110 km sections 7-8 and 8-9.
    shunt: capacitor susceptance at buses 7 and 9 (p.u. on 100 MVA).
    """
    Y = np.zeros((11, 11), complex)

    def branch(a, b, y, charging=0j):
        a, b = a - 1, b - 1
        Y[a, a] += y + charging / 2
        Y[b, b] += y + charging / 2
        Y[a, b] -= y
        Y[b, a] -= y

    for a, b in TRANSFORMERS:
        branch(a, b, MBASE / (1j * transformer_x))
    for a, b, km in LINES:
        n = tie_circuits if km == 110 else 1
        branch(a, b, n / (km * Z_KM), 1j * n * km * B_KM)
    Y[6, 6] += 1j * shunt[0]
    Y[8, 8] += 1j * shunt[1]
    return Y


def textbook_power_flow(Y):
    """Example 12.6 operating point: G1 slack, G2-G4 at 700, 719 and 700 MW."""
    load = np.zeros(11, complex)
    load[[6, 8]] = [9.67 + 1j, 17.67 + 1j]
    p = -load.real.copy()
    p[1:4] = [7, 7.19, 7]

    def volts(x):
        return np.r_[1.03, 1.01, 1.03, 1.01, x[10:]] * np.exp(1j * np.r_[0, x[:10]])

    def mismatch(x):
        v = volts(x)
        s = v * np.conj(Y @ v)
        return np.r_[s.real[1:] - p[1:], (s.imag + load.imag)[4:]]

    sol = root(mismatch, np.r_[np.zeros(10), np.ones(7)], tol=1e-11)
    residual = float(np.max(abs(sol.fun)))
    if residual > 1e-8:
        raise RuntimeError("Power flow did not converge.")
    return volts(sol.x), load, residual


def symmetric_power_flow(Y):
    """Symmetric zero-transfer case of Klein et al.: 1400 MW generation and
    1367 MW load per area, G1 and G3 both at zero angle."""
    load = np.zeros(11, complex)
    load[[6, 8]] = 13.67 + 1j
    free = np.array([1, 3, 4, 5, 6, 7, 8, 9, 10])

    def volts(z):
        angle = np.zeros(11)
        angle[free] = z[:9]
        return np.r_[1.03, 1.01, 1.03, 1.01, z[9:]] * np.exp(1j * angle)

    def mismatch(z):
        v = volts(z)
        s = v * np.conj(Y @ v) + load
        p = s.real[free].copy()
        p[:2] -= 7
        return np.r_[p, s.imag[4:]]

    sol = root(mismatch, np.r_[np.zeros(9), np.ones(7)], tol=1e-11)
    residual = float(np.max(abs(sol.fun)))
    if residual > 1e-9:
        raise RuntimeError("Power flow did not converge.")
    return volts(sol.x), load, residual


class Case:
    """Network with voltage-dependent loads, solved for given machine EMFs.

    Active load varies as V**active_exponent (1: constant current,
    2: constant impedance); reactive load is constant impedance.
    """

    def __init__(self, Y, load, V, H, active_exponent):
        self.Y, self.load, self.V, self.H = Y, load, V, H
        self.exponent = active_exponent
        # Initial guess for the network solver (bus voltages, dq currents).
        current = (Y @ V)[:4] / MBASE
        rotated = current * np.exp(-1j * np.angle(V[:4] + 1.7j * current))
        self.guess = np.r_[V.real, V.imag, rotated.real, -rotated.imag]

    def network(self, x, xd, xq, ra, return_voltage=False):
        """Solve the network for rotor angles x[0:4] and the EMFs x[8:12]
        (q axis) and x[12:16] (d axis) behind xd, xq and ra."""
        d, _, eq, ed = x.reshape(4, 4)
        rot = np.exp(1j * d)

        def residual(z):
            v = z[:11] + 1j * z[11:22]
            iq, id_ = z[22:26], z[26:30]
            scale = abs(v) / abs(self.V)
            power = self.load.real * scale ** self.exponent + 1j * self.load.imag * scale ** 2
            kcl = self.Y @ v + np.conj(power / v)
            kcl[:4] -= MBASE * rot * (iq - 1j * id_)
            vdq = v[:4] / rot
            return np.r_[kcl.real, kcl.imag,
                         vdq.real + ra * iq + xd * id_ - eq,
                         -vdq.imag + ra * id_ - xq * iq - ed]

        sol = root(residual, self.guess, tol=1e-10)
        if np.max(abs(residual(sol.x))) > 1e-8:
            raise RuntimeError("Network equations did not converge.")
        z = sol.x
        v = z[:11] + 1j * z[11:22]
        iq, id_ = z[22:26], z[26:30]
        pe = (v[:4] * np.conj(rot * (iq - 1j * id_))).real
        return (iq, id_, pe, v) if return_voltage else (iq, id_, pe)


def textbook_case():
    Y = admittance()
    V, load, residual = textbook_power_flow(Y)
    case = Case(Y, load, V, H_BOOK, active_exponent=1)
    case.pf_residual = residual
    return case


def symmetric_case(transformer_x=0.15, tie_circuits=1, shunt=2.0):
    Y = admittance(transformer_x, tie_circuits, (shunt, shunt))
    V, load, residual = symmetric_power_flow(Y)
    case = Case(Y, load, V, np.full(4, 6.5), active_exponent=2)
    case.pf_residual = residual
    return case


def mode_label(shape):
    """Label a G1..G4 speed mode shape as 'inter-area', 'area-1', 'area-2'
    or, when both areas swing internally to a similar degree, 'local'."""
    s = np.asarray(shape, complex)
    s = s / s[np.argmax(abs(s))]
    a, b = abs(s[0] - s[1]), abs(s[2] - s[3])
    if abs(s[:2].mean() - s[2:].mean()) > max(a, b):
        return "inter-area"
    if min(a, b) > 0.5 * max(a, b):
        return "local"
    return "area-1" if a > b else "area-2"


def tie_power(v):
    """Active power (MW) from bus 7 towards bus 8 on the two 110 km circuits."""
    return 100 * (v[6] * np.conj((v[6] - v[7]) * 2 / (110 * Z_KM))).real
