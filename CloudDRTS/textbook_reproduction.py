"""Reproduce Kundur [1] Example 12.6 (Tables E12.3 and E12.4) eigenvalue by eigenvalue.

Machine: Kundur's equivalent circuit (model 2.2: field + one d damper, two q
circuits), stator transients and speed in the stator equations neglected,
equal mutual inductances, standard parameters converted with the classical
expressions of [1]. Without saturation this model has the same eigenvalues as
the GENROU model in detailed_machine.py (checked below).

Saturation as in the small-signal analysis of [1]: the operating point uses the
total saturation factor Ksd = psi_at/(psi_at + psi_I); the linear model uses the
incremental factor 1/(1 + Bsat psi_I) on both axes. Here this is a magnetizing
relation that is affine in each axis about the operating point, so its exact
Jacobian is that linear model.

Excitation (Table E12.4): (ii) thyristor KA 200, TR 0.01; (iii) with TGR
TA 1, TB 10; (iv) with PSS KSTAB 20, TW 10, T1 0.05, T2 0.02, T3 3, T4 5.4.
Case (i), the self-excited dc exciter, is not modelled.
"""
import json
from scipy.linalg import eig
from common import np, table, OUT
import kundur as kd
from literature_checks import Case
from detailed_machine import Detailed

W0 = 2*np.pi*60
LL, LD, LQ, LD1, LQ1, L2, RA = .2, 1.8, 1.7, .3, .55, .25, .0025
TD10, TD20, TQ10, TQ20 = 8., .03, .4, .05
ASAT, BSAT, PSIT1 = .015, 9.6, .9

def classical_parameters():
    """Fundamental parameters from standard parameters, classical expressions of [1]."""
    Lad, Laq = LD-LL, LQ-LL
    Lfd = Lad*(LD1-LL)/(Lad-(LD1-LL)); L1d = 1/(1/(L2-LL)-1/Lad-1/Lfd)
    L1q = Laq*(LQ1-LL)/(Laq-(LQ1-LL)); L2q = 1/(1/(L2-LL)-1/Laq-1/L1q)
    return dict(Lad=Lad, Laq=Laq, Lfd=Lfd, L1d=L1d, L1q=L1q, L2q=L2q,
                Rfd=(Lad+Lfd)/(W0*TD10), R1d=(L1d+Lad*Lfd/(Lad+Lfd))/(W0*TD20),
                R1q=(Laq+L1q)/(W0*TQ10), R2q=(L2q+Laq*L1q/(Laq+L1q))/(W0*TQ20))

def saturation_curve(psi_at):
    return np.where(psi_at > PSIT1, ASAT*np.exp(BSAT*(psi_at-PSIT1)), 0.)

class CircuitMachines:
    """Four units on the textbook network; saturation 'none' or 'incremental'."""
    def __init__(self, case, saturation='none'):
        if saturation not in ('none', 'incremental'):
            raise ValueError('saturation must be none or incremental')
        self.case, self.p = case, classical_parameters()
        p = self.p
        V = case.V[:4]; I = (case.Y@case.V+np.conj(case.load/case.V))[:4]/9
        if saturation == 'incremental':
            # Total saturation at the operating point (air-gap flux behind Ra + jXl).
            psi_at = abs(V+(RA+1j*LL)*I); psi_I = saturation_curve(psi_at)
            ksd, kinc = psi_at/(psi_at+psi_I), 1/(1+BSAT*psi_I)
        else:
            ksd = kinc = np.ones(4)
        Lads, Laqs = ksd*p['Lad'], ksd*p['Laq']
        delta = np.angle(V+(RA+1j*(LL+Laqs))*I); r = np.exp(-1j*delta)
        iq, id_ = (I*r).real, -(I*r).imag
        vq, vd = (V*r).real, -(V*r).imag
        psid, psiq = vq+RA*iq, -(vd+RA*id_)
        self.psiad0, self.psiaq0 = psid+LL*id_, psiq+LL*iq
        self.iad0, self.iaq0 = self.psiad0/Lads, self.psiaq0/Laqs
        ifd = self.iad0+id_
        self.Lkd, self.Lkq = kinc*p['Lad'], kinc*p['Laq']
        self.efd0 = p['Lad']*ifd           # Efd in the unsaturated Lad-reciprocal system
        self.tm = psid*iq-psiq*id_
        self.x = np.r_[delta, np.zeros(4), p['Lfd']*ifd+self.psiad0, self.psiad0, self.psiaq0, self.psiaq0]
        if max(abs(self.rhs(self.x))) > 1e-8:
            raise RuntimeError('Inconsistent operating point.')

    def electrical(self, x):
        p = self.p; d, w, pfd, p1d, p1q, p2q = x.reshape(6, 4)
        Ldk = 1/(1/self.Lkd+1/p['Lfd']+1/p['L1d']); Lqk = 1/(1/self.Lkq+1/p['L1q']+1/p['L2q'])
        psi2d = Ldk*(self.psiad0/self.Lkd-self.iad0+pfd/p['Lfd']+p1d/p['L1d'])
        psi2q = Lqk*(self.psiaq0/self.Lkq-self.iaq0+p1q/p['L1q']+p2q/p['L2q'])
        iq, id_, _, v = self.case.network(np.r_[d, w, psi2d, -psi2q], LL+Ldk, LL+Lqk, RA, True)
        psiad, psiaq = psi2d-Ldk*id_, psi2q-Lqk*iq
        currents = ((pfd-psiad)/p['Lfd'], (p1d-psiad)/p['L1d'], (p1q-psiaq)/p['L1q'], (p2q-psiaq)/p['L2q'])
        torque = (psiad-LL*id_)*iq-(psiaq-LL*iq)*id_
        return currents, torque, abs(v[:4])

    def rhs(self, x, efd=None):
        p = self.p; efd = self.efd0 if efd is None else efd
        (ifd, i1d, i1q, i2q), torque, _ = self.electrical(x)
        return np.r_[W0*x[4:8], (self.tm-torque)/(2*self.case.H), W0*p['Rfd']*(efd/p['Lad']-ifd),
                     -W0*p['R1d']*i1d, -W0*p['R1q']*i1q, -W0*p['R2q']*i2q]

def system(m, excitation='manual'):
    """State function and equilibrium for manual, thyristor, tgr or pss excitation."""
    if excitation == 'manual':
        return m.rhs, m.x
    KA, TR, TA, TB = 200., .01, 1., 10.
    KS, TW, T1, T2, T3, T4 = 20., 10., .05, .02, 3., 5.4
    vt0 = m.electrical(m.x)[2]; vref = vt0+m.efd0/KA
    tgr, pss = excitation == 'tgr', excitation == 'pss'
    x0 = np.r_[m.x, vt0, m.efd0/KA if tgr else [], np.zeros(12 if pss else 0)]
    def f(x):
        xs, vm = x[:24], x[24:28]; k = 28; u = vref-vm; extra = []
        if tgr:
            z = x[k:k+4]; efd = KA*(z+TA/TB*(u-z)); extra.append((u-z)/TB)
        elif pss:
            xw, x1, x2 = x[k:k+4], x[k+4:k+8], x[k+8:k+12]
            y2 = KS*xs[4:8]-xw; y3 = x1+T1/T2*(y2-x1); y4 = x2+T3/T4*(y3-x2)
            efd = KA*(u+y4); extra += [y2/TW, (y2-x1)/T2, (y3-x2)/T4]
        else:
            efd = KA*u
        return np.r_[m.rhs(xs, efd), (m.electrical(xs)[2]-vm)/TR, *extra]
    return f, x0

def eigen(f, x0, h=1e-6):
    A = np.column_stack([(f(x0+e)-f(x0-e))/(2*h) for e in np.eye(len(x0))*h])
    return eig(A, left=True)

def rotor_modes(m, excitation):
    f, x0 = system(m, excitation)
    if max(abs(f(x0))) > 1e-8:
        raise RuntimeError('Inconsistent equilibrium.')
    values, left, right = eigen(f, x0)
    participation = abs(left.conj()*right)
    mechanical = participation[:8].sum(0)/participation.sum(0)
    ix = sorted([i for i in range(len(x0)) if 1 < values[i].imag < 12 and mechanical[i] > .3],
                key=lambda i: values[i].imag)
    if len(ix) != 3:
        raise RuntimeError(f'Expected three rotor modes, found {len(ix)}.')
    labels = [kd.mode_label(right[4:8, i]) for i in ix]
    return dict(zip(labels, values[ix]))

# [1] Table E12.3 (manual excitation): one entry per real eigenvalue or complex pair;
# the near-zero pair 1,2 (-0.76e-3 +- j0.22e-2) is compared separately.
TABLE_E123 = [-.096, -.111+3.43j, -.117, -.265, -.276, -.492+6.82j, -.506+7.02j, -3.428, -4.139,
              -5.287, -5.303, -31.03, -32.45, -34.07, -35.53, -37.89+.142j, -38.01+.038j]
# [1] Table E12.4: rotor modes (inter-area, area 1 local, area 2 local).
TABLE_E124 = {'thyristor': (.031+3.84j, -.490+7.15j, -.496+7.35j),
              'tgr': (.123+3.46j, -.450+6.86j, -.462+7.05j),
              'pss': (-.501+3.77j, -1.826+8.05j, -1.895+8.35j)}

def main():
    case = Case(kd.Y, kd.load, kd.V0, kd.H, 1)
    rows = []
    # 1) Equivalence with GENROU (detailed_machine.py) without saturation.
    circuit = np.sort_complex(eigen(*system(CircuitMachines(case, 'none')))[0])
    genrou = np.sort_complex(Detailed(case, 'none').spectrum())
    equivalence = float(max(abs(circuit-genrou)))
    # 2) Table E12.3, full spectrum.
    m = CircuitMachines(case, 'incremental')
    spectrum = eigen(*system(m))[0]
    near_zero = sorted(spectrum, key=abs)[:2]
    rest = [z for z in spectrum if abs(z) > 1e-3 and z.imag >= 0]
    e123 = []
    for ref in TABLE_E123:
        ours = min(rest, key=lambda z: abs(z-ref)); rest.remove(ours)
        e123.append(abs(ours-ref))
        rows.append(['E12.3', 'manual', ref.real, ref.imag, ours.real, ours.imag, abs(ours-ref)])
    # 3) Table E12.4 rotor modes.
    e124 = {}
    for excitation, refs in TABLE_E124.items():
        modes = rotor_modes(m, excitation)
        diffs = []
        for label, ref in zip(['inter-area', 'area-1', 'area-2'], refs):
            ours = modes[label]; diffs.append(abs(ours-ref))
            rows.append(['E12.4', f'{excitation}/{label}', ref.real, ref.imag, ours.real, ours.imag, abs(ours-ref)])
        e124[excitation] = max(diffs)
    table('textbook_reproduction', ['table', 'case', 'ref1_real', 'ref1_imag', 'real', 'imag', 'abs_difference'], rows)
    manual = rotor_modes(m, 'manual')
    summary = dict(
        genrou_vs_circuit_no_saturation_max_eigenvalue_difference=equivalence,
        table_E123_max_abs_difference=max(e123), table_E123_entries=len(e123),
        table_E123_near_zero_pair=[[z.real, z.imag] for z in near_zero],
        table_E124_max_abs_difference=e124,
        manual_modes={k: dict(sigma=v.real, frequency_Hz=v.imag/(2*np.pi), damping_pct=-100*v.real/abs(v))
                      for k, v in manual.items()},
        not_modelled='Table E12.4 case (i), self-excited dc exciter.')
    (OUT/'textbook_reproduction.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
