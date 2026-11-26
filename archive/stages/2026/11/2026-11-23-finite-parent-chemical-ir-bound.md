# Finite-parent-cutoff chemical Bose-IR bound — positive on-shell weak negative control (2026-11-23)

## Question

Does a nonzero-chemical-potential **massive** thermal \(\chi_1\leftrightarrow\chi_2+\phi\) weak collision reference at \(\mu_\phi\uparrow m_\phi\) automatically probe the physical infrared \(p_\phi\to0\), Bose zero mode, or an unbounded Bose stimulated rate? **No**, if the parent momentum is finitely truncated below the extreme boost needed to reach that region.

The result builds on [exact parent support](../../../../kinetics/PHI-BOSE-PARENT-CUTOFF-KINEMATIC-IR-AUDIT.md) and the [new independent massive chemical weak reference](../../../../kinetics/PHI-BOSE-CHEMICAL-WEAK-REFERENCE-AUDIT.md). Here the quantum kernel is **positive forward**, with \(T_1=T_2=T_\phi=1\) GeV, \(\mu_2=0\), and fixed \(\Delta\mu=\mu_1-\mu_2-\mu_\phi=+0.06\) GeV; as \(\mu_\phi\) rises from 0 to \(m_\phi-\epsilon\), we shift \(\mu_1=\mu_\phi+0.06\) to keep the same microscopic chemical affinity. Masses \(m_1=.1\), \(m_2=.02\), \(m_\phi=86\times10^{-9}\) GeV, coupling \(y=.001\), g1=4 are **illustrative** Yukawa choices, NOT cosmological parameter fits.

## Rigorous finite-support bound

For any fixed physical decay shell energy \(E=E_2+k\), the positive net source is

\[
J(\mu_\phi) =
(1-f_1(E;\mu_1))\,f_2(E-k;0)\,b_\phi(k;\mu_\phi)
\,[\exp(\Delta\mu/T)-1].
\]

At fixed \(\Delta\mu>0\) and \(T\), shifting \(\mu_\phi\) and \(\mu_1\) together by \(0<\delta\le m_\phi\) preserves all factors except \(1-f_1\) and \(b_\phi\). Since
\[
\frac{\partial}{\partial \mu_1}\ln(1-f_1)=-f_1/T\in[-1/T,0],
\]
the parent unblocking ratio is between \(e^{-\delta/T}\) and 1. For the scalar,
\[
\frac{b_\phi(k;\delta)}{b_\phi(k;0)}
=\frac{\exp(k/T)-1}{\exp((k-\delta)/T)-1}
\]
is decreasing in \(k\) when \(\delta>0\). If every event has \(k\ge k_{\min}>\delta\), the **positive-integrand** collision source \(R_\phi>0\) therefore obeys
\[
\boxed{
e^{-\delta/T}
\;\le\;
\frac{R_\phi(\delta)}{R_\phi(0)}
\;\le\;
\frac{\operatorname{expm1}(k_{\min}/T)}
     {\operatorname{expm1}((k_{\min}-\delta)/T)}
}.
\]

No Jacobian, decay width, spin multiplicity or parent distribution normalization is needed for this inequality beyond positivity and the stated fixed-thermal-\(\Delta\mu\) assumptions. It is **not** an upper bound on a full untruncated nonthermal FRW history, other 2→2/number-changing channels, or a separately evolving \(p_\phi=0\) condensate. With a hard \(p_{\rm parent}\le p_{\max}\) cutoff at 2,8,16 GeV, the exact boosted kinematic minimum is attained at \(p_{\max}\) within this tested range.

## Independent numerical confirmation

Using [analytic 1D massive chemical source](../../../../kinetics/phi_bose_chemical_weak_reference.py) at \(\delta=m_\phi-10^{-12}\) GeV, \(T=1\) GeV and the independent exact boosted kinematic lower energy:

| parent pmax (GeV) | exact kmin (GeV) | \(R_\phi(\delta)/R_\phi(0)\) | rigorous upper bound | actual rise |
|---:|---:|---:|---:|---:|
| 2 | 0.00119925093757947 | 1.000005953704 | 1.000071758750 | 5.954 ppm |
| 8 | 0.000299988288329003 | 1.000008057002 | 1.000286799745 | 8.057 ppm |
| 16 | 0.000149998547511642 | 1.000008078143 | 1.000573704142 | 8.078 ppm |

The common lower bound is \(e^{-\delta/T}=0.999999914001\). For pmax16, the physical scalar rest energy \(m_\phi=8.6\times10^{-8}\) GeV is roughly **1744 times below** the smallest accessible daughter energy; the *entire* sampled scalar IR window is excluded by the parent cutoff. A near-critical chemical potential therefore changes this test's integrated rate only by \(\mathcal O(10^{-5})\), with a rigorous maximum rise of ~0.0574% for the tested pmax16; this is not evidence of saturation or absence of a condensate in the actual Universe.

**Reproduction (no new solver, no synthetic CI claim):**

\`\`\`python
from kinetics.phi_bose_chemical_weak_reference import MPH, analytic, kinematics
from math import exp, expm1
for pmax in (2., 8., 16.):
    delta = MPH - 1e-12
    kmin = kinematics(pmax)[1]
    r0 = analytic(pmax, 1., 1., .06, 0., 0.)["Rphi"]
    r1 = analytic(pmax, 1., 1., .06+delta, 0., delta)["Rphi"]
    print(pmax, kmin, r1/r0, expm1(kmin)/expm1(kmin-delta))
\`\`\`

The numerical scan is a local, deterministic analytic-quadrature calculation; it does not add a new independent physical collision solver. To assess genuinely critical 86-eV Bose production, extend parent UV with controlled thermal-tail errors and include independently derived low-momentum/zero-mode channels, repeated angular/qmin/dt convergence, SM rates and FRW energy accounting. The current positive weak gate does not establish physical \(f_\phi(q)\), BEC, relic \(\Omega_\phi\), or CLASS-ready predictions.
