# Rigorous parent-UV truncation certificate for the massive chemical Bose weak source (2026-11-23)

## Question and scope

The existing exact massive Yukawa 1↔2 weak source at nonzero chemical potentials [reference](phi_bose_chemical_weak_reference.py) was compared at parent cutoffs 2,8,16 GeV, but a finite-cutoff comparison alone does not certify the *omitted* parent thermal tail. Here we obtain an **analytic positive bound** on that missing number rate for a restricted thermal forward reaction. This is a cutoff certificate for a *prescribed* thermal distribution, not a physical late-time nonthermal spectrum or a strong solver.

Assume \(T_1=T_2=T_\phi=T>0\), fermionic parent/daughter, real scalar \(m_\phi>0\), \(\mu_\phi<m_\phi\), and constant positive chemical affinity \(\Delta\mu=\mu_1-\mu_2-\mu_\phi>0\). Use the same on-shell Yukawa rest width \(\Gamma\), parent degeneracy \(g_1\), and boosted isotropic rest-frame angle as the reference. The parent radial number source at momentum p is

\[
\frac{dR_\phi}{dp}=\frac{g_1\Gamma m_1}{2\pi^2}\frac{p^2}{E_1}
\left\langle f_1(E_1)[1-f_2(E_2)][1+b_\phi(k)]
(1-e^{-\Delta\mu/T})\right\rangle_{\cos\theta^*}.
\]

All factors are nonnegative. Since \(E_1=\sqrt{p^2+m_1^2}\ge p\), \(p^2/E_1\le p\),
\(f_1(E_1)\le e^{(\mu_1-E_1)/T}\le e^{(\mu_1-p)/T}\),
\(1-f_2\le1\), and on shell \(k\ge m_\phi\), the Bose factor obeys
\(1+b_\phi(k)\le 1+[\exp((m_\phi-\mu_\phi)/T)-1]^{-1}\).
The normalized rest-angle average of a constant is exactly one. Hence the **strict infinite-UV number-source remainder** satisfies

\[
\boxed{0\le R_\phi(\infty)-R_\phi(P)
\le \frac{g_1\Gamma m_1}{2\pi^2}
e^{\mu_1/T}\left(1+\frac{1}{\operatorname{expm1}((m_\phi-\mu_\phi)/T)}\right)
(1-e^{-\Delta\mu/T})\,T(P+T)e^{-P/T}.}
\]

This is valid without an angular grid or a numerical extrapolation; the elementary integral is \(\int_P^\infty p e^{-p/T}dp=T(P+T)e^{-P/T}\). It is *not* a certificate for a reverse-dominated net rate, unequal temperatures, nonthermal tails, energy-weighted moments, or a zero-mode source. The bound intentionally becomes useless as \(\mu_\phi\uparrow m_\phi\); a tighter finite-P kinematic bound was separately derived [here](../archive/stages/2026/11/2026-11-23-finite-parent-chemical-ir-bound.md).

## Numeric benchmark, independently checked by direct positive tail quadrature

Use \(m_1=.1\), \(m_2=.02\), \(m_\phi=86\times10^{-9}\), \(y=.001\), \(g_1=4\), \(T=1\), \(\mu_1=-.09\), \(\mu_2=.05\), \(\mu_\phi=-.2\) (all energies GeV), \(\Delta\mu=.06\). The reference integral at parent P=60 is \(R_\phi(60)\approx 3.474671097336273\times10^{-12}\ {\rm GeV}^4\); its own infinite remainder is bounded by the same closed expression evaluated at P=60. The numerical tail column is obtained from a *direct positive integral over [P,60]*, not by subtracting two nearly equal accumulated rates.

| P (GeV) | \(R_\phi(P)\) (GeV⁴) | direct omitted [P,60] (GeV⁴) | strict omitted [P,∞] upper bound (GeV⁴) | observed fraction of \(R_\phi(60)\) |
|---:|---:|---:|---:|---:|
| 2 | 1.9897147452494e−12 | 1.4849563520868e−12 | 6.6435952599147e−12 | 42.7366% |
| 8 | 3.4645263492995e−12 | 1.0144748036802e−14 | 4.9403478634212e−14 | 0.291963% |
| 12 | 3.4744110993118e−12 | 2.5999802443305e−16 | 1.3070146187389e−15 | 0.00748267% |
| 16 | 3.4746649796009e−12 | 6.1177353664750e−18 | 3.1304594788077e−17 | 0.000176067% |
| 20 | 3.4746709604937e−12 | 1.3684256674709e−19 | 7.0827274868399e−19 | 0.00000393829% |
| 24 | 3.4746710943765e−12 | 2.95979717619696e−21 | 1.5443414166223e−20 | 0.0000000851821% |

The bound is conservative (~4.5–5.2 times the directly measured tail over these cutoffs). The important distinction is that the **rigorous fractional bound relative to the truncated rate** at P=16 is \(3.13046\times10^{-17}/3.47466\times10^{-12}\approx9.01\times10^{-6}\), whereas the measured finite tail is \(\approx1.76\times10^{-6}\) of the near-complete rate. Thus the positive thermal **number** source is UV controlled at P=16 to better than 0.001% by a mathematical inequality, under exactly these trial parameters.

The same **P=2** source misses ~42.7% of the total, a substantial negative control. A fixed thermal cutoff benchmark does not justify freezing a nonthermal parent at P=16; in a cosmological evolution its tail and temperature must be checked afresh.

## Reproduction and next gate

[Source](phi_bose_chemical_uv_bound.py), [focused tests](../tests/test_phi_bose_chemical_uv_bound.py). Commands:

\`\`\`bash
OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_chemical_uv_bound
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_chemical_uv_bound.py
\`\`\`

The table was calculated with SciPy independent interval quadrature locally; public GitHub CI of the new regression is to be verified. The strict analytic bound is the new scientific result. The test does **not** certify φ daughter IR resolution, physical strong cell density of states, three energy moments, microscopic rates, evolving bath/FRW, 86-eV relic, Bose condensate, or CLASS. Next: obtain an analogous energy-weighted tail bound and converge parent radial/angle/φ grids for a *shared physical finite-volume strong collision* before any cosmological interpretation.
