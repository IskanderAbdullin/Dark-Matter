# Massive Yukawa Bose/FD weak reference at nonzero chemical potentials (2026-11-23)

## Physical question and scope

The preceding [real-dilog exact thermal weak source](PHI-BOSE-THERMAL-EXACT-WEAK-REFERENCE-AUDIT.md) treated a massive 86-eV scalar at zero chemical potentials and unequal parent/daughter temperatures. The [older thermal-affinity reference](THERMAL-CHEMICAL-AFFINITY.md) derived an analytic chemical source for a **massless** scalar at common temperature. Neither by itself supplies a chemical-potential-resolved, **massive** boosted on-shell number **and energy** reference with independent quadrature. We construct one for prescribed occupations, not for a strong evolution.

Benchmark (illustrative **GeV** units): Dirac parent mass m1=0.1, daughter m2=0.02, **real** scalar mφ=86e-9, Yukawa coupling y=0.001, parent spin/antiparticle multiplicity g1=4, daughter T2=Tφ=T0=1, μ2=+0.05, μφ=−0.2 (<mφ); parent T1=1 or 0.9 and μ1 varied. All species remain fixed thermal occupations. No thermal finite SM, FRW or physical chemical potential history is asserted.

## General closed daughter primitives

At each parent energy E and boosted scalar energy k∈[kmin,kmax], E2=E−k, define

\[
z=\exp[(k-\mu_\phi)/T_0],\qquad
a=\exp[(E-\mu_2-\mu_\phi)/T_0].
\]

Then the exact inverse quantum occupation product is

\[
 f_2(E-k)\,b_\phi(k)\,dk
 =\frac{T_0\,dz}{(z-1)(z+a)}.
\]

The integral of the positive inverse factor follows from
\[
 I_0(E)=\frac{T_0}{1+a}
  \left[\ln\frac{z-1}{z+a}\right]_{z_{\min}}^{z_{\max}}.
\]

The **physical scalar energy** integral is
\[
 I_1(E)=\mu_\phi I_0(E)+\frac{T_0^2}{1+a}
 [P(z_{\max})-P(z_{\min})],
\]
where
\[
 P(z)=-\operatorname{Li}_2(1-z)
      -\ln z\ln(1+z/a)-\operatorname{Li}_2(-z/a).
\]

For each exact on-shell event the affinity is
\[
 {\cal A}= -\frac{E-\mu_1}{T_1}
 +\frac{E-\mu_2-\mu_\phi}{T_0},\qquad
 F-R=(1-f_1)f_2 b_\phi\,\operatorname{expm1}({\cal A}).
\]

At T1=T0 the entire sign is fixed by Δμ=μ1−μ2−μφ. Common-temperature chemical equilibrium Δμ=0 has **zero net source on every microscopic shell**, irrespective of the shape of individual thermal spectra. The positive parent-angle measure is
\( [g_1p^2dp/(2\pi^2)]\,\Gamma_{\rm rest}m_1/E\,(d\cos\theta^*/2)\);
the exact daughter-energy Jacobian is \(d\cos\theta^*=dk/b\), with \(b=pp^*/m_1\). A separate 1D adaptive parent quadrature using I0,I1 is compared to a **different explicit 2D Gauss parent × log-k angle code**, which directly evaluates F/R and all three physical energies on each event. The scalar energy integral is independently checked against adaptive direct daughter integration.

## Verified local numerical reference

| parent pmax (GeV) | Δμ / temperatures | exact 1D Rφ (GeV⁴) | exact 1D Qφ (GeV⁵) | high-grid relative R error | high-grid relative Qφ error |
|---:|---|---:|---:|---:|---:|
| 2 | +0.06, T1=1 | +1.9897147452e−12 | +7.7311510048e−13 | +2.8e−14 | −2.9e−12 |
| 8 | +0.06, T1=1 | +3.4645263493e−12 | +2.4128996922e−12 | −1.5e−12 | −9.5e−13 |
| 16 | +0.06, T1=1 | +3.4746649796e−12 | +2.4488853230e−12 | −1.31e−9 | −1.13e−11 |
| 2 | −0.06, T1=1 | −1.9264518526e−12 | −7.4543230743e−13 | +2.9e−14 | −3.0e−12 |
| 8 | −0.06, T1=1 | −3.3230242600e−12 | −2.2963619899e−12 | −1.6e−12 | −1.0e−12 |
| 16 | −0.06, T1=1 | −3.3325726197e−12 | −2.3302524978e−12 | −1.36e−9 | −1.2e−11 |
| 2 | Δμ=0, T1=0.9 | −4.3665885859e−12 | −1.9455914217e−12 | +4.9e−15 | −6.1e−11 |
| 8 | Δμ=0, T1=0.9 | −1.1899779971e−11 | −1.1109348106e−11 | −1.3e−13 | −1.1e−11 |
| 16 | Δμ=0, T1=0.9 | −1.2004664881e−11 | −1.1484757414e−11 | −1.1e−10 | −1.0e−11 |

At exact common-temperature chemical equilibrium μ1=−0.15, T1=T0=1, the residual |Rφ| is below 8e−27 GeV⁴ at pmax=2/8/16 (finite floating-point affinity), with vanishing microscopic entropy to floating precision. On every forward/reverse case, independently accumulated Q1+Q2+Qφ is below ~7e−27 GeV⁵ at high grid 72×48; each forward/reverse entropy rate is nonnegative. Independent **direct** scalar energy daughter integrals agree with the real-dilog primitive for p=0.1/1/2/8/16 GeV under 3e−10 relative tolerance.

Numerical convergence of the **forward** 2D weak Qφ against 1D analytical reference, at three (parent radial, angle) grids:

| pmax | 24×12 | 48×24 | 72×48 |
|---:|---:|---:|---:|
| 2 | 2.64e−8 | 2.96e−12 | 2.95e−12 |
| 8 | 7.03e−5 | 8.69e−9 | 9.46e−13 |
| 16 | 7.73e−5 | 2.04e−8 | 1.13e−11 |

These are absolute relative errors, NOT statistical confidence intervals. Increasing the parent cutoff from 2 to 16 changes the forward number rate by ~74.6%, so a single pmax=2 agreement is not UV convergence. The chemical balance test at μφ=−0.2 GeV regulates the extreme Bose IR and does NOT test near-critical μφ→mφ, Bose condensation or the physically needed zero-mode.

## Reproduction, checks and decision

- [Standalone source](phi_bose_chemical_weak_reference.py), [three regression tests](../tests/test_phi_bose_chemical_weak_reference.py). Commands: \`OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_chemical_weak_reference.py\` and \`OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_chemical_weak_reference\`.
- **Local independent focused tests: 3/3 passed** in 0.36 s in a Python/SciPy environment (direct k-integral, reverse/equilibrium/entropy, three pmax×three angle/radial grids).
- These tests verify a **physically normalized weak prescribed-spectrum** reference only. They are NOT proof of a physical finite-volume strong operator: daughter phase-space cell capacities, connected event graph, finite-step positivity/entropy, long-time microscopic rates, realistic SM, FRW, Ωφ, BEC and CLASS remain open.

**NEXT_ACTION:** apply both zero-μ and nonzero-μ number/energy weak benchmarks to any proposed positive shared-capacity χ1↔χ2+φ strong operator, including reverse-dominated and exact-chemical-equilibrium cases. Do not promote the arbitrary resonant synthetic lattice to physical kinetics; independently prove cell DOS normalization, shell microreversibility, graph nullity, time convergence and IR convergence first.
