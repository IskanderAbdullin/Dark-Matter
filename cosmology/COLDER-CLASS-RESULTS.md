# Verified CLASS scan: artificial colder 86-eV Bose spectra with exactly the same number density

**Run:** [successful public CLASS job 44890489385](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44890489385), 2026-11-21, CLASS wrapper classy 3.4.3.0. Both the baseline thermal 86-eV scalar and the toy four-case shape scan were actually run in the Einstein–Boltzmann solver, with complete raw JSON outputs in the run artifact. This file permanently preserves primary high-level results and their assumptions.

## Model — no claim of physical production

Start from the baseline ONE REAL scalar relativistic thermal Bose relic (m=86 eV, gstarS decoupling 55, chemical potential zero before decoupling). Define a purely illustrative family:

\[
f_\alpha(p)=\alpha^{-3}
\left\{\exp[p/(\alpha T_b)]-1\right\}^{-1},
\qquad 0<\alpha\le1.
\]

Here T_b is the reference boson momentum-temperature scale. The factor alpha^(-3) holds the *present particle number* and thus omega h² constant, while momentum coldness is reduced by alpha and early radiation density by alpha. This is **NOT** an equilibrium Bose distribution, an actual stimulated freeze-in calculation, evidence of an early condensate, or a microscopic dark-sector Lagrangian. Huge mode occupancy at small alpha may invalidate neglecting collective/collisional dynamics, even if a collisionless CLASS benchmark is mathematically definable.

In the CLASS file convention (with one scalar quantum-cell factor included), f0(q)=alpha^(-3)/[(2pi)^3(exp(q)-1)], q=p/(alpha*T_b), T_ncdm=alpha*T_b/Tgamma, and deg_ncdm=1. No omega_ncdm is forced; CLASS infers it from the PSD.

## Actual linear Einstein–Boltzmann outcomes

One fixed cold reference at present omega_nonbaryon*h²=0.11912939486 (plus omega_b=0.02237) and the same primordial amplitude/tilt is compared with each collisionless noncold scalar. CLASS returned omega_ncdm*h²=0.11913012078 for each alpha to the precision shown. The early extra effective radiation fraction follows exactly the alpha scaling in the code.

| Toy alpha | CLASS N_eff early | Extra Delta N_eff | k at T(k,z=5)=0.5, 1/Mpc | T(k~1/Mpc,z=5) |
|---:|---:|---:|---:|---:|
| 1 (thermal) | 3.10881512 | 0.06481512 | **0.43275347** | 0.010281 |
| 0.5 | 3.07640756 | 0.03240756 | **0.75023757** | 0.307512 |
| 0.3 | 3.06344454 | 0.01944454 | **1.13675861** | 0.618455 |
| 0.1 | 3.05048151 | 0.00648151 | **2.87030768** | 0.949297 |

Here T=sqrt(P_nonthermal/P_CDM) denotes linear **matter fluctuation amplitude**, not the boson temperature. If T=0.5, linear power is P/P_CDM=0.25. The k~1 row refers to the nearest actual sampled k=0.9862059 inverse Mpc.

At z=0 the corresponding amplitude-half points are also computed and archived in the JSON; the z=5 values are shown because they are directly relevant to early structure formation.

**Meaning:** simply decreasing primordial momentum by a factor ten at fixed m=86 eV and present abundance shifts the benchmark half-amplitude scale from k≈0.43 to k≈2.87 inverse Mpc. This is substantial but still represents a finite suppression scale; whether it meets high-redshift galaxy/Lyman-alpha/CMB constraints cannot be decided without the relevant likelihood and nonlinear baryonic modeling. More importantly, no physical mechanism for creating this alpha-scaled Bose occupation has been established.

## Verification and limitations

1. Analytic number conservation and early-radiation scaling have been independently integrated over the positive PSD, with <1e-4 relative moment error.
2. The same CLASS source and numerical convention as the independently converged alpha=1 thermal case were used. CLASS background omega_ncdm is constant across all four experiments and N_eff follows the alpha scaling to <1e-7 absolute.
3. The alpha=1 CLASS half-amplitude crossing agrees with the independently rerun thermal benchmark to within a few parts per million.
4. The original thermal model has an independent coarse/fine q-grid comparison near its k-half crossing (0.037% change). **The alpha=0.1–0.5 spectra have not yet undergone the same independent q-grid-resolution and CLASS precision study**, especially at deeply suppressed k, and should be treated as preliminary benchmark outputs.
5. This is a single-species collisionless initial-PSD exercise. A physical theory needs a Boltzmann production kernel for f(q), internal conserved charges, quantum scattering/condensate dynamics and an actual formation/evolution mechanism for late galactic near-critical noncondensed Bose halos.

## Next theoretical computation

A compelling microphysical program is to take the momentum-dependent freeze-in mechanisms of Yin (2023) and Sakurai–Yin (2025), compute their actual normalized f(q) including inverse processes and Bose enhancement, pass those spectra unchanged through CLASS, and independently derive whether any late GP–Boltzmann–Poisson halo has time to reach the JCAP 2021 near-critical f(E) and sustain a finite tidal boundary.

Do not interpret a toy alpha that happens to resemble a cold spectrum as a derivation from those mechanisms. Exact sources:
- CLASS official https://github.com/lesgourg/class_public/blob/master/explanatory.ini and https://github.com/lesgourg/class_public/blob/master/source/background.c .
- Benchmark code in this repository [cold_psd_scan.py](cold_psd_scan.py), [thermal_bose_class.py](thermal_bose_class.py).
- Standard thermal WDM Lyman-alpha constraint for a *different* benchmark: https://doi.org/10.1103/PhysRevD.109.043511 (5.7 keV lower mass under its own assumptions; cannot transpose directly to this bosonic PSD).
- Stimulated production references https://doi.org/10.1007/JHEP05(2023)180 and https://doi.org/10.1007/JHEP03(2025)202 .
