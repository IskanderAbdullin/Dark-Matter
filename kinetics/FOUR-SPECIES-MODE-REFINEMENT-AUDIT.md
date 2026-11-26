# Fermi-logit finite-energy-mesh convergence in coupled χ1χ2–φ–lepton 2↔2 (2026-11-22)

## Hypothesis and independent gate

The immediately preceding [orthogonal phase-space-cutoff/mode-mesh audit](FOUR-SPECIES-CUTOFF-MESH-SEPARATION-AUDIT.md) found that an energy-mode grid with 256 cells still underestimates the *analytic-occupation, same-quadrature* relativistic 5D Pauli scattering entropy by 5.73% under physically admissible oscillatory fermion occupations. Further energy-mode refinement is necessary before a strong time history can be called a numerical solution even of this fixed synthetic kinetic model.

Prespecified check: derive a mathematical bound for affine interpolation of fermion log-odds; independently verify it for three differently curved χ1, χ2 and shared-lepton logit functions; evaluate four or more successively halved uniform physical energy-mode spacings at **unchanged** 5D radial, relative-angle, CM log-angle, azimuth and incoming cutoff, using a separately computed non-gridded-occupation weak integral for Qχ1, Qχ2 and collisional H. Do not confuse spatial mode-grid order with full 5D quadrature order, nor assume signed Qχ errors converge monotonically.

## Exact interpolation remainder

For a twice continuously differentiable logit ℓ(E)=ln[f(E)/(1−f(E))] on [E_i,E_{i+1}] with spacing h, the ordinary linear interpolation error satisfies

    |ℓ(E)−I_hℓ(E)| ≤ h²/8 * sup_{[E_i,E_{i+1}]} |ℓ''(E)|.

For the declared synthetic fermion occupations:
    ℓχ1(E)=−E/1.1+0.3 sin(4E), sup|ℓχ1''|=4.8;
    ℓχ2(E)=−E/.9−0.2 sin(3E), sup|ℓχ2''|=1.8;
    ℓℓ(E)=−E+0.25 cos(4E), sup|ℓℓ''|=4.0.

The logistic expit derivative is ≤1/4, so the absolute occupation error due solely to logit interpolation also obeys |Δf|≤h² sup|ℓ''|/32. This is an analytic interpolation bound, not a certified heat/entropy integral error. The test validates direct node-logit interpolation on ≥10,000 deterministic energies at N=64,256,1024 for all three species, with no fitted constants.

## Measured conditional strong-to-weak rate error

Hold the **shared physical energy mesh endpoint** at mχmax + 2cut = 40.1 and the **incoming physical cutoff** at 12, with relativistic full finite-vector χℓ Pauli scattering and 1 common evolving lepton distribution. Mχ1=.1, Mχ2=.02, mℓ=0, G=.01, mV=.1; p/k segmented 12×3, nμ6, outgoing mapped CM ncos12, φ-azimuth n12. All are synthetic dimensionless input conditions; the separate real-φ Yukawa operator remains present in the coupled model, but the rate comparison here isolates *elastic* instantaneous weak Q/S.

The independent analytic-Fermi-occupation *same-quadrature* 5D weak reference is Qχ1=−1.2672292278252052e−11, Qχ2=+8.75060244836362e−12, Hsum=+1.1236203110330876e−10. For the identical quadrature using **Fermi-logit energy node interpolation + paired conservative deposition**, compare:

| energy intervals per species | spacing h (χ1) | Qχ1 relative error | Qχ2 relative error | collisional H relative error |
|---:|---:|---:|---:|---:|
| 256 | 0.156250 | −0.5030% | −0.2297% | **−5.7321%** |
| 512 | 0.078125 | −0.0192% | −0.0745% | **−1.4830%** |
| 1024 | 0.039063 | −0.0178% | +0.0049% | **−0.4003%** |
| 2048 | 0.019531 | +0.0021% | −0.0056% | **−0.0966%** |
| 4096 | 0.009766 | −0.00054% | −0.00151% | **−0.02447%** |

When halving h, successive magnitudes of H error shrink by measured ratios **3.865, 3.705, 4.146, 3.946**, consistent with order ≈2 **for this fixed smooth-logit interpolation and quadrature**. The analytic χ1 logit remainder bounds at these h are respectively 0.01464844,0.00366211,0.000915527,0.000228882,0.000057220. All tested node-event entropy identities and stoichiometric/energy weak-source residuals stay at ≲10⁻²³ absolute.

**Negative heat control:** relative Qχ2 error *changes sign* between 512→1024→2048; Qχ1 error also changes sign between 1024→2048. Therefore an apparently monotone entropy convergence trend does NOT imply that each signed scattering heat has a monotone Richardson extrapolation at the same finite incoming quadrature. The 4096-mode H error below 0.03% is *only* the Fermi energy-logit interpolation component of the full uncertainty; radial, angular, cutoff, physical lepton chemistry and time evolution could dominate.

## Independent cutoff tail at fixed fine energy mesh

To avoid falsely attributing missing high-energy modes to phase-space truncation, hold the **same 1024-bin mesh cutoff20** and compare incoming physical cuts8,12,16,20 with a fixed maximum radial segment width ΔE=1. For the oscillatory occupation:

| physical incoming cut | projected χ1 Q | projected χ2 Q | projected entropy H |
|---:|---:|---:|---:|
| 8 | −1.186766906e−11 | +8.484575879e−12 | 1.075224478e−10 |
| 12 | −1.267003978e−11 | +8.751029775e−12 | 1.119122391e−10 |
| 16 | −1.272106772e−11 | +8.768388704e−12 | 1.121178920e−10 |
| 20 | −1.272268463e−11 | +8.768695981e−12 | 1.121258725e−10 |

From cut16→20, the χ1 heat change is ≈0.0127% on this fixed 5D angular/radial rule; from cut8→20 it is ≈7.2%, so cut8 was meaningfully incomplete. Still no convergence certificate across nlocal/nμ/ncos/nφ here; do NOT promote the apparent cutoff stabilization into a physical γ_el.

## Reproducibility, negative/positive scope

- [Finite-grid refinement runner](four_species_mode_refinement_gate.py), [three focused tests](../tests/test_four_species_mode_refinement_gate.py), [versioned machine output](results/four-species-mode-refinement-2026-11-22.json), [dated stage](../archive/stages/2026/11/2026-11-22-four-species-mode-refinement.md).
- Reproduce: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_four_species_mode_refinement_gate.py (3 local tests passed, ~13.84s); OPENBLAS_NUM_THREADS=1 python -m kinetics.four_species_mode_refinement_gate --full.

**Scientific progress, not cosmological discovery:** a mathematically derived interpolation error bound and conditional near-second-order convergence of true elastic entropy toward a separately evaluated analytic-f 5D weak integral. **Remaining:** actual simultaneous strong χ1,χ2,φ/SM time convergence, lepton species/gauge spin counts, full 5D quadrature error and explicit φ IR/zero Bose occupation, Friedmann/redshift, parameter viability, 86-eV ΩDM/BEC/frozen f(q) and CLASS. Preserve CLASS-ready=false.
