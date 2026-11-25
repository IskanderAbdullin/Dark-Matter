# Matched-width independent kinetic kernels: numerically converged candidate
2026-11-25. Continuation of PARENT-GRID-AND-ANGLE-2026-11-25.md.

## Completed numerical evidence
Hold M1=86 TeV, M2/M1=.9, mphi=86 eV, x0=.1, kappa_SM1=47, kappa2=.03, g=4 (symmetric two-pair interpretation), published Ghisoiu-Laine rates and bath energy conservation. All results below omit scalar self-interactions and dark thermal self-energies.

At kappa=.006, full parent-angle solutions through x=80 yield Omega h^2=.11543645, .11527558, .11515891 at n=240,320,480. At n240 changing the two heavy grids to uniform gives .11543869. The separate-rectangle n120 full result is .11419809. Its n240 and 320 explicit full integrations hit their 600 s budgets at x=13.79224 and 11.36040; archive them as incomplete, not relic predictions.

At kappa=.00685:
- Parent-angle full n240: .1203497953; n320: .1202148353. The fixed width was chosen from the preceding .006 and .0072 runs; no mesh-dependent abundance rescaling.
- Independent separate-rectangle same-grid hybrid n240: .1197704085; n320: .1198516658. Local handoff at x=2.5, relative scalar number source mismatch .01427% and .01433%, parent shape TV ~1.84e-5 and temperature change ~4.3e-8. The handoff matches total energy and retains identical daughter grid and collision operator, unlike prior cross-kernel closure.
- The full and hybrid methods agree to <0.5% at matched width; mesh drifts <0.12%. Hot number fractions q/T_actual>.1: full .00338425/.00340495; hybrid .00345041/.00344541.
- A stricter full n240 calculation (no occupation floor, rtol=2e-6 and fermion atol=1e-13, exact dense Jacobian) gives .1203496203. It agrees with approximate sparse Newton Jacobian run to 1.46e-6 relative. Saved time-dependent scalar PSD. Approximate Jacobian drops only the global H-factor derivative, not the H-factor in the equations; no general speedup or guaranteed positivity is claimed.
- Same-grid closure finite-difference Jacobian error at n240/320 ~3.94e-10/1.06e-9; maximum hybrid FRW residual ~4.1e-8. Targeted kernel tests 5 passed.

At kappa=.0072, independent n240 hybrid .12146194 versus full angle .12201862; n320 angle .12189247. Useful independent consistency check.

## Conditional linear transfer
CLASS-community 26 synchronous gauge, actual frozen full parent-angle PSD, positive three-component file partition, 120/240-node strategy 2 checks. No PSD abundance normalization. At z=5, k=10/Mpc, P/P_CDM=.98390730 (n240) and .98351104 (n320). Minimum sampled ratio through 20/Mpc=.96956276/.96904752. Alternate partition drift <=.18%/.058%. Exact fully resolved perturbation hierarchy has NOT been tested on this candidate; default fluid approximation is used. No observational likelihood or claim of viability.

## Priority change
This establishes a robust numerical solution of the specified vacuum-decay, scalar-collisionless kinetic model. It is not yet a physically validated particle model. Test contact self-interactions, in-medium scalar dispersion, Bose-enhanced 2<->2 redistribution, possible condensate and number-changing effects before optimizing portal rates or interpreting CLASS as the prediction of a halo model. Compare pressure/velocity transfer in interacting and collisionless limits; retain negative findings.
 
Raw outputs including incomplete runs, stricter full PSD history, and both actual CLASS spectra: results/converged-independent-kernels-2026-11-25.json.
