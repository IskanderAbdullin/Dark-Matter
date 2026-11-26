# Finite-cell recoil projection: converged negative detailed-balance control

2026-11-21. **Synthetic micro-event measure, not the physical 5D χe collision operator.** This is a separate consistency counterexample, following [BIN-LOCAL-DETAILED-BALANCE.md](BIN-LOCAL-DETAILED-BALANCE.md).

## Pre-registered question and reference
Does pairing a *continuous* time-reversal symmetric positive event measure, obeying exact microscopic local detailed balance at each event, automatically yield the discrete center-energy ratio exp[(Ei-Ej)/T] when both finite bins are integrated? The independent reference is exact pointwise microscopic balance, not an assumed coarse balance. Prediction: generally no; bin-center projection bias should tend to zero as width tends to zero.

Choose continuous initial/final energies E∈[Ei−h/2,Ei+h/2], E′∈[Ej−h/2,Ej+h/2], Ei=2, Ej=1, T=1 (all in arbitrary common energy units), a strictly positive symmetric event weight W(E,E′)=exp[−0.17(E+E′)](1+0.12(E−E′)²), and paired forward/reverse weights W exp[±(E−E′)/(2T)]. Their pointwise ratio is **exactly** exp[(E−E′)/T]. Integrate both over the same rectangle using tensor Gauss–Legendre n×n and form R=Kforward/Kreverse. The coarse equilibrium requirement instead demands Rcenter=exp[(Ei−Ej)/T]. This construction is an algebraic counterexample to naive projection, **not** a derivation of physical χe matrix elements, recoil kinematics, degeneracy, or rates.

## Reproducible result
Run `python kinetics/finite_bin_recoil_projection_control.py` (NumPy required). Relative defect R/Rcenter−1 at n=32:

| h | defect |
|---:|---:|
| 0.8 | +0.0227700461053 |
| 0.4 | +0.00570907623235 |
| 0.2 | +0.00142824952813 |
| 0.1 | +0.000357122796190 |
| 0.05 | +0.0000892844613802 |

At every width, n=8,16,32,64 agree within 2×10⁻¹²; halving h reduces the defect approximately fourfold. Thus quadrature convergence alone **does not** cure finite-bin detailed-balance error. A center-energy FD distribution then yields spurious pair current unless the coarse transition rates are constructed with a declared equilibrium closure or exact energy resolution.

## Next decisive physical experiment
Replace synthetic W by actual χe Lorentz-invariant vector-contact 5D event weight, including χ and electron state capacities, electron FD occupation and final Pauli, exact on-shell recoil, paired forward/reverse sampling, and adapted IR radial cells. Compare raw projected rates, center-balance corrected rates, independent strong/weak heat, and grid-width extrapolation; quantify distortion introduced by correction. Only then construct Pauli-bounded finite-bath evolution. No physical Γel, 86-eV Ω, BEC or CLASS readiness is established here.
