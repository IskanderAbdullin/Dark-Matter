# Endpoint-weight error in full Z′ energy-shell quadrature: independent weak comparison (2026-11-22)

## Hypothesis and controlled comparison

The exact daughter angular-cell intersection in
[VECTOR-NEUTRINO-SHELL-CELL-OVERLAP-AUDIT.md]
(VECTOR-NEUTRINO-SHELL-CELL-OVERLAP-AUDIT.md)
normalizes each vector's flat angular decay probability to
floating accuracy, but its continuum vector weak ENERGY source
remains biased by +2.229% (h=.2) and +0.592% (h=.05).
The vector energy quadrature previously assigned the
**full h** density-of-states weight to the boundary
energy E_j=Emax, whereas the independent reference
integrates E only up to Emax. This may cause
a leading-order boundary error that hides true
angular accuracy.

Predeclared test: compare three **identical energy/daughter
center grids and exact angular overlap kernels**:
(A) original vector midpoint wV_j=gV E_j p_j h/(2π²);
(B) standard half-cell/trapezoid upper endpoint weight;
(C) positive exact mass-dependent radial 0th density
weight integrated over the actual E cell clipped
to [mV,Emax]. Do not change microscopic quantum
F−R, degree-of-freedom counting, shell energies or
implicit bounded event dynamics. Compare N/E/S
against independent continuum 144×144 Gauss weak
integral for h=.4,.2,.1,.05 and retain negative tests.

## Exact radial 0th phase-space measure

For spin-1 massive vector (gV=3) and energy E≥mV,

    dNV phase space = gV/(2π²) E sqrt(E²−mV²) dE.

The primitive is

    ∫ dE E sqrt(E²−mV²)
       = (E²−mV²)^(3/2)/3.

For vector midpoint E_j=jh use
L_j=max(mV,E_j−h/2),
U_j=min(Emax,E_j+h/2), and

    wV,j^exact = gV/(6π²)
       [ (U_j²−mV²)^(3/2) − (L_j²−mV²)^(3/2) ].

Only vector sites E_j>mV are included; the
lowest threshold cell at p=0 is not represented
and its truncated contribution remains an
identified approximation. The full angular
weights W_jil are scaled by **the same positive
ratio wV,j^exact / wV,j**, preserving reaction
conservation and exactly normalized shell sum.
The upper-trapezoid scheme instead multiplies
only the Emax shell by 1/2.

These alternatives adjust vector *number*
density-of-states quadrature, but still deposit
the energy of an entire E bin at E_j. They are
**not** exact physical cell-integrated heat
moments or a fully converged strong operator.

## Independent numerical results

Use the same prescribed independent nonthermal
fν, fantiν and 0.75 BE vector fixture as the prior
stage, mV=.4, Γ0=T=1, Emax=4. Relative errors
below use a **separate continuum p,CM-angle
Gauss 144×144** reference:
dNV=0.0058088005769424,
dEV=0.0115779873168571,
entropy=0.00144991988902934.

| h | Original full-endpoint E error | Upper half E error | Exact radial-cell E error | Exact-cell N error |
|---:|---:|---:|---:|---:|
| .4 | +4.2202% | −0.7938% | −0.6307% | −2.6938% |
| .2 | +2.2292% | −0.2765% | −0.2288% | −1.0101% |
| .1 | +1.1547% | −0.0980% | −0.0866% | −0.3813% |
| .05 | +0.5922% | −0.0342% | **−0.0320%** | **−0.1417%** |

The E-source error at h=.05 shrinks from +0.5922%
to −0.0320% in magnitude, an approximately
**18.5-fold reduction**, while the exact angular
shell probability remains 1 per vector. Number
error is not eliminated: at h=.05 it is −0.1417%,
similar in magnitude to the prior +0.1640%, but
its sign changes. Entropy error also becomes
−0.0706% instead of +0.2816%.
The three small finite-h exact-cell E errors
h=.2→.1→.05 decrease in magnitude
~0.2288%→0.0866%→0.0320%, consistent
with significant boundary-bias removal, not
a proof of universal convergence order.

**Independent negative finding:** removing the
upper-boundary overcount does not guarantee
N-source or entropy accuracy; the lower
mV threshold cell, daughter k weighting,
energy-center collocation, angular variation
within a cell and tail Emax remain.
The coarse h=.4 number error actually
worsens from −0.7119% with the original
weights to −2.6938% with clipped exact cells.
Thus more exact individual quadrature factors
do not imply monotonically better all-moment
errors for a fully discretized kinetic system.

## Microreversible conservation and reproducibility

Source:
[vector_neutrino_shell_radial_endpoints.py]
(vector_neutrino_shell_radial_endpoints.py).
Tests:
[test_vector_neutrino_shell_radial_endpoints.py]
(../tests/test_vector_neutrino_shell_radial_endpoints.py).
[Compact JSON]
(results/vector-neutrino-shell-radial-endpoints-2026-11-22.json).
[Permanent stage]
(../archive/stages/2026/11/2026-11-22-vector-neutrino-radial-endpoints.md).

    python -m pytest -q tests/test_vector_neutrino_shell_radial_endpoints.py
    python -m kinetics.vector_neutrino_shell_radial_endpoints

Tests demand independent energy-source improvement at
three resolutions, a **remaining nonzero number error**,
per-vector flat angular normalization and exact
eventwise C,L,E and quantum entropy after a positive
implicit step. These are scoped math/engineering
checks; real GitHub CI must be separately confirmed.

## Remaining scientific gates

This *deterministic reduced model* uses synthetic
nonthermal distributions, Γ0=1 arbitrary time unit,
fixed mV and a truncated set of Z′↔ννbar events.
It does not solve the actual thermal ν-e±γ
plasma, gauge-neutral real φ + χ1/χ2 kinetic
production, χ elastic Γel, FRW entropy/H,
full physical energy dispersion/quadrature
across all modes, or zero-mode BEC.
Precise 86-eV f(q)/ΩDM/CLASS-ready cannot
be inferred. Continue with independent
angular/cell integrated heat, cutoff/timestep
convergence and then the physical gauge/SM
model. All negative cases stay archived.
