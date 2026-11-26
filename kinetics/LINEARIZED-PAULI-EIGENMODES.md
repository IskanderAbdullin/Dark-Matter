# Two nonthermal linearized full-Pauli collision modes, not a full kinetic spectrum

**Scope (2026-11-21):** the χe contact EFT, massless electron/positron FD bath, 5D exact-CM scattering with both χ and electron final-state Pauli, scalar spin-averaged |M|² of the already tested [portal](MICROPHYSICS-MODEL.md). [Code](linearized_pauli_modes.py), [independent tests](../tests/test_linearized_pauli_modes.py), [numerical benchmark](results/linearized-pauli-two-modes-2026-11-21.json).

**Question:** The previously derived single fixed-N thermal-energy moment has γ_T≈3.0318×10⁻7 GeV. Is that a universal kinetic relaxation rate for *all* nonthermal χ perturbations? The hypothesis "all modes relax at γ_T" is falsifiable with an independent two-mode quantum collision Dirichlet form.

## Quantum equilibrium linearization, unit-normalized

Fix a common T of both electron and χ equilibrium f=FD((E−μ)/T), and introduce infinitesimal χ logit perturbations δℓχ(p)=∑_a ε_a h_a(E), where ℓ=ln(f/(1−f)). Then δf=f(1−f)δℓ. Elastic χ number conservation has exactly one null direction h0(E)=1. Define the symmetric susceptibility

```text
S_ab = gχ/(2π²) ∫_0^∞ dp p² f(1−f) h_a(E)h_b(E)  [GeV³].
```

The matrix of quantum dissipation obtained by linearising the *same* forward-minus-reverse 5D Pauli χe contact kernel is

```text
D_ab = (gχ ge/2) ∫ d³p d³k/(2π)^6 v_M
                 ∫dΩ* [dσ/dΩ*] F_eq
                 [h_a(E')−h_a(E)][h_b(E')−h_b(E)]  [GeV⁴],

F_eq = fχ(E) fe(k)[1−fχ(E')][1−fe(k')].
```

Here dσ/dΩ* already contains the tree spin-average; gχ ge count initial internal states exactly once. The χe collision operator linear response defines `S dε/dt=−D ε` **only after Galerkin projection** onto the finite basis. D is symmetric positive semidefinite, and D(1,h)=0 exactly before quadrature. We numerically solve generalized symmetric eigenproblem Dv=λSv. Its eigenvalues are **Ritz values restricted to this basis**, not the exact slowest eigenvalue of the complete infinite-dimensional Boltzmann kernel. Adding basis functions can produce slower Ritz modes; numerical eigenvalue convergence over *q grids alone* does not prove spectrum convergence.

Two particle-number-orthogonal isotropic trial basis functions (using equilibrium susceptibility weight w=p²f(1−f)):

```text
x=E/T, mean=<x>_w, var=< (x−mean)² >_w,
h1=x−mean,
h2=(x−mean)²−var−[< (x−mean)³ >_w/var]*(x−mean).
```

Then `S12=0`, `S11,S22>0`, and both basis functions are orthogonal to number h0. An independent continuous q integral checks those orthogonality conditions.

## Numerical benchmark and independent nonlinear audit

**Reproducible model**, *not* cosmological parameter inference: Mχ=.1 GeV, Tχ=Te=1 GeV, μχ=μe=0, G=.001 GeV⁻², gχ=ge=4, electron-only massless bath, cutoff pmax=kmax=32T. With n=24 nodes **in each** p,k,µ,uCM,φCM dimension, the two-mode generalized rates are

```text
λ_slow,Ritz = 2.1253713826764962×10^-7 GeV
λ_fast,Ritz = 7.323034821208818×10^-7 GeV
γ_T,Rayleigh = D11/S11 = 3.031821176217816×10^-7 GeV.
```

The n=20/cutoff32 rates are respectively 2.12535755×10⁻7 and 7.32304192×10⁻7; n=16/cutoff28 gives 2.12563771×10⁻7 and 7.32301225×10⁻7. The Rayleigh γ_T differs from the independent fixed-N FD heat capacity and ∂Q/∂T estimate by <2×10⁻5 on matched n24/cutoff32; the original nonlinear **five-dimensional full-Pauli** heat source finite-differenced under ±0.001 χ logit h1 agrees with −T D11 to <1×10⁻5 on n12/cutoff24. That checks the Dirichlet matrix normalization against a separate implementation of the nonlinear collision integral, not merely by its positive-square construction.

**Scientifically important negative result:** the 2D-space slow Ritz rate is **29.90% lower than** the one-dimensional thermal energy Rayleigh γ_T. A scalar *single energy-moment rate* therefore cannot be identified with even this two-mode nonthermal relaxation spectrum. It cannot be inserted as universal `gamma_el(p,T)` in the old BEST BGK without an explicit many-mode calibration. In this same restricted basis the off-diagonal D correlation is ≈0.447, quantitatively demonstrating coupling between thermal and nonthermal distortions.

The positivity of D is a rigorous property of the provided equilibrium integrand given microreversibility; a positive numerical two-mode D is **not** proof that a full finite-volume time integrator preserves 0≤f≤1 or entropy. This model has no finite T-dependent g*, self-consistent SM bath, physical high-T UV mediator, ν/e±/muons beyond the e-only kernel, or dynamical 1↔2/86-eV freezeout. No physical DM relic, BEC or CLASS-ready is claimed.

## Next gate

1. Enlarge the orthonormal logit basis (third/fourth IR and hard-tail modes) and estimate the low eigenvalue *basis truncation* error, separately from five-dimensional quadrature.
2. Verify charge-null eigenvector and entropy quadratic form of the independently evaluated full strong Cχ on finer q grids.
3. Implement true finite-volume positivity-preserving quantum collision dynamics with finite electron/SM heat, no reset; then couple the exact Yukawa 1↔2 term only after gauge/UV and time-dependent cosmology validated.

Reproduction: `python -m pytest -q tests/test_linearized_pauli_modes.py` and `python kinetics/linearized_pauli_modes.py`; see the JSON for detailed n/cutoff scan.
