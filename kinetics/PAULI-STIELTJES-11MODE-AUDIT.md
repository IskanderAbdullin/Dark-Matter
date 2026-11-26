# Independent Stieltjes 2..11-mode relativistic Pauli χe collision spectral audit

**Purpose:** independently extend the [2..5-mode quantum collision Galerkin audit](PAULI-RITZ-BASIS-AUDIT.md) to more nonthermal modes, without mistaking the well-resolved 5D angular/momentum quadrature for *basis-size* convergence. Microphysics is still the conditional massless electron-only four-fermion vector EFT with Mχ=.1 GeV, T=1 GeV, G=.001 GeV⁻², gχ=ge=4, electron FD bath fixed. No full kinetic spectral gap or 86-eV relic is established.

## Numerically independent basis construction

[stieltjes_pauli_modes.py](stieltjes_pauli_modes.py) obtains FD-susceptibility-weighted **orthonormal polynomials** of x=E/T through the stable three-term Stieltjes recurrence `q_(n+1)=((x−alpha_n)q_n−beta_(n−1)q_(n−1))/beta_n`, with 450-point independent thermal momentum quadrature on p∈[0,90 max(M,T)] and exact fFD(1−fFD) weight. This avoids catastrophic subtraction of high-degree monomials. Constant q0 is an **exact collision-number null direction** and excluded from the nonconserved Ritz matrix. Quantum full-Pauli Dirichlet kernel D and susceptibility S are evaluated independently of the earlier polynomial Gram–Schmidt implementation; rates from 2 and 5 nonconstant modes reproduce it to <10⁻10 relative at the same 5D collision grid.

The 5D integration still uses exact two-body CM boosts, Møller flux, Pauli F_eq and spin-averaged vector-contact dσ/dΩ. A positive 11-mode matrix is a **finite Galerkin approximation** and cannot, by itself, prove a strictly positive continuum spectral gap. At fixed basis/momentum-domain, n=20→24 is a numerical convergence check, not a basis convergence check.

## New negative result: no stable low kinetic rate through 11 modes

For the same benchmark and n=24 Gauss points per dimension, cutoff pmax=kmax=32T:

| Non-number logit modes | Lowest Ritz relaxation (GeV) |
|---:|---:|
| 2 | 2.125371383×10⁻⁷ |
| 3 | 1.725865686×10⁻⁷ |
| 5 | 1.329341851×10⁻⁷ |
| 7 | 1.111531407×10⁻⁷ |
| 9 | 9.645584196×10⁻⁸ |
| 11 | **8.564757723×10⁻⁸** |

Thus the 11-mode rate is ≈**59.70% lower than** the already nonconverged two-mode rate, and ≈**71.75% lower than** the conditional fixed-N thermal-energy Rayleigh rate γ_T=3.031821176×10⁻⁷ GeV. Meanwhile for **the same 11-mode basis**, n20→n24 changes the smallest Ritz rate from 8.561217213×10⁻⁸ to 8.564757723×10⁻⁸ GeV (≈0.0413%). That **strong contrast** demonstrates that good 5D collision quadrature does not resolve the physical slow spectrum. Strong coupling scaling all Ritz values ∝G² and basis-orthogonality ~10⁻15 provide necessary consistency checks.

This does **not** prove that the true infinite-dimensional collision operator is gapless, nor that it has an eigenvalue of 8.56×10⁻⁸ GeV: the numbers are variational Ritz estimates in nested polynomial subspaces, subject to assumed regularised FD inner product, contact EFT, limited cutoff, and absence of SM dynamics. Higher-order modes, far-UV compact support and continuum tails could continue to lower the variational bound. Proving a full gap/no-gap needs independent functional analysis and phase-space-tail assumptions.

[Source](stieltjes_pauli_modes.py), [tests](../tests/test_stieltjes_pauli_modes.py), [compact results](results/pauli-stieltjes-11mode-2026-11-21.json). The direct old phenom BGK γ=10⁵ cannot be identified with any of these conditional rates.

## Next physical research stages

1. Analyze high-energy and low-momentum localized test perturbations rather than only global polynomials; scan cutoff/basis and derive analytically controlled bounds on Dirichlet form / susceptibility, with no unwarranted continuum spectral-gap statement.
2. Map **real microscopic** χe Pauli event integrals to finite-volume fermionic capacity transitions, preserving energy and particle number; [reduced finite-bath Pauli network](FINITE-BATH-QUANTUM-PAIR-REFERENCE.md) is a conservation reference with *arbitrary, not derived* transition Γ.
3. Implement evolving finite SM radiation plus UV-gauge-consistent ν/e plasma, thermal masses and 1↔2 stimulated Bose, and only then long-time cosmological freezeout and observational 86-eV dark-matter tests. All physical CLASS/BEC/Ω claims remain false until independently validated.
