# Galerkin basis-size audit: apparent slow Pauli kinetic rate is NOT converged

This follow-up to [the two-mode linearised full-Pauli χe quantum collision study](LINEARIZED-PAULI-EIGENMODES.md) tests the **basis truncation dimension** separately from 5D momentum/angular quadrature. The kernel remains the conditional electron-only vector contact EFT with χ mass=.1GeV, common T=1GeV, μ=0, G=.001GeV⁻², spin degeneracies gχ=ge=4. This is a numerical model sensitivity test, not a cosmological DM prediction.

## Independent construction and falsification

Use a second implementation [higher_order_pauli_modes.py](higher_order_pauli_modes.py): orthonormal polynomials of x=E/T weighted by the physical χ susceptibility w(p)=gχ p²fFD(1−fFD)/(2π²). Twice-reorthogonalized polynomial Gram–Schmidt, 280-node independent physical-p quadrature, produces h0=constant plus h1..h5, all number-orthogonal and mutually orthonormal. For each degree, compute **the same actual 5D full-Pauli positive quadratic form** D_ab from the tree χe collision element and outgoing azimuth; generalized Ritz matrix gives approximate relaxation rates. Discrete S identities and h0 null checked independently. n2 matches the previously **different direct two-mode analytic centered basis** to ≲1e−10 relative on the same CM grid. G→2G scales *all* Ritz rates 4-fold.

The Rayleigh-Ritz min–max theorem predicts the smallest Ritz eigenvalue **must not increase** as the trial space expands. This is not a statement that it converges to a nonzero positive physical spectral gap without proof of a coercive bound or full operator analysis.

## Actual numerical negative result

<escape>For n=24 Gauss points in each of five collision dimensions and cutoff p,k=32T, the smallest generalised Ritz relaxation rates [GeV] are:</escape>

| Nonconserved logit modes | Smallest Ritz rate, GeV |
|---:|---:|
| 2 | 2.125371382675×10⁻7 |
| 3 | 1.725865685791×10⁻7 |
| 4 | 1.490332342631×10⁻7 |
| 5 | 1.329341851403×10⁻7 |

This fifth-mode value is **37.45% below the two-mode estimate** and about **56.15% below** the original one-dimensional thermal-energy rate γ_T≈3.03182×10⁻7 GeV. So the two-mode slow Ritz value from the previous paper-style report, even after n/cutoff refinement, **must NOT be called the "slowest kinetic relaxation rate"**. It is only the lowest eigenvalue *within its two-mode trial subspace*.

At five modes and cutoff32, n20 gives 1.329598340929e−7 and n24 1.329341851403e−7 (a relative 1.93×10⁻4 shift); n12/cutoff32 gives 1.38540395252e−7 (4.2% higher). Thus collision-quadrature convergence is largely under control at n20→24, while **basis-size convergence is manifestly not achieved**. The underlying full kinetic eigenvalue could be lower still, and an actual spectral gap may be zero under continuum assumptions. This is a genuinely limiting negative result, not a number to fit to a Bose dark-matter scenario.

Source: [code](higher_order_pauli_modes.py); [tests](../tests/test_higher_order_pauli_modes.py) covering basis invariants, previous independent two-mode comparison, monotone Ritz behaviour, n refinement and coupling scaling; [compact measured results](results/pauli-ritz-basis-scan-2026-11-21.json). No full Boltzmann evolution or finite SM-bath γ established.

## Next scientific phase

Investigate **eigenvalue basis-size asymptotics** with better adapted modes (IR/UV weighted functions, localized momentum packets and physically supported χ energy perturbations) to determine whether slow tails lead to arbitrarily weakly relaxed modes. Report any evidence of no strictly positive spectral gap with precise assumptions, but do NOT claim it after five modes. Independently time-integrate a small source conserving particles and Pauli bounds, verifying whether logarithmic/exponential thermalisation modes from Ritz analysis actually control the slow response. Gauge-complete SM/ν and UV-EFT/mediator must be treated before any realistic FRW 86-eV relic.
