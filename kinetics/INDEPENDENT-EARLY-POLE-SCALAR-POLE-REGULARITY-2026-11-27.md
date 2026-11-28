# A stepped approximation can get the local pole curvature sign wrong

The [continuous conserved ramp](INDEPENDENT-EARLY-POLE-SCALAR-CONTINUOUS-POLE-RAMP-2026-11-27.md) and the 8/16/32/64-subcell approximations represent the same intended local rounding, while preserving the native two-cell number and energy in each case. For the squared-pole piece alone, measure the even finite difference G(x)=F(x)+F(−x)−2F(0) at matched invariant offsets; report G/x and G/x² rather than a one-sided fit. The reconstructed archive is fully fixed during each fold.

At q=4.1e−6 on fine, in the exterior-partner net at x=1e−6:

| Occupation representation | G/x², angular order 18 | G/x², upper angular order |
| --- | ---: | ---: |
| Original archive steps | +5.196e−16 | +5.192e−16 (order 32) |
| Ramp approximated with 64 subcells | +9.587e−15 | +9.701e−15 (order 32) |
| Exactly continuous ramp | −8.953e−14 | −8.954e−14 (order 28) |

The 64-cell ramp also gives a positive value at x=2e−6 with both orders. For the exact ramp the finite difference divided by x² stays negative and nearly constant over those offsets. The crude local curvature of the *staircase approximation* therefore does not approximate the local curvature of the intended continuous profile at those offsets, even when a quadrature order test on the staircase agrees in sign. At larger offsets the stepped sign can change. This is a comparison at fixed tested offsets only: it does not assert a mathematical noncommutation of limits for every sequence, or convergence of the archived profiles to a physical continuum state. Raw samples below the resolved ranges are retained as a numerical cancellation warning, not interpreted.

To interpret the width convention, a unit cusp G=x gives the local constant-subtracted Breit-Wigner integral ½ log(1+R²/gamma²) over [0,R]. The finite principal-value-squared kernel gives the same logarithm plus gamma²/(R²+gamma²)−1. A smooth quadratic G=x² instead gives R−gamma atan(R/gamma) for Breit-Wigner; the difference from the principal-value result vanishes with gamma. Independent log-coordinate quadrature tests these closed expressions. Subtracting the *global* constant delta rather than the *local window area* also changes the central integral by −2 atan(gamma/R)/gamma per unit pole constant. These manufactured checks separate regularity and bookkeeping; they do not identify a real width, a physical cusp, or a new collision term. The finite-width kernel conventions are discussed by [Ala-Mattinen et al.](https://doi.org/10.1103/PhysRevD.108.096034).

From `kinetics/`:
```sh
for g in fine third; do
 OPENBLAS_NUM_THREADS=1 python -W error early_pole_archive_scalar_pole_regularity_reference.py "$g"
done
```

Raw for both archived/stepped profiles: `results/early-pole-archive-scalar-pole-regularity-2026-11-27.json.gz`. The previously published continuous-ramp artifact contains the corresponding exact profile values.
