# Finite-width FD inverse portal and reduced gain duration

Combines exact vacuum vector BW with FD massless electron angular weights from the pole gate, for initially empty dark final states. No full nonlinear collision evolution. Swap positive s/E' integrals; at fixed E and invariant s, b=(s-2m²)/2, B=sqrt(b²-m⁴), partner bounds are E'_±=(E b±p B)/m², lower rationalized as (b²+p²m²)/(E b+p B). Within these limits integrate e^(-E'/T) times angular FD weight. Partner absolute cap m+60T, relative cap 50T, parent energy relative cap 60T; no certified tail enclosure. Positive segmented Gauss around finite resonance and thermal scales, adaptive outer integral, no delta function at pole. Initial electron-only source, massless e, no electron chemical potential, fixed vacuum width. Other SM and medium corrections omitted.

| T | alpha | MB full BW | FD n12 | FD n16 | FD n24 |
|---|---|---|---|---|---|
| 0.4 | 0.012 | 1.002222929335 | 1.002195334795 | 1.002195334782 | 1.002195334777 |
| 0.6882487043 | 0.012 | 1.000204681492 | 1.000193199536 | 1.000193199550 | 1.000193199550 |
| 5 | 0.0001 | 1.000059779585 | 1.000063697317 | 1.000063697316 | 1.000063697316 |

The independent swapped MB ratio agrees with previous direct BW ratios within 8e-9 absolute. FD ratios n12/16/24 agree to <1.5e-11 absolute in these three tests; this is quadrature sensitivity, not rigorous error bounds. FD changes inversion by -2.76e-5 at .4 GeV, -1.15e-5 at .688 GeV, +3.92e-6 at 5 GeV. Pole-only FD at .4 GeV gave ratio1.0001757067, full FD1.0021953348, confirming off-pole importance.

Apply previous analytic reduced duration bound to this FD injection. Conditions still fixed T/source, spontaneous chi1 decay, no chi2 depletion/inverse/Bose feedback. Eight samples:

| T | alpha | J1/J2 | positive gain exponent upper (reduced) |
|---|---|---|---|
| 0.4 | 0.012 | 1.002195334782 | 3.200970e-13 |
| 0.6882487043 | 0.0001 | 0.995343512294 | 0.000000e+0 |
| 0.6882487043 | 0.001 | 1.000088693361 | 1.340868e-11 |
| 0.6882487043 | 0.012 | 1.000193199550 | 9.793591e-13 |
| 1 | 0.001 | 1.000140123241 | 4.063957e-9 |
| 2 | 0.001 | 1.000092522522 | 1.029757e-7 |
| 5 | 0.0001 | 1.000063697316 | 2.773468e-5 |
| 10 | 0.0001 | 1.000035130448 | 7.838068e-6 |

At T=.6882487043 alpha=.0001 source inversion is negative. Largest sharpened upper over these eight FD samples 2.77347e-5. Not a continuous parameter maximum and not a bound for the full coupled dark network. Physical depletion of chi2 can alter the conclusion; proceed to actual abundance/depletion dynamics rather than assuming selective mother reservoir.

Reproduce python -m kinetics.phi_bose_portal_fd_bw_gate. Locally executed. No remote CI claim or physical BEC/relic/CLASS.
