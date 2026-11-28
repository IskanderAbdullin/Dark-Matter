# Bound the cold parent-energy threshold over the photon-pole uncertainty

The [moment-only transverse enclosure](INDEPENDENT-EARLY-POLE-ENDPOINT-TRANSVERSE-MOMENT-ENVELOPE-2026-11-28.md) bounds the frozen transverse photon frequency using only the per-cell charged number and energy. Propagate it into a **kinematic** radial-energy threshold for the integrated cold band at k=.0015, x=.05, q in [1e−6,1.5e−5], with the scalar mass fixed to each card's value. No Pauli factors or matrix element enter the new threshold bound.

For transfer time component −u, u=ω−E_s, and spatial modulus v, the minimum parent energy E is the root

```text
F(E;u,v) = v sqrt(E²−M²) − u E − (M²−m_f²+u²−v²)/2 = 0.
```

Here u>0, v>u, and the constant in parentheses is positive throughout the band and frequency range. The latter ensures that the required radial cosine has negative sign on the threshold branch. F starts negative at E=M and is strictly increasing for E>M. Since F_u=−E−u and F_v=sqrt(E²−M²)+v, E increases with u and decreases with v. Thus the frequency and angular extrema use the two frequency endpoints and the collinear cases v=k±q. For v=k+q, increasing q decreases E. For v=k−q, increasing q increases E: at the root sqrt(E²−M²)>uE/v, and the verified inequality u/v>q/E_s gives `F_q=(E+u)q/E_s−sqrt(E²−M²)−v<0`. Both whole-band extremes are therefore evaluated at the **upper** q edge, without inferring continuous support from a finite angular grid.

| Evolved card | Whole-band threshold range E/T | q upper edge, cosine zero | Width from ω alone at cosine zero |
| --- | ---: | ---: | ---: |
| Fine, one charge subcell | 13.877–21.201 | 16.645–16.871 | .22649 |
| Fine, two charge subcells | 14.840–23.201 | 18.050–18.114 | .06390 |
| Third, one charge subcell | 14.299–22.034 | 17.255–17.395 | .13947 |
| Third, two charge subcells | 15.844–25.698 | 19.558–19.601 | .04339 |

At this fixed angle and q, the one- and two-subcell ranges do not overlap on either grid: the difference between complete evolved cards is larger than the intrabin optical uncertainty of either card. The four cards also differ in scalar mass and charged moments, so this comparison does not attribute the shift solely to one input. The entire retained band remains on the past spacelike transfer branch in each case. The parent threshold across the band is below the archived upper parent energy near 32, but the occupation near the threshold may be small and is not determined by kinematics. In particular these are not bounds on the rate, the Pauli-weighted response, or the physical photon population.

The [receipt](results/early-pole-endpoint-cold-threshold-enclosure-2026-11-28.json) checks the closed-form threshold against an independent scalar root and samples interior points as a guard on the analytical endpoints. As in the pole enclosure, this is a floating-point evaluation of analytical bounds, not a directed-rounding certification. From `kinetics/`:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_endpoint_cold_threshold_enclosure.py \
  FINE_X05.json FINE_D2_X05.json THIRD_X05.json THIRD_D2_X05.json
```
Use the decompressed x=.05 endpoint states from the [one-/two-subcell comparisons](INDEPENDENT-EARLY-POLE-ENDPOINT-CHARGE-AND-COLD-BAND-2026-11-28.md).
