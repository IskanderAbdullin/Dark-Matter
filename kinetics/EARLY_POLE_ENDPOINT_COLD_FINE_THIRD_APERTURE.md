# The third-grid cold branch opens at higher incoming parent energy

Use the actual fitted fine/third D2 poles from the complete held cold-band result. Both models share the native reaction masses for this diagnostic. At fixed scalar momentum `q`, varying the scalar gap `m_s` and photon energy `omega` changes transfer energy `t = sqrt(q² + m_s²) − omega < 0`. Varying the angle changes the transfer magnitude `v`. For the past-directed spacelike branch, set `A = M² + t² − v² − m² > 0`, with `M` and `m` the parent and partner masses.

The physical opening solves
```
F(E) = v sqrt(E^2-M^2) + E t - A/2 = 0.
```
Independently verify this equation at the stored root. At the root, the derivatives of F with respect to E, t and v are strictly positive: the first exceeds v+t, the second is E−t, and the last is sqrt(E²−M²)+v. Thus the minimum admissible parent energy decreases with increasing t or v. For fixed q it increases monotonically over cosine from −1 to +1. The third chart has the more negative t at the sampled q, so its threshold is higher at **every** cosine there. The two endpoint values bound the complete angular range.

At cosine zero, the exact threshold changes as follows, in the numerical parent energy-over-temperature coordinate:

| Scalar momentum | Fine (E_{min}) | Third (E_{min}) | Symmetric scalar-gap part of shift | Symmetric photon-energy part |
| ---: | ---: | ---: | ---: | ---: |
| 2.5e−6 | 18.12431 | 19.63396 | +1.22049 | +0.28916 |
| 3.5e−6 | 18.12287 | 19.63227 | +1.22029 | +0.28911 |
| 4.7e−6 | 18.12052 | 19.62950 | +1.21995 | +0.28903 |
| 8.2e−6 | 18.10975 | 19.61682 | +1.21840 | +0.28867 |

The arrows are a kinematic explanation for sensitivity near this aperture. They are not a substitution for a rate calculation: residue, charge occupations, angular measure and matrix elements still enter, and no assumed thermal multiplier is applied here.

```sh
cd kinetics
python -W error -m unittest test_early_pole_endpoint_cold_native_aperture_rectangle
python -W error early_pole_endpoint_cold_native_aperture_rectangle.py \
 results/early-pole-cold-fine-third-d2-native-cdf-2026-11-28.json
```
