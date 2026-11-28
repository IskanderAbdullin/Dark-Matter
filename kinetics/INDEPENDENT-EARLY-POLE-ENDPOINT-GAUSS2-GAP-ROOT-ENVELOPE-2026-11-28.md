# A fixed-reference moment envelope bridges to a frozen Hartree root bound

Prescribe every shell's count N and mean energy at a **reference mass** m₀. The [cellwise Jensen/chord bound](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-ENERGY-MOMENT-ENVELOPE-2026-11-28.md) supplies a total inverse-energy interval [L,U] at m₀ for any nonnegative occupation with these moments and supports. Now freeze each candidate shape as a measure in q. Its inverse-energy integral I(m) decreases with m and the Hartree residual `m²−m_bare²−(λ/2)I(m)` is strictly increasing for λ>0. If the new root lies above m₀, I there is at most I(m₀)≤U. If it lies below m₀, I there is at least I(m₀)≥L. Thus all candidate roots are bounded by
```text
lower = min(m₀, sqrt(m_bare² + λ L/2))
upper = max(m₀, sqrt(m_bare² + λ U/2)).
```

The moment constraints are imposed once at the reference mass; the proof does not assume unchanged shell energies as m moves. Compare two numerical fixtures: an early genuine warm-Gauss checkpoint and a late S8 stored histogram lifted onto the positive volume nodes at fixed shell count. In each case the reference is its Gauss nodal root. Fit the independent [positive energy-matching reconstruction](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-GAUSS2-ENERGY-RECONSTRUCTION-2026-11-28.md), freeze the reconstructed q distribution, and solve its own mass root.

| Example | Reference nodal mass | Relative width of all frozen roots | Relative shift of the chosen positive representative |
| --- | ---: | ---: | ---: |
| Early, x=.005002 | 2.7384042e−10 | .016108 | +.00342594 |
| Late lifted S8, x=.05 | 8.7318281e−4 | 6.873e−10 | +1.138e−12 |

These are analytic bounds evaluated in floating point; the independent root solve uses a relative numerical margin of 5e−12. No late full Gauss-node history was run here. The early smoke state does not predict a physical onset clock. The strikingly different widths apply only to the prescribed moments of these two numerical states. They do not bound a dynamically varying distribution, a rate aperture, or cosmological observables.

Run `early_pole_endpoint_gauss2_gap_root_envelope.py` with the exact Gauss energy-reconstruction smoke state and then the S8 x=.05 state. The [receipt](results/early-pole-endpoint-gauss2-gap-root-envelope-2026-11-28.json) contains both inverse-energy and root intervals.
