# Shell counts and mean energies bound pressure without choosing a within-cell shape

A [positive scalar quadrature](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-POSITIVE-CELL-QUADRATURE-2026-11-28.md) reproduces the moments of one assumed constant-in-shell histogram. Here allow **arbitrary nonnegative occupation measure** within each saved scalar shell. Assume the total count N and, in the stronger test, the mean energy Ē in every shell are prescribed at the same fixed mass m. The kinematic energy support is [A,B], with A=sqrt(a²+m²) and B=sqrt(b²+m²) for the radial edges a,b.

For pressure numerator P(E)=E−m²/E and v²(E)=1−m²/E², both functions increase and are concave. Their shell expectations lie between the endpoint chord at Ē and the function value at Ē. Inverse energy 1/E is decreasing and convex, so its expectations have the opposite Jensen/chord order. For inverse energy, the exact chord surplus is
```text
Chord(1/E)(Ē) − 1/Ē = (Ē−A)(B−Ē)/(A B Ē) .
```
The pressure width is m² times this surplus. Evaluate these differences in the stable excess energy r=E−m=q²/(E+m), and sum the intervals weighted by the prescribed shell counts. Count-only bounds use just the extreme endpoint weights.

The four matched x=.05 charts use the common S2 mass. In the example controls the per-shell mean energy is computed from the archived constant-in-shell reference; the resulting bound covers all shapes with those same counts and energy means, not just that example.

| Chart | Pressure width with count only | Pressure width with count and energy | v² width with count and energy | Inverse-energy width with count and energy |
| --- | ---: | ---: | ---: | ---: |
| Native | 35.696% | 1.111e−7 | 1.681e−5 | 1.389e−9 |
| S2 | 35.643% | 1.103e−7 | 1.705e−5 | 1.377e−9 |
| S4 | 35.617% | 1.101e−7 | 1.700e−5 | 1.375e−9 |
| S8 | 35.605% | 1.101e−7 | 1.698e−5 | 1.375e−9 |

Widths are relative to the corresponding total moment of each example, not to particle count. Holding each mean energy tightly constrains these three smooth moments, while allowing unresolved spatial spectra. This does not provide the actual mean energies or their time evolution, control a moving mass gap globally, establish a physical tail, or constrain apertures and CLASS transfer. The inequalities are analytic; numeric endpoints are evaluated in floating point with explicit roundoff tolerances, not with directed interval rounding.

The [receipt](results/early-pole-endpoint-scalar-energy-moment-envelope-2026-11-28.json) checks a distinct positive three-atom measure and the cellwise constant reference at two quadrature orders. From `kinetics/`, run `early_pole_endpoint_scalar_energy_moment_envelope.py` with the matched native 4096-step-start, S2, S4 and S8 x=.05 state archives.
