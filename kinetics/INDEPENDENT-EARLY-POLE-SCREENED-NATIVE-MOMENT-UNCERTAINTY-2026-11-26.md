# Native-shell shape uncertainty with count, energy and gauge moment fixed

The [continuous-shell control](INDEPENDENT-EARLY-POLE-SCREENED-HEAT-VS-REDISTRIBUTION-2026-11-26.md)
preserved native cell number but not its energy. This test asks a
narrower question: does the instantaneous screened scalar response
change materially when even the archived number, bare energy and
reciprocal gauge moment are held fixed in each native daughter
shell?

For a native shell with constant archived occupation `f0`, add a
bounded perturbation proportional to
`epsilon f0 (1-f0) (1-t²) P(t)` where `t` spans the shell from
−1 to +1. Subtract from a cubic, quartic or quintic probe its
three projections with respect to the weighted measure
`p² (1-t²) dp`. The constraints span `1, E, 1/E` for the exact
daughter bare mass. Computing them in centered energy coordinates
avoids a nearly singular difference between reciprocal energy and
its constant and linear terms. The result gives zero differences
in count, energy, and the field's `J ∝ ∫ p² f/E dp`, independently
in each altered source shell. The perturbation vanishes at native
edges and respects Pauli bounds. Only daughter shells intersecting
[1.39, 1.48] are changed; the parent, exterior field and scalar
archive are untouched.

An independent quadrature audits the three moment differences.
Use one common third hard aperture, with four subdivisions only
in cells intersecting the narrow production region. This is
needed because a two-moment cell aligned exactly with a full
native shell is blind to the perturbation's first two moments.
Both signs of an amplitude .75 are tested per polynomial. These
are deliberately chosen uncertainty witnesses, not inferred
probability bounds or reconstructed true distributions.

The independent per-shell moment checks close to better than
`5.4e-14` relative to the redistributed moment. The maximum
change in occupation across these probes is about `0.067`.
The total cold-source tangent in units of `10^6` per reference
exposure squared is:

| Probe | Amplitude | Fine archive | Third archive |
| --- | ---: | ---: | ---: |
| Unchanged | 0 | -1.963 | +1.821 |
| 3 | -0.75 | -0.043 | +0.408 |
| 3 | +0.75 | -3.858 | +3.326 |
| 4 | -0.75 | -0.778 | +3.856 |
| 4 | +0.75 | -3.151 | +0.007 |
| 5 | -0.75 | -0.036 | +0.936 |
| 5 | +0.75 | -3.844 | +2.981 |

These variants can nearly erase the local response or amplify it
substantially without changing its sign unambiguously within this
two-sided family. The heat-only part remains negative for all
the variations. The sensitivity therefore cannot be attributed
only to the earlier continuous reconstruction's changed native
energies. The active mesh still projects to two moments per
new cell, so this is a demonstrated ambiguity of the permitted
native inputs in the finite model, not a rigorous bound over
every allowed sub-shell function.

The comparison is local in time: it tests the derivative of
the cold scalar source along the screened elastic direction.
It is not an evolved spectrum, a full-transport sign check, or
a claim that any of these within-shell profiles is realized.

Reproduce from `kinetics/` with one BLAS thread. Run these sequentially
to avoid retaining multiple refined packet graphs at once.

```sh
python early_pole_native_shell_shape_invariants_reference.py --degree 3
python early_pole_cold_screened_affinity_reference.py --mesh third --subdivision 4
python early_pole_cold_screened_affinity_reference.py --mesh third --subdivision 4 --shape-degree 3 --shape-amplitude -.75
python early_pole_cold_screened_affinity_reference.py --mesh third --subdivision 4 --shape-degree 3 --shape-amplitude .75
```

Raw: `results/early-pole-screened-native-moment-uncertainty-2026-11-26.json.gz`.
