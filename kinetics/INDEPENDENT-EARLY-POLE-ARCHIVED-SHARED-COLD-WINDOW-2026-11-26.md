# A matched cold window separates scalar binning from the hard-grid response

An original scalar bin is not interchangeable between the two
archives: the bin containing q=4.1e−6 and its initial occupation
differ considerably. Instead transfer the original counts by exact
shell-volume overlap onto the same interval 1e−6≤q≤1e−5 on both
archives. Evolutions use independent scalar bin counts and one
Hartree gap, keeping the exterior scalar counts as true energy
spectators. Both rates use the chosen vacuum reference
normalizations on the represented domain. This is still a
fixed-scale finite-box comparison with only the hard selected pair
subset and no matched medium vertex.

In the narrow active hard box, using four versus eight common scalar
bins changes the short scalar gain by at most .011% relative for
both archives and both channel switches. That does not establish
hard-cell convergence. The baseline initial slope from the raw
fermion bins provides a sharper test: independent hard-energy
aperture integration at each scalar mode differs from the coarse
local Fermi reconstruction by 4.68% on the fine archive and 21.14%
on the third archive in normalized L1. Splitting the live hard
cells at actual archived edges near the resonance reduces these
initial discrepancies below 1e−5 on both. The independent
raw reference uses the first-order-equivalent square-root pole
model, not the reconstructed rates.

| Hard support and near-aperture projection | Suppression of short scalar gain, fine | Third |
| --- | ---: | ---: |
| Narrow, coarse near-aperture cells | 3.77% | 3.65% |
| Wider upper cut 4/5, coarse near-aperture cells | 7.94% | 7.69% |
| Narrow, actual archive splits near aperture | 1.57% | 1.78% |
| Wider upper cut 4/5, actual splits near aperture | 3.51% | 4.18% |

Each uses the same archived initial state as its matched decay-only
control over reference exposure 8e−5. With the refined aperture,
the fine incident quadrature check gives 1.62% for the narrow
case. Extending the lower charged support from .45 to .15 while
keeping the total initial hard field domain identical changes the
refined fine wider result from 3.51% to 3.73%. The sign and direction
of the upper support effect are shared here, but its magnitude
depends materially on how the narrow fermion aperture is represented.
No scalar-momentum refinement by itself fixes that error.

The combined energy, count and entropy identities remain checked.
The coarse and refined runs do not establish a common physical
scalar yield, longer evolution or a CLASS-ready spectrum. They
show which numerical representation was masking a narrow transfer.

Run from kinetics/ with
early_pole_archived_shared_cold_window_reference.py --grid fine or third
(--n 4/8, --wide 0/1, --refine-aperture; use --floor .15 and
--field-floor .15 for the low-side comparison).
The independent initial-slope comparison is
early_pole_archived_cold_initial_flux_reference.py.
Raw: results/early-pole-archived-shared-cold-window-2026-11-26.json.gz.
