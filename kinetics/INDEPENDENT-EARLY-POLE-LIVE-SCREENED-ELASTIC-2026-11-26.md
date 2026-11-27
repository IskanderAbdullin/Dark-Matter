# Live screened elastic packets: stiffness and an unresolved cold-gain sign

The [frozen direct-Q calculation](INDEPENDENT-EARLY-POLE-SCREENED-ELASTIC-TRANSFER-2026-11-26.md)
finds heat into the daughter on both archived spectra. Here repeat
with live Fermi count and centered energy moments, the same scalar
Hartree mode and optional photon pair channels from the previous
[matched cold-window test](INDEPENDENT-EARLY-POLE-ARCHIVED-SHARED-COLD-WINDOW-2026-11-26.md).
The elastic packets use the frozen archived polarization below
`Q=.08` and the vacuum kernel above. They replace rather than add to
the previous t-hard elastic component. The phase-space invariant
cut remains `s≥1`, the common scalar window is fixed, and the live
fermion radial window is [.45, 5]. The comparison is at fixed scale
and uses the reference clock, not a late cosmological solution.

Each packet carries one parent and one daughter in both directions.
Remap the final daughter radial coordinate to close the common
variational energy. Check that the same momentum transfer magnitude
can connect the initial and final daughter on that packet; monitor
how remapping affects the originally integrated invariant s support.
At the initial archive the transfer-squared-weighted difference in
its allowed angular fraction is at most `6.1e-5`. Across the runs
the count/energy and entropy ledgers close at numerical precision.
The guaranteed nonnegative local entropy contribution is written
using the forward/reverse affinity to avoid cancellation at tiny Q.
With exact common-temperature affinities, the packet residual is near
roundoff even with different chemical potentials, and a hotter parent
transfers energy to the daughter.

The finer packet quadrature gives initial daughter heat 273.20
(fine) and 261.95 (third). Against the independently integrated
number-preserving frozen states in the identical [.45, 5] box
(278.84 and 274.29), the live static comparison is within 2.1% and
4.6%, respectively. Differences include the momentum reconstruction,
the massless weight and dressed energy remap. The split is about
half screened and half vacuum. This agreement does not imply
agreement of the response of the much narrower scalar aperture.

## A stiffness and sign check

In a separately manufactured thermal state, freeze the remap and
linearize the reciprocal packets. The positive collision Gram matrix
annihilates both fermion counts and total variational energy. Its
susceptibility-normalized spectral radius is `6.28e8` (fine) or
`9.94e8` (third) per reference log-scale interval at the finer
packet quadrature. Using only the same graph's vacuum-transfer
part gives around `1.5e6`. These are numerical relaxation scales
for this projection; they are not physical full-plasma rates. An
unrestricted explicit trial step left the admissible moment domain.
Starting explicitly at `1e-11` and allowing adaptation resolves only
the early transient shown here.

At exposure `8e-8`, with matched initial states and photon-pair
dynamics enabled:

| Archive | Scalar gain without elastic | With screened + vacuum elastic | Relative change |
| --- | ---: | ---: | ---: |
| Fine | 6.29061e-9 | 5.25234e-9 | -16.51% |
| Third | 5.80557e-9 | 6.55010e-9 | +12.82% |

Without photon pairs the percentages are -16.50% and +12.82%.
A separate strict-hard-only fine control at the same exposure gives
a change of just +0.0029%, not the previous longer exposure result.
The soft portion cannot be added to that run without accounting
for the altered packet support and common conserved terms.
All sampled growth increments still have the same positive sign:
the relative effect of elasticity differs, not the sign of overall
scalar production. The stiff runs take 520 or 412 RHS evaluations
at this duration, and absolute energy drift is below `5e-16`;
the integrated entropy mismatch is below `4e-12`. The two archive
shapes are also mapped to different native-aligned hard grids, so
a physical sign cannot be inferred. The next grid control must use
a common hard mesh and a controlled representation of the initial
occupations rather than comparing the two native reconstructions
as if they differed only in scalar-bin size.

No final scalar freeze-out, transfer function or complete transport
rate follows. Screening is held frozen; no other gauge (2↔2)
channels or expansion are included in this test.

Run from `kinetics/` with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`:

```sh
python early_pole_frozen_elastic_direct_q_reference.py --check-azimuth
python early_pole_frozen_elastic_direct_q_reference.py --grid fine --form cell-conservative --level 1 --window 5e-5 10 --bounds .45 5
python early_pole_archived_cold_screened_elastic_reference.py --grid fine --level 1
python early_pole_cold_screened_elastic_stiffness_reference.py third --level 1
python early_pole_archived_cold_screened_elastic_reference.py --grid fine --evolve
python early_pole_archived_cold_screened_elastic_reference.py --grid fine --evolve --strength 0
```

Raw: `results/early-pole-live-screened-elastic-2026-11-26.json.gz`.
