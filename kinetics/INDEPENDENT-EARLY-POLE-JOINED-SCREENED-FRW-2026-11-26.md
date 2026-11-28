# Joined screened channels in comoving expansion: energy and clock checks

This reference puts the five finite tested channels
(photon pairs, hard Compton, hard flavor conversion, unlike
screened elastic and same-flavor screened pole) into the
[comoving cold-window model](INDEPENDENT-EARLY-POLE-ARCHIVED-COLD-FRW-SPECTATOR-2026-11-26.md).
The cold scalar and its exterior Hartree contribution are
retained. The collision supports and vacuum matrix weights
are fixed at their initial comoving coordinates; the
screened responses are held at the initial archived state
in this expansion test. It is not a full cosmological
collision or plasma calculation.

In the two slow half steps, physical expansion advances
and the independent trace is integrated. The elastic
implicit step stays at the midpoint scale factor and uses
both its full variational gradient for remapping and the
exact **finite** difference of comoving energy. This
difference is computed by integrating changes of the
occupation functions against redshifted energies and by
taking the difference of the quadratic gauge term as a
bilinear form. The fixed exterior matters for the gauge
cross term and absolute Hubble energy, but cancels in the
finite occupation increment at a fixed scale factor.

The two energy tests are independent: a fourth-order
finite difference in scale at fixed occupation matches
the trace to better than `1.9e-11` on both archived grids,
and the exact finite occupation-increment formula agrees
with the difference of absolute energies to better than
`7e-16` across the tested fixtures. The conservation
ledger in the coupled calculation subtracts the integral
of the comoving trace; it is not a test of constant
comoving energy in the presence of massive redshift.

The ordinary reference clock uses visible radiation.
Massless dimensionless (2↔2) rates per log scale
receive the expected relative factor `exp(eta)`; the
cold decay uses its existing redshifted shell kinematics.
The second clock independently includes the conditional
dark energy in (H), sampled at the actual slow stages and
pre-elastic midpoint. This extra clock is a controlled
numerical comparison, not a fitted cosmology.

For a fine-grid reference exposure `2e-6`:

| Setting | Cold scalar gain |
| --- | ---: |
| Fixed scale, radiation clock | 1.5212995e-7 |
| Expansion, radiation clock | 1.5213271e-7 |
| Expansion, conditional full H | 1.4808550e-7 |

The short reference expansion changes the gain by about
+0.0018%; the optional full H clock lowers it by about
2.66% relative to the expanding radiation-clock result.
The comoving trace integral in the expanding runs is
about `7.53e-10`; the redshift-subtracted energy defect
is less than `9e-14`, with the entropy mismatch less
than `2e-12`. All comparisons are confined to early
reference intervals and to the declared active/exterior
split. They neither establish freeze-out nor yield an
observable structure transfer function.

Reproduce from `kinetics/` with one BLAS thread:

```sh
python early_pole_joined_screened_frw_reference.py fine --fixture
python early_pole_joined_screened_frw_reference.py fine
python early_pole_joined_screened_frw_reference.py fine --expands
python early_pole_joined_screened_frw_reference.py fine --expands --full-H
```

Raw: `results/early-pole-joined-screened-frw-2026-11-26.json.gz`.
