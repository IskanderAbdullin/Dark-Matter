# Crossed hard meshes and the local heat/redistribution split

The [live screened run](INDEPENDENT-EARLY-POLE-LIVE-SCREENED-ELASTIC-2026-11-26.md)
gave opposite short scalar-gain shifts on its native archive grids.
The two archived spectra can now be projected independently onto
either of the same two hard aperture meshes. The charged count and
centered energy in each new cell are obtained by integrating the
archived occupations. The scalar bin support, reference exposure,
elastic cut, fixed polarization and optional photon-pair channel are
held to the conventions of the preceding comparison; each source
archive retains its own scalar state and polarization.

With photon pairs enabled, at exposure `8e-8` the relative gain
shift from turning on this elastic subset is:

| Archived state | Fine hard mesh | Third hard mesh |
| --- | ---: | ---: |
| Fine | −16.51% | −11.17% |
| Third | −4.37% | +12.82% |

In this experiment the third archive changes sign when only its
hard-cell representation changes. Both archives suppress the
gain when placed on the fine hard mesh; on the third hard mesh
their relative signs differ. All the gains are positive. Without
photon pairs, the same comparison gives −16.50%, −11.17%, −4.38%,
and +12.82%. Using the same hard mesh is necessary but not a
physical convergence criterion.

## Which part of the elastic direction drives the initial response?

At a fixed archived state write each signed reciprocal packet
rate as `W Λ A`, where `A` is the dimensionless FD affinity and
`Λ` is the positive logarithmic mean of forward and inverse
occupation factors. Let `s` be the signed vector of cell moment
deposits for a forward packet, and `e` its variational daughter
energy change. Define:

- `G = Σ W Λ e²`, the instantaneous heat-direction mobility;
- `Qdot = Σ W Λ e A` and `Sdot = Σ W Λ A²`;
- `F_heat = (Qdot/G) Σ W Λ e s`, and `F_shape = F_total − F_heat`.

Their sum is the actual elastic direction. Each conserves both
fermion counts and the total variational energy; the shape
remainder changes daughter energy by zero to roundoff. The split
is orthogonal in the instantaneous packet mobility, and
`Qdot²/(G Sdot)` is the fraction of entropy production aligned
with daughter heat. A synthetic exact Fermi state with different
inverse temperatures and chemical potentials gives fraction one
and `Qdot/G = beta_d − beta_P` on the same graph. This is a local
diagnostic, not a separate claim that either projected vector by
itself generates globally positive occupations at finite time.

At finer packet resolution, for the step-occupation projections
only 0.15–0.48% of the archive entropy production aligns with
daughter heat. More decisively, the local derivative of the cold
scalar source along that heat direction is small and negative
in all four crossed runs. The zero-net-heat remainder drives the
larger derivative and owns the sign change. Its daughter-cell
component accounts for most of the effect. Cell attribution finds that a small group of daughter cells near
1.41–1.47 dominates. In the third/third case the two daughter
cells [1.4315, 1.44] and [1.44, 1.4667] contribute about +5.52 and
−1.78 in the table units, so an unresolved sign in this aperture
is not surprising. The derivative is
checked independently with forward moment displacements and
central displacements of admissible Fermi logits rather than
evaluating an inadmissible backward moment step.

The following entries are derivatives of the same instantaneous
scalar source with respect to reference exposure (units of
`10^6` for the tabulated number), not finite-time changes:

| Archive / hard mesh | Full elastic direction | Daughter-heat direction | Zero-net-heat remainder |
| --- | ---: | ---: | ---: |
| Fine / Fine | −2.752 | −.009 | −2.743 |
| Fine / Third | −1.033 | −.013 | −1.020 |
| Third / Fine | −.125 | −.020 | −.105 |
| Third / Third | +3.622 | −.034 | +3.656 |

## Representation control

Within the same original archive, replace the step distribution
by the derivative of a monotone cumulative count in `p³/3`.
This preserves each full native shell's particle number but moves
occupations within shells; the exterior is held fixed in the
short live model. It is a controlled alternative initial
representation, not an assertion that the actual distribution
is smooth. On the unchanged native apertures this changes the
initial fine tangent from −2.752 to +.200, and third from +3.622
to −1.360. In all four smoothed comparisons, daughter heating
remains positive, the heat-only scalar tangent remains negative,
and the shape remainder dominates. Consequently neither the
sign of heat exchange nor agreement of integrated counts determines
the evolution of the sharply supported scalar source.

At the shorter fine exposure `2e-8`, increasing the packet
quadrature changes the integrated step-state relative result from
−12.78% to −12.68%. This paired check limits *packet* quadrature
error for the short fine control, not the hard-mesh or occupation
representation uncertainty isolated above. A separate local mesh
sweep on the shared third aperture with the smooth representation
subdivides cells intersecting [1.395, 1.475]. Going from two to
four subdivisions keeps both tangents negative, from −.308 to
−.261 (fine archive) and −1.174 to −1.213 (third archive) in the
table units. This does not establish the physical sign or full
grid convergence.

To reproduce from `kinetics/` (one BLAS thread):

```sh
python early_pole_cold_screened_affinity_reference.py --level 1 --cells
python early_pole_cold_screened_affinity_reference.py --level 1 --projection cell-conservative
python early_pole_cold_screened_affinity_reference.py --grid fine --mesh third --projection cell-conservative --subdivision 4
python early_pole_archived_cold_screened_elastic_reference.py --grid third --aperture-grid fine --evolve
python early_pole_archived_cold_screened_elastic_reference.py --grid third --aperture-grid fine --evolve --strength 0
```

Raw: `results/early-pole-screened-heat-vs-redistribution-2026-11-26.json.gz`.
Individual cell derivatives: `results/early-pole-screened-heat-cell-attribution-2026-11-26.json.gz`.
No finite-temperature conclusion about the final scalar relic or
its transfer function follows.
