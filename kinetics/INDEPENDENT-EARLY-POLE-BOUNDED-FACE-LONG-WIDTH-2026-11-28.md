# Local time refinement does not remove the live thin-cell width dependence

The [bounded extent pilot](INDEPENDENT-EARLY-POLE-BOUNDED-FACE-EXTENTS-2026-11-28.md) advanced each archive with changing mass and expansion by Δlog(a)=1e−6. Now extend that interval fiftyfold with the same bounded midpoint method, and compare two thin face fractions independently of the time step: 1/4096 and 1/32768. The initial step histograms are identical before evolution within an archive, the same charged dimension is used for both widths, and the scalar and original charge edges are unchanged.

| Archive | Final Δ scalar number, 1/4096 | Final Δ scalar number, 1/32768 | Narrow/regular ratio |
| --- | ---: | ---: | ---: |
| Fine | 8.2218713e−7 | 3.0494640e−7 | .37090 |
| Third | 5.3600669e−7 | 9.3431421e−8 | .17431 |

Refine each separately through 128, 256 and 512 steps. The last relative time changes are below 8.16e−6 on all four runs. All tested charged cells respect Pauli bounds; total paired number and energy continuity close independently (maximum absolute energy residual below 3.24e−15). Nevertheless the increment falls strongly as the physical face fractions of the numerical grid shrink. The initially exact thin-face step reconstruction is therefore **not a converged live solver**. The time refinement only controls each particular discrete system. Investigate the narrow-cell kinetics or a different conservative representation before using it for a long history.

## Separate held-chart entropy control

With geometric occupation and hole factors, the change of Bose/Fermi entropy weighted by the collision source equals the sum of (forward−inverse) log(forward/inverse) over its signed events. Each term is nonnegative for interior occupations. The two evaluations on the archived inputs agree to 2.89e−15 relative. Over the short held chart interval Δlog(a)=1e−6 all recorded bounded Euler and convex two-evaluation substeps increase entropy at 16, 64 and 256 steps. This is a sampled finite-update check, not an unconditional discrete entropy theorem. It has no bearing on the observed live width discrepancy.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for g in fine third; do
  python -W error early_pole_live_bounded_face_extents_pilot.py "$g" --span 5e-5 --base 128 --midpoint
  python -W error early_pole_live_bounded_face_extents_pilot.py "$g" --span 5e-5 --base 128 --midpoint --pad .000030517578125
done
python -W error early_pole_bounded_face_entropy_pilot.py
```
See the [two-width receipt](results/early-pole-bounded-face-long-width-2026-11-28.json) and [entropy receipt](results/early-pole-bounded-face-entropy-2026-11-28.json).
