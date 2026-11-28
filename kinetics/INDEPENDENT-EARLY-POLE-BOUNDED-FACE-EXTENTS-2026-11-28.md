# Advance the thin charged cells without clipping their occupations

The [padded stiffness check](INDEPENDENT-EARLY-POLE-PADDED-STIFFNESS-2026-11-28.md) found very fast charged responses in otherwise accurate frozen rates. An unrestricted implicit test overshot the Pauli bound even on a fixed chart. This independent pilot bounds the proposed *reaction extents*, not the occupations after a step.

Each signed event already supplies weights for a scalar, two daughter cells and two parent cells. For a proposed batch of events, sum candidate incoming and outgoing number separately on each destination cell. An outgoing portion is capped by available particles; a charged incoming portion is capped by available holes. The actual event multiplier is the smallest permitted ratio across all its deposits. Thus outgoing removal cannot exceed the inventory, incoming gain cannot exceed capacity, and the same multiplier preserves all linear event invariants. There is no post-step clipping. The two-evaluation fixed-chart variant averages the initial state and two successive bounded Euler advances, preserving the same convex occupation bounds.

For the frozen archive-started interval Δlog(a)=1e−6 with face fraction 1/4096, five subdivisions check fixed-chart conservation and time refinement. At 256→512 steps the scalar increment changes by 1.09e−4 (fine) or 1.49e−4 (third) for bounded Euler, and by 2.22e−8 or 5.36e−8 for the convex two-evaluation variant. The latter last-two refinement ratios are about 4.24 and 4.01. These ratios describe the observed scalar-number increment at this finite interval, not a general accuracy theorem for a limiter near active bounds. After each run the paired number and fixed-energy differences are checked independently.

The separate short evolving-chart pilot uses the same face grid, updates the Hartree mass and the radiation correction, and evaluates the fixed-clock collisions between two half increments of scale factor. An optional midpoint predictor evaluates the extents at the predicted state while applying capacities from the actual base state. The code measures energy changes by taking occupation and quasiparticle-energy differences before summing, independently integrating the collisionless trace on each half interval. At Δlog(a)=1e−6 the measured 64→128-step scalar increment changes are 4.37e−4 and 5.98e−4 for live bounded Euler, versus 4.50e−7 and 8.83e−7 for the midpoint predictor (fine and third). Midpoint's observed refinement ratios are 4.03 and 4.02. All sampled occupation bounds hold; the maximum paired number and continuity residuals are 3.52e−16 and 1.45e−15 in the stated normalization. These residuals compare independent budgets, not an imposed repair. Neither pilot replaces the default time driver or certifies a full early cosmological trajectory.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_bounded_face_extents_pilot.py
for g in fine third; do
  python -W error early_pole_live_bounded_face_extents_pilot.py "$g"
  python -W error early_pole_live_bounded_face_extents_pilot.py "$g" --midpoint
done
```
Compare the [held-chart receipt](results/early-pole-bounded-face-extents-2026-11-28.json) and [short evolving-chart receipt](results/early-pole-live-bounded-face-extents-2026-11-28.json).
