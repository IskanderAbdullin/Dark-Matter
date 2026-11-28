# Average an integrable charged-cut archive edge over external energy

The [closed radial result](INDEPENDENT-FINITE-MASS-QED-CHARGED-CLOSED-RADIAL-PV-2026-11-27.md)
found a pointwise logarithm when a positive bare charged loop coincides
with a nonzero step in the archived occupation. Instead of using the
value at that single external energy, integrate the same reference over
a finite symmetric external panel centered on the declared step. The
normalized parabolic weight is 3(1−u²)/4 on −1≤u≤1, where the
energy displacement is H u for panel half-width H.

If J is the coefficient of the radial cutoff logarithm at the edge,
the external-energy profile contains −J log|E−E_edge|. Its normalized
panel average is exactly J(4/3−log H) for a frozen J. Subtract the
logarithm before one quadrature and add the analytic expression.
For an independent calculation integrate the unmodified profile
on nonuniform subpanels, with no external sample on the edge.
Each sampled radial integral uses the full archived radial domain.

The fine daughter and third parent cases give finite averages on all
four panels. The independent unsubtracted external quadrature agrees
within 1.4e−6 in the raw coefficients. Reducing the panel half-width
by four raises the step average by +1.001480 and +0.853419,
respectively. The singular term alone predicts +1.001878 and
+0.853900. The remaining smooth-background response is below 4.9e−4.
Under the same width change the simple continuous interpolation
control changes by less than 5e−5. The finite step vs continuous
difference at a fixed width is much larger than the quadrature
error, and the continuous interpolation is not moment preserving.

The isolated pointwise archive-edge divergence is therefore
integrable in external energy. An outer panel and its width are part
of the stated mathematical reference, not a physical regulator or a
calibrated rate. The test leaves thermal cut signs, pole shifts,
and matched medium corrections unspecified.

Replay from kinetics; this performs many full radial solves:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_charged_outer_edge_projection_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-charged-outer-edge-projection-2026-11-27.json).
