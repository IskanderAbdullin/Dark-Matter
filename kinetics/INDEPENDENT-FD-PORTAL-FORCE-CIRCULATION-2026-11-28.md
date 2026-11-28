# A full FD differential flow has nonzero force circulation away from equilibrium

Hold the visible T fixed and retain the same timelike charged-lepton
phase space and affinity as in the differential reciprocity test.
Integrate (F,Q/T_SM) counterclockwise around three force rectangles,
using the original independently evaluated rates for the full-FD boundary.
Separately integrate the analytic interior curl with the local blocking
derivative. Both positive and negative circulations are resolved; the
pure hot-force and cold-force rectangles have opposite signs as predicted
by the signed affinity.

Replace only the dark blocking factor with a convex blend of one and
the original dark blocking value. Keep visible blocking and the affinity
unchanged. At fixed state the interior differential curl is proportional
to the retained dark blocking fraction. In all three rectangles a
half fraction gives exactly half the full result within float accuracy;
at zero fraction the interior result vanishes and the boundary residual
is consistent with quadrature. This is a controlled mechanism toggle,
not a trajectory comparison between ideal gases.

Full-FD boundary and interior agree within 1.51e−13 of the absolute
boundary scale at the reported mesh; the relative change of the resolved
curl between the two independent orders is below 5.2e−8. Differential
reciprocity at a selected point is not a global potential for these
force variables. This does not contradict detailed balance at the
common equilibrium state or nonnegative entropy at each sampled state.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_force_circulation_reference.py
```
