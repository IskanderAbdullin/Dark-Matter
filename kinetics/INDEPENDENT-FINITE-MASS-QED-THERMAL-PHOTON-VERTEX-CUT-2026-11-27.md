# The bare thermal photon cut of the scalar vertex is leading-soft only

In the one-loop scalar vertex, replace just the photon propagator
by its free null thermal cut. Keep both charged external pole spinors
at their vacuum masses. This selects one contribution in Feynman
gauge, not the complete thermal vertex or external residues.

For a photon direction n=(1,n-hat), the two signed insertions
k=±omega n have denominators −2P·k and −2b·k. Their product is
even. The on-shell projected scalar numerator is

    N(k) = 4(P·b − P·k − b·k + k²)
           + 2(M+m) (M b·k + m P·k)/(P·b+Mm).

The projected slash term is linear in k. For a null photon the
quadratic piece is zero, so symmetrizing positive and negative
thermal photon frequencies leaves precisely 8 P·b in the
numerator. At this perturbative order, this cut with the specified
external spinors has no finite remainder after its own leading
eikonal subtraction. The diagonal external-leg contributions
would still have to be treated in a complete rate; the separate
cross-term coefficient here is not the full radiator.

Explicit spin traces and azimuthal plasma-angle integration
confirm the result on both archived recoil geometries, including
forward directions. The all-angle residuals of the null test are
below 8.5e−18 in the scaled integral. For contrast a constructed
off-null four-vector (|k|/omega=.9) gives a nonzero remainder
scaling as omega² at fixed ratio; reducing the energy by a factor
of ten reduces the residual by a factor very close to 100.
That constructed ratio is not a photon spectral function.

This is an exact identity for the stated **bare photon cut**.
It does not remove the thermal charged cuts, residue corrections,
vertex matching in the resummed soft sector, or widths. It
cannot be substituted for a complete finite-temperature decay rate.

Replay in kinetics with the two archived spectra:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_thermal_photon_vertex_cut_reference.py
```
[Raw receipt](results/finite-mass-qed-thermal-photon-vertex-cut-2026-11-27.json).
