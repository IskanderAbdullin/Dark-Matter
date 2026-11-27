# A reciprocal finite-mass log kernel on hard momentum shells

The [equilibrium control](INDEPENDENT-FINITE-MASS-QED-LOG-THERMODYNAMIC-2026-11-27.md)
isolates the nonfactorizable one-photon exchange term. Here we project
its charge-even, bare-on-shell form against two distinct fixed isotropic
fermion shapes and a shared photon shape on the hard domain q in
[0.25,8]. The chosen shapes and masses are controlled inputs, not an
evolved cosmological plasma. The factorized part is the same as in the
[earlier reciprocal proxy](early_pole_reciprocal_gauge_energy_functional_reference.py).

For two on-shell fermions of the same mass, let
`D=E_p E_k−pk cos(theta)`. The two equal-occupancy signed-energy
branches have spin-trace divided by exchange denominator proportional
to `(−8s+16m²)/(2(m²−s))`, with s equal to D and −D. With average
normalized to one at zero mass, their angular integral is
`1−m² log[(p+k)/|p−k|]/(2pk)`. One way to simplify the endpoint
ratio is to set p=m sinh(u), k=m sinh(v); the endpoint invariants are
m² cosh(u±v), and the resulting ratio of half-angle tangents is
(p+k)/(p−k). Twelve direct 65-digit angular quadratures, including
nearly coincident unequal momenta, match the reduction within 5.01e−51.
The coefficient also reduces to the equilibrium fermion pole term of
[Grohs and Fuller, Eq. (43)](https://arxiv.org/abs/1706.03391).

Use a common symmetric shell matrix for the energy and for its
derivative, so the two frequency branches do not create mismatched
interaction energies. The log kernel is weakly singular on the
same-shell diagonal; integrate in two triangles, using a fourth-power
distance substitution. Adjacent shells have a common-corner
singularity handled by splitting the normalized rectangle and factoring
out its small radial scale. No log value is discarded on the diagonal.
Use fixed comoving q and hold f fixed. For scaled mass x=a*m, the
log-kernel scale derivative multiplies its integrand by
`2−x²/(p²+x²)−x²/(k²+x²)`. Thus the trace of the interaction term
follows by differentiating a⁴rho, which is distinct from differentiating
the equilibrium pressure with temperature-dependent occupations.

For the reference small-mass shapes (0.05,0.003), at 64 shells the
log part is 14.8% of this interaction trace, while it changes the
largest individual fermion shift by only 0.131% relative to the
factorized value. The larger-mass stress case (0.5,0.25) reaches
about 9.69% for the largest shift. In both cases the log energy for
piecewise-constant shapes is within 0.45% of the smooth hard-shape
reference at 64 shells and approaches second order under refinement.
The fixed-shell 24/40 point matrix difference is below 2.91e−11
relative to the log-matrix scale, so the shell-compression error and
diagonal quadrature error are distinct. The shared midpoint first
variation closes the tested finite energy increments within 5e−23
and the fixed-comoving trace finite difference agrees within 9.1e−17
in absolute project test units.

This is only the real bare-on-shell hard exchange candidate and its
algebraic energy/trace projection. It is not a proof of quasiparticle
kinetics, of dynamical screening, of a matched finite-mass collision
operator or of a physical decay aperture. The next check is to drive
this same energy gradient on the archived hard shell shapes, without
calling the smooth test an archive measurement.

Replay:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_shell_log_reference.py
```
[Raw receipt](results/finite-mass-qed-shell-log-2026-11-27.json).
