# The finite-mass logarithmic QED correction: equilibrium control

The [reciprocal gauge-energy proxy](early_pole_reciprocal_gauge_energy_functional_reference.py)
retains the factorized piece of the leading equilibrium interaction
while inserting finite bare energies. At nonzero mass there is also
an integrable logarithmic term. A small change in the integrated
pressure need not imply an equally small change in the trace or
in a pointwise fermion energy shift.

Scope: one unpolarized Dirac species, one massless Abelian gauge
boson, equal particle/antiparticle occupations, zero chemical potential,
common equilibrium temperature, and the order-g² vacuum-on-shell
coefficient. This does not assert that the project distributions are
thermal or that a logarithmic correction alone matches their
nonequilibrium self-energies.

* For equilibrium pressure use the factorized and logarithmic terms of
  [Bennett et al., Eq. (4.7)](https://arxiv.org/abs/1911.04504),
  equivalently [Eqs. (D.26) and (D.30) of Escudero et al.](https://arxiv.org/html/2511.04747v2).
  Integrate the symmetric log kernel on the triangle r<p with the
  substitution r/p=1−t⁴, resolving the diagonal without setting its
  value to zero. Use composite outer momentum panels (including the
  near-rest scales) and check both mesh and momentum cutoff.
* Independently reproduce the same pressure from one half the integral
  of the single-particle shifts, using the electron and photon pole
  formulae of [Grohs and Fuller, Eqs. (43)–(44)](https://arxiv.org/abs/1706.03391).
  For the pole integral the diagonal logarithm is integrated with
  logarithmic quadrature weights on each side of the pole. For the
  p→0 coefficient at fixed m>0 the limiting integral uses
  `log[(p+k)/|p−k|]/p → 2/k` after integration.
* Check the equilibrium thermodynamic trace with
  `(rho−3p)/(g²T⁴)=−x d[p/(g²T⁴)]/dx`, where x=m/T.
  This derivative includes the response of equilibrium occupations;
  it is not an expansion derivative at fixed comoving occupation.

The eight mass ratios from 0 to 5 recover the massless coefficient
−5/288, match pressure from the independently integrated poles within
8.0e−17 absolute, and match the thermodynamic trace derivative to
2.5e−14 absolute (all dimensionless). Here the logarithmic term is
always positive in pressure, reducing the magnitude of the negative
factorized contribution. At m/T=0.05 its size is 0.0183% of the
factorized pressure magnitude and 5.54% of the factorized trace
magnitude; in the p→0 coefficient it reduces the factorized fermion
mass-squared shift by 0.181%. At m/T=2 those pressure and pole ratios
rise to 4.89% and 21.7%; the logarithmic trace correction there has
the opposite sign from its value at m/T=0.05. These are coefficient
comparisons, not errors on the project’s kinetic trajectory.
In particular the p→0 comparison is not a soft-resummed
quasiparticle result.

Next test: retain the symmetric log kernel when defining a reciprocal
candidate for the two distinct charge-symmetric project populations,
then check its first variation and fixed-comoving trace separately.
A nonequilibrium derivation and matching to the supported real-time
collision operators remain separate requirements.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_log_thermodynamic_reference.py
```
[Raw receipt](results/finite-mass-qed-log-thermodynamic-2026-11-27.json).
