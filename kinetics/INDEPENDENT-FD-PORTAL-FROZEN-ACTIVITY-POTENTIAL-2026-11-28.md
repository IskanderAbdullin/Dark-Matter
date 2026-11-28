# A convex test-force potential does not make the physical force field conservative

At a fixed local bath state, let the forward and reverse event prefactors
be a+ and a−, with ratio exp(d), and let gamma be their geometric mean.
For the mathematical test force v define the integrated potential

    Phi_state(v) = 4 integral gamma [ cosh((v_alpha + z*v_beta)/2) - 1 ].

At v equal to the actual force, its gradient equals the original finite
number and heat flux vector. Its test-force Hessian is the folded
positive reaction matrix with weights (a+ + a−)/2. The factor gamma is
held fixed only for these test derivatives. Cosh-type dissipation
potentials for jump processes are discussed in [Liero, Mielke, Peletier
and Renger, *On microscopic origins of generalized gradient structures*](https://arxiv.org/abs/1507.06322).
The local relation here is derived independently from the two directional
rates, rather than assuming the complete cosmological evolution has
that mathematical form.

The physical derivative permits the activity to vary with the bath.
The difference between that derivative and the fixed-activity Hessian
is independently reproduced by integrating the analytic derivative of
gamma. This term carries the measured antisymmetry. At equilibrium
it vanishes. At fixed T_SM=.05 GeV and the sampled cold-force pair
(−.6,−.16), its largest component is about 78.1% of the largest physical
Jacobian component. The eight controls give a scaled split defect below
3.0e−15 and an independent test-gradient/Hessian finite-difference
defect below 1.59e−13 on the reported mesh. Both quadratures reproduce
the same convex Hessian.

The test-force construction is local to the stated finite timelike
lepton portal. It neither contradicts the physical force circulation
nor proves a trajectory, expansion stability, or physical abundance.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_frozen_activity_potential_reference.py
```
