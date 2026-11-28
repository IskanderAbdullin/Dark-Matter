# Differential reciprocity and number/heat zeros in a finite FD two-bath displacement

Use the existing independent finite-width charged-lepton vacuum benchmark,
with visible temperature held fixed. Write alpha = 2 eta_D,
beta = 1 − T_SM/T_D, and z = E_lab/T_SM. The number/heating flow vector
is (F,Q/T_SM). Entropy production remains alpha*F + beta*Q/T_SM >= 0.

The local dark blocking factor is B_D=(1−f1)(1−f2). At each physical
event, the difference between its beta derivative and z times its
alpha derivative is

    (e1 - e2) * (f2 - f1) * B_D / 2,

where e1 and e2 are the two dark energies divided by T_SM. It is
nonnegative because the FD occupation decreases with energy. The
differential antisymmetry in the integrated off-diagonal elements
weights this factor by exp(alpha+beta*z)−1, along with the common
positive kernel and visible weight. If the affinity has one nonzero
sign throughout support, the antisymmetry has the same sign.

At the common thermal state the Jacobian is exactly the symmetric
equilibrium moment matrix. Away from that state the blocking
derivative need not be symmetric; this does not conflict with local
detailed balance or nonnegative entropy production. At T_SM=.05 GeV
the off-diagonal difference divided by their geometric mean is
+1.08585% for (alpha,beta)=(.6,.16) and −3.27759% for (−.6,−.16).
Independent fourth-order differences of the original rates and
entropy verify the analytic Jacobian on all nine states. Two
independent angle/energy orderings also agree.

The differential-reciprocity zero is distinct from the number and
timelike heat zeros at each of the tested finite displacements.
At beta=.12, the at-F=0, at-Q=0, and differential-reciprocity values of
alpha are approximately (−1.12234, −1.28205, −1.03794) at .30 GeV,
and (−1.93345, −1.94853, −2.05361) at .6882487043 GeV.
The slopes at common temperature independently reproduce the nonlinear
roots nearest equilibrium. The reciprocity slope is the negative
mean z of the positive blocking-difference weight at equilibrium;
the number and heat slopes use the original zeroth, first and second
equilibrium energy moments.

In the scanned .30−.68825 GeV interval the common-state reciprocity
slope crosses the number slope near .63860864 GeV and the heat slope
near .64414456 GeV. Both crossings are individually bracketed by the
same two nearby temperature cuts and independently resolved on the
two quadrature orders. The test finds these sign-changing slope
intersections in the sampled domain; it does not certify uniqueness
under arbitrary temperature or couplings.

At .64135 GeV and beta=±.002, two direct quadrature orders find the
reciprocity zero strictly between number and heat zeros. For positive
beta, (F,Q/T_SM) ≈ (−3.310e−22,+3.691e−21) in the existing rate units,
with positive entropy ≈1.870e−23. For negative beta both flux signs
reverse and entropy is again positive. Consequently, equality of the
two off-diagonal derivatives at an isolated nonequilibrium state
does not imply that either flux vanishes.

These are conditional local charged-lepton timelike benchmarks, not
elastic heat, a crossing temperature in an evolving universe, a
freeze-out abundance, or a proof of general linear stability.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_differential_reciprocity_reference.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_differential_reciprocity_reference.py --intermediate
```
