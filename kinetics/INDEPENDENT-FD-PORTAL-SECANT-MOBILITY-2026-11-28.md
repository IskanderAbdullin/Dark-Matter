# Positive nonlinear secant response does not imply reciprocal differential response

Use the same charged-lepton FD model, fixed visible temperature, force
coordinates u=(alpha,beta), and flows (F,Q/T_SM). The affinity at one
event is d=alpha+beta*z, where z=E_lab/T_SM. Divide the signed local
rate by d, with the continuous positive limit at d=0. This weight
is the logarithmic mean of the two event directions. Then

    M(u) = integral [ J(d)/d ] * [1,z] * [1,z]^T
    flows(u) = M(u)*u
    entropy(u) = u^T*M(u)*u.

The positive matrix is symmetric at every state, because it is a
weighted Gram matrix. It is a secant from equilibrium, not the
differential Jacobian. Differentiating M(u)*u adds state dependence
from the energy and blocking weights. This allows a nonzero force
curl, even when entropy stays nonnegative. At the common detailed
balance state the extra derivative term is zero and the tangent
and secant matrices coincide.

The independently paired numerical controls reconstruct the
actual rates and entropy to less than 8e−17 in the stated scaled
norm, with positive eigenvalues. At T_SM=.05 GeV and forces (.6,.16),
the largest tangent-minus-secant component is about 1.67 times
the largest secant component. Near the selected local reciprocal
but nonequilibrium point at .64135 GeV the two matrices are close
but not equal. This finite model is not a trajectory or expansion
stability test.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_secant_mobility_reference.py
```
