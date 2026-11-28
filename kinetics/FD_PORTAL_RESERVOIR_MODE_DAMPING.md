# Physical equilibrium modes of the local FD collision steps

At fixed visible temperature, write the local force as `u=(alpha,beta)`. At detailed balance the linear decay generator in these coordinates is `B=A^{-1} diag(2,T_SM) J`, where `A` is the Jacobian from forces to dark number/energy moments, and `J` is the full equilibrium FD flux Jacobian. The code checks positive rates and independently cross-checks them on two phase-space meshes.

For a positive mode with rate lambda and `z=dt*lambda`, the symmetric moment-chord update has linear damping factor `(1-z/2)/(1+z/2)`. Backward Euler has factor `1/(1+z)`; the exact smooth linear flow has `exp(-z)`. Thus the symmetric method changes the mode's sign for `z>2`, and its magnitude approaches 1 as `z` increases. A discrete entropy identity alone does not imply strong stiff-mode damping. The executable uses signed nonlinear perturbations about the same FD equilibrium, not a separately manufactured quadratic model. Central perturbations at two amplitudes isolate the linear response and check cross-mode leakage.

In these three samples, the fast-to-slow rate ratio ranged from 4.44 to 258.94. At fast-mode `z=80` the numerical response was -0.951220 for the symmetric step and +0.0123457 for backward Euler; exact linear damping is about 1.80e-35. Independent-mesh decay rates agreed within 2.1e-9 in relative units.

Run from `kinetics`:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_reservoir_mode_damping_reference.py
```

This is a fixed-reservoir local stability diagnostic. It does not prescribe global timesteps or imply a cosmological result.

The distinction between dissipation and stiff damping is known for discrete-gradient schemes; see [Dissipative numerical schemes on Riemannian manifolds with applications to gradient flows](https://arxiv.org/abs/1804.08104). The rates and nonlinear mode tests here are independently computed from the FD portal and are not taken from that paper.
