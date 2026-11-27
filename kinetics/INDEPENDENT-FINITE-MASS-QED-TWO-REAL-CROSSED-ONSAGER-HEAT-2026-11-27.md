# Two infrared-finite real channels give a positive number–heat response at equilibrium

Extend the [fixed-background photon heat comparison](INDEPENDENT-FINITE-MASS-QED-HARD-REAL-PHOTON-HEAT-2026-11-27.md)
with an independent linear-response control. Here all species start
at a common equilibrium temperature, not the nonthermal occupations
used in that comparison. Use the same masses and scalar momentum.

Perturb the scalar log inverse Bose odds by $X$ and the inverse
photon temperature by $Y$, keeping the charged distributions fixed.
For a scalar-producing event with a photon of energy $k$, use
$\sigma=+1$ for the charged pair channel (outgoing photon) and
$\sigma=-1$ for the parent-plus-scalar channel (incoming photon).
Its log forward/reverse odds change by $X+\sigma kY$.
Consequently, the linear model for scalar number and photon energy
has a matrix equal to the positive event equilibrium turnover
times the outer product of $(1,\sigma k)$, integrated over real
phase space. Each channel is symmetric and strictly positive
definite in these two forces; the spectrum contains more than one
photon energy. Their off-diagonal entries have opposite signs.

The independent centered perturbation of the exact FD/BE gain-loss
terms agrees with this matrix to 3.8e−8 relative. Angular endpoint
panels resolve the collinear features; the final two quadratures
differ by less than 7.9e−7 relative in the combined entries.
In the fine fixture at temperature 1, the pair susceptibility
$C_{NQ}=+11352.38$ is partly offset by the scalar absorption
$C_{NQ}=-1277.11$. The combined normalized number–heat correlation
is about .71 at temperature 1 and .59–.60 at .4 on the two mass
fixtures. Hence a single relaxation coefficient cannot replace
these two independent moments in this test.

These are two **unscaled one-charge real** channels at one fixed
scalar mode. No other collision channel, dressed quasiparticle
rate, entire two-plasma transport solution, or cosmological history
is inferred.

Replay from kinetics:
~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_real_crossed_onsager_heat_reference.py
~~~
[Fine receipt](results/finite-mass-qed-real-crossed-onsager-heat-fine-2026-11-27.json),
[third receipt](results/finite-mass-qed-real-crossed-onsager-heat-third-2026-11-27.json).
