# Ordinary scalar chemical equilibrium can supply count without a cold spectrum

[Relativistic real-scalar chemical kinetics](https://arxiv.org/abs/1906.07659)
distinguishes chemical equilibrium from a thermal shape
with a nonzero effective chemical potential. This note
uses the narrower **normal-phase ideal Bose chemical
equilibrium** limit, $\mu_\phi=0$, of a gas whose
number-changing reactions remain active before a
hypothetical sudden decoupling. It does not evaluate
whether the project's quartic or other amplitudes can
reach that state.

## Equal-count ordering at arbitrary mass over temperature

Write $\theta=a_f T_f$, $r=m_\phi/T_f$, and
$f_r(q)=[\exp\sqrt{q^2/\theta^2+r^2}-1]^{-1}$.
For fixed comoving number $n_c$, solve
$$
 \theta=\left[\frac{2\pi^2n_c}{I_2(r)}\right]^{1/3},
 \qquad I_j(r)=\int_0^\infty du\,
                 \frac{u^j}{\exp\sqrt{u^2+r^2}-1}.
$$
Since $I_2$ decreases with $r$, $\theta$ increases.
For any $r_2>r_1$, the difference between the
**squares** of the two exponential arguments is
$(r_2^2-r_1^2)-q^2(\theta_1^{-2}-\theta_2^{-2})$.
It starts positive and crosses zero once.
Thus $f_{r_2}$ is smaller below the crossing and
larger above it, with the same total radial number.
The fraction below **every** comoving momentum
cut cannot increase with $r$. Equivalently,
the relativistic equilibrium shape minimizes every
increasing momentum cost in this fixed-count family.
This is a statement about the frozen comoving PSD,
not about the amount of time a later-decoupled
species spends free streaming. It does not compare
arbitrary chemical potentials, entropy histories,
ongoing collisional perturbations or coherent modes.

At $\Omega_\phi h^2=.060$, $m_\phi=86$ eV and the
adopted present entropy, even the limiting relativistic
equilibrium spectrum has $\langle q/T_{\gamma,0}\rangle=.890$
and only 1.73% of its count below $q/T_{\gamma,0}=.1$.
For $r=1$ the cut fraction is .218%, and for $r=3$
it is .0202%. A large occupation right at rest is not
a substitute for a substantial cold count.

## Energy comparison at the modeled 80 keV anchor

A hypothetical equilibrium gas with the scalar
half-abundance count needs dark temperature 27.16 keV.
Its energy is 28.3% of the recorded aggregate dark
energy, or 3.93 times the recorded scalar energy.
If all recorded dark energy were converted to an
ideal $\mu_\phi=0$ gas without expansion or losses
during that conversion, it could supply a hypothetical
frozen scalar $\Omega_\phi h^2=.1545$. The number
deficit of the older restricted reaction network
is therefore not a shortage of total dark energy:
it is a restriction on available routes and net
scalar counts. This instantaneous scenario assumes
unspecified efficient conversion and relinquishes
the original stable-fraction and scalar count
ledger; it is not a solution.

At the required half-abundance thermal shape,
forcing collisionlessness immediately from the
80 keV anchor gives a 5.54 Mpc mean ballistic
path in the fixed reference background. This
checks the warmth of a particular hypothesized
shape; a path is not a matter power likelihood.

~~~sh
python -W error scalar_chemical_equilibrium_cold_order_reference.py \
 EOS.late.dat \
 results/warm-prompt-reheating-post-2mev-elastic-2026-11-26.json.gz
~~~

[Number, energy and momentum checks](results/scalar-chemical-equilibrium-cold-order-2026-11-26.json).
