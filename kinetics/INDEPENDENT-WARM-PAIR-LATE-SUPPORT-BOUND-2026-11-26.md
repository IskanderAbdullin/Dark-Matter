# A late support-conditioned bound on stable pair depletion

Use the 1.206 MeV retained state and inclusive count from the
[32+28 time-step warm history](INDEPENDENT-WARM-COUPLED-B10-2026-11-26.md).
Let the stated tree-level signed $\chi\bar\chi\leftrightarrow\phi\phi$
reaction be the only process changing the inclusive stable count.
Assume incident stable particles satisfy $p_{\rm lab}\le P/b$,
with $P=17.466$ MeV and $b=a/a_{\rm anchor}$; the photon clock
satisfies $bT_\gamma\ge T_{\rm anchor}$; and previously excluded
particles do not re-enter. Let $F$ bound the scalar occupation
at forward final momenta for the entire continuation, not just at
the anchor. The fixed-count proxy also assumes unchanged visible
comoving entropy.

For equal incoming masses and equal outgoing masses, $t,u\le0$,
so the product of the heavy-exchange denominator magnitudes is at
least $M^4$. In the [tree pair amplitude](INDEPENDENT-SCALAR-PAIR-ARCHIVED-RATES-2026-11-26.md),
apply $|x\sqrt{1-x^2}|\le 1/2$, $k^2\le m^2+P^2$ and
$p_{\rm cm}\le P/b$. This gives
$\overline{|\mathcal M|^2}\le C p_{\rm cm}^2$,
$C=1.4570\times10^{-9}\ {\rm GeV}^{-2}$. The separate numerical
angle-grid check uses at most 0.732 of this analytic bound.

The identical final scalars give
$$
\sigma v_{\rm M}
=\frac{k}{64\pi\sqrt{s}E_1E_2}\int_{-1}^1dx\,\overline{|\mathcal M|^2}
\le\frac{CP^2}{64\pi m^2b^2}.
$$
Their Bose factors give at most an additional $(1+F)^2$.
For inclusive comoving stable number $N$, the incident charge
densities satisfy $n_\chi n_{\bar\chi}\le (N/(2b^3))^2$.
Each forward reaction removes two fermions; dropping inverse
production can only weaken the lower bound. Photons alone
give $H\ge H_\gamma/b^2$. Integrating the resulting inequality
for $1/N$ gives
$$
\frac{N(\infty)}{N(1)}
\ge\left[1+\frac{N(1)CP^2}{384\pi m^2H_\gamma}(1+F)^2\right]^{-1}.
$$

The coefficient before $(1+F)^2$ is 0.279727. At $F=1$, the
conditional stable proxy is at least 27.11, versus 0.120;
at $F=0.001$ the floor is 44.87. Reaching .120 requires that
the permitted occupation supremum reach roughly 40.3 (or a
hypothesis fail). This is necessary, not sufficient. The
resolved peak at the anchor is 7.06e-4; this does not establish
an occupation bound for the entire future.

Additional number-loss channels, changing momentum support,
stronger Bose stimulation, exterior re-entry or nonstandard
photon and entropy histories can evade this restricted result.
It does not exclude other reheating histories or determine an
exact asymptotic yield.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error warm_scalar_pair_late_support_bound_reference.py \
    results/warm-coupled-pair-b10-28-shape-2026-11-26.json \
    results/warm-coupled-pair-b10-count-audit-2026-11-26.json
~~~
