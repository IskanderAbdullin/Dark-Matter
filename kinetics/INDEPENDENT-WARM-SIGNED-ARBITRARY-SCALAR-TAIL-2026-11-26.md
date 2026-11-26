# A signed-pair depletion bound without a scalar occupation cap

The [previous analytic bound](INDEPENDENT-WARM-PAIR-LATE-SUPPORT-BOUND-2026-11-26.md)
put a uniform cap on the outgoing scalar occupation. The
[finite-support strengthening](INDEPENDENT-WARM-PAIR-SUPPORT-CLOSURE-2026-11-26.md)
additionally treated the exterior as unoccupied. Here, scalar
occupation may be arbitrarily large and have any momentum
support. It must carry only a finite part of the conserved
total dark number. The bound instead assumes a cap $F_\chi$
on each stable fermion occupation and a shrinking incident
stable support $p_{\rm lab}\le P/b$. As before, the only
channel allowed to change stable number is the tested signed
tree-level pair collision. No entropy injection or additional
number-changing channels are included.

For stable occupations $f_1,f_2$ and scalar occupations
$g_3,g_4$, the exact stable-loss factor is
$$
 f_1 f_2 (1+g_3)(1+g_4)-g_3g_4(1-f_1)(1-f_2)
 = f_1f_2(1+g_3+g_4)+g_3g_4(f_1+f_2-1).
$$
If $F_\chi\le 1/2$, the quadratic part cannot increase
stable loss, no matter how large the scalar occupation is.
This argument requires the **signed** collision kernel and
does not apply to a forward-only simulation with exact Bose
enhancement. We keep the linear terms and bound them using
scalar number, not occupancy.

For $p=P/b$ and $E_p=\sqrt{m_\chi^2+p^2}$, the outgoing
scalar energies supported by two incident stable particles
lie in $[E_-,E_+]=[E_p-p,E_p+p]$. Use the invariant measure
$d\Pi=d^3p/[(2\pi)^3 2E]$. The occupied scalar leg has
$\int d\Pi\,g\le n_\phi/(2E_-)$, and the other leg has
$\int_{E_-}^{E_+}d\Pi\le E_p\,p/(2\pi^2)$. Integrating the
two stable legs with fixed total four-momentum gives standard
two-body phase space $p_{\rm cm}/(4\pi\sqrt{s})\le p/(8\pi m_\chi)$.
These estimates need no scalar isotropy, smoothness or pointwise
upper bound. With the analytic tree envelope
$\overline{|\mathcal M|^2}\le C\,p^2$ from the previous bound,
and $f_1f_2\le F_\chi^2$, the physical loss density from
the linear terms is at most
$$
 L_{\rm lin}\le \frac{F_\chi^2 C\,p^4 E_p\,n_\phi}
                        {4\pi^3 m_\chi E_-}.
$$
The constant includes a factor of four for incoming spins,
one half for two identical scalars, two lost particles per
event, and two equal contributions from $g_3+g_4$.

Use $n_\phi\le N_{\rm dark,0}/b^3$ and the previous photon
floor $H\ge H_\gamma/b^2$. Then the integral of
$b^3L_{\rm lin}/H$ over $d\ln b$ is at most
$$
 D=\frac{F_\chi^2 C P^4 E_{P}N_{\rm dark,0}}
         {8\pi^3 m_\chi E_{-,0}H_\gamma}.
$$
Let $N_{\rm vac,lo}=N_{\chi,0}/(1+A)$ be the solution's
vacuum-only analytic lower bound, with $A$ the previous
quadratic-loss coefficient without Bose stimulation.
The positive deficit relative to that comparison solution
cannot grow faster than the additive linear term: a smaller
stable density reduces the quadratic loss. Hence
$N_\chi(\infty)\ge\max(0,N_{\rm vac,lo}-D)$.

| Anchor visible T | Exact maximum retained fermion f at anchor | Assumed future stable cap | Stable fixed-count Ωh² lower bound |
| --- | ---: | ---: | ---: |
| 1.206 MeV | 5.725e-6 | 1e-5 | 41.41 |
| 1.206 MeV | 5.725e-6 | 3e-5 | 13.53 |
| 0.414 MeV | 3.060e-6 | 3e-5 | 55.16 |
| 0.414 MeV | 3.060e-6 | 1e-4 | 41.09 |

A separate gross-rate integration tests two incident cutoffs and
two fermion occupancies. Scaling the trial scalar population by
0, 1 and 2 isolates the linear and quadratic forward terms.
The linear term is below 0.00111 of this analytic envelope in
those examples; the quadratic term is always smaller than its
matched reverse term. These numerical tests supplement the
inequalities and are not used to set the bound.

The 0.414 MeV bound stops excluding 0.120 once the hypothesized
cap exceeds about 1.91e-4. The cap and incident support must
hold **throughout** the continuation. The anchor measurements
do not prove either assumption under unmodeled kinetic reheating
or inverse production of stable particles at high momenta.
This is a restricted pair-channel result, not a model-independent
cosmological exclusion.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error signed_pair_arbitrary_scalar_tail_bound_reference.py --verify-collision \
  --anchor results/warm-coupled-pair-b10-28-shape-2026-11-26.json \
           results/warm-coupled-pair-b10-count-audit-2026-11-26.json chained_endpoint_cases \
  --anchor results/warm-coupled-pair-b30-32-shape-2026-11-26.json \
           results/warm-coupled-32-28-32-audit-2026-11-26.json stages
~~~
