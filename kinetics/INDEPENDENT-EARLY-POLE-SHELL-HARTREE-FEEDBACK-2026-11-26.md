# A newly populated cold shell can move its own pole

The [fixed mass clock](INDEPENDENT-EARLY-POLE-DECAY-MODE-CLOCK-2026-11-26.md)
does not justify assigning the equilibrium thermal
mass to an initially empty scalar sector. In a
restricted on-shell closure add only a positive
vacuum-subtracted real-scalar quartic tadpole,
using the standard $\lambda\phi^4/4!$ convention
and its [thermal normalization](https://academic.oup.com/ptep/article/2021/2/023B09/6020275):
$$
 m_{\rm pole}^2=m_{\rm vac}^2+
 \frac{\lambda}{2}\int \frac{d^3p}{(2\pi)^3}
                         \frac{f(p)}{\sqrt{p^2+m_{\rm pole}^2}} .
$$
This leading Hartree closure is evaluated
self-consistently as the shell fills. It uses
neither a vacuum-to-equilibrium mass jump nor the
fermion static curvature as a pole. An equilibrium
massless BE bath would give $\lambda T^2/24$;
an initially empty shell gives no occupation
tadpole. Damping and nonlocal terms are not included.

The shell shares a single occupation over momenta
from .5 to 1.5 times $p_*=4.1\times10^{-6}T$.
Compute the number-weighted clock from the finite
pole interval, keep the parent FD reservoir fixed
and impose a radiation-era $T\propto a^{-1}$.
The rest width is scaled with the illustrative
scalar Yukawa Born amplitude relative to the
massless-daughter value; no physical bath coupling
is inferred. For a comparison target only, put 92%
of a chosen final full-abundance yield in this
very narrow band. This is **not** the upstream
fraction measured in that particular band; the
upstream public hot cut used $p/T=.1$.

| Held scenario | Fraction of illustrative narrow shell target at $M/T=.05$, or earlier reach |
|---|---:|
| No self tadpole, empty fermion daughter | reaches target at $M/T=.00720$ |
| $\lambda=7.65\times10^{-8}$, empty massless fermion daughter | .1724 |
| $\lambda=1.085\times10^{-6}$, empty massless fermion daughter | .01216 |
| $\lambda=7.65\times10^{-8}$, empty fermion pole $b/T=.003$ | .2105 |
| $\lambda=7.65\times10^{-8}$, held massless FD daughter with $\mu/T=-3\times10^{-5}$ | $9.51\times10^{-16}$ |

The last row is a nearly replenished daughter-bath
extreme, not a trajectory of the asymmetric portal.
In the second row the final scalar pole is
$m/T=.001128$, much larger than the shell momentum,
while the on-shell decay is still open. The
negative feedback acts against the explosive
early solution even without closing phase space.
The calculation is not a closed number-and-energy
system, does not evolve other momenta or the
daughter distribution, and cannot establish
an actual in-medium mass, spectrum or relic.
A shape-independent check within the *same* positive
tadpole assumptions is separate.

~~~sh
python -W error early_pole_shell_hartree_feedback_reference.py
~~~

[Shell histories and gap equation controls](results/early-pole-shell-hartree-feedback-2026-11-26.json).
