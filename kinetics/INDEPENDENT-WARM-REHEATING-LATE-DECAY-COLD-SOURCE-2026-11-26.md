# How cold can an unscattered late two-body decay be?

For an unpolarized vacuum $H\to\chi\phi$ decay,
the scalar has rest-frame energy
$e_*=(M^2-m_\chi^2+m_\phi^2)/(2M)=48$ MeV.
For a parent Lorentz factor $\gamma$,
$E_\phi=\gamma e_*+p_*\sqrt{\gamma^2-1}\cos\theta$,
with uniform rest-frame $\cos\theta$ and
$p_*=\sqrt{e_*^2-m_\phi^2}$.
Given a lab energy cut $m_\phi<E_c<e_*$, all
parents with rapidity below
$\operatorname{acosh}(e_*/m_\phi)-\operatorname{acosh}(E_c/m_\phi)$
are excluded. Even allowing every remaining parent
an optimal boost, the cold fraction is bounded by
$$
 P_{\rm max}=\frac{E_c^2-m_\phi^2}
 {2p_*[p_*+\sqrt{e_*^2-E_c^2}]}.
$$
Indeed the unsaturated uniform angular fraction is
$1/2+(E_c-\gamma e_*)/(2p_*\sqrt{\gamma^2-1})$;
its derivative has the sign of $e_*-E_c\gamma$.
The maximum is at $\gamma=e_*/E_c$.
This applies to any parent energy spectrum with
the same isotropic vacuum rest-frame distribution,
not only to thermal parents.

At the visible 80 keV anchor, $q/T_{\gamma,0}=.1$
corresponds to lab momentum 8.24 keV. To put
a single daughter below that cut requires a parent
energy at least 291 GeV; even the best boost gives
a fraction below $7.38\times10^{-9}$.
The comoving cut translates to a smaller physical
cut for any later decay, so the anchor is the
most favorable value. The example threshold in
$q/T_{\gamma,0}$ is illustrative, not a data cut.

For a separate count illustration, grant
$\rho_{\rm anchor}/M$ heavy parents at the
anchor, irrespective of their actual composition;
add both the allowed charged-lepton energy grant
and the prior heavy-pair regeneration budget.
If every such heavy decays and there is **no repeated
inverse prompt formation**, the corresponding
scalar yield below the .1 cut is at most
$4.88\times10^{-14}$. For scale, this is about
$1.93\times10^{-11}$ of the full count required
for the scalar half-abundance target. Inverse
cycling can make the number of decay events exceed
this count; it is not included in the illustration.
The single-event vacuum angular result also does
not apply to decay angles reweighted by quantum
occupation, polarized parents, or collisional
cooling of daughters. It is not a bound on earlier
or coherent scalar production.

~~~sh
python -W error warm_prompt_reheating_late_decay_cold_source_reference.py \
 EOS.late.dat \
 results/warm-prompt-reheating-post-anchor-scalar-ceiling-2026-11-26.json \
 results/warm-prompt-reheating-cold-lepton-energy-allowance-2026-11-26.json
~~~

[Exact daughter kinematics and illustrative count ledger](results/warm-prompt-reheating-late-decay-cold-source-2026-11-26.json).
