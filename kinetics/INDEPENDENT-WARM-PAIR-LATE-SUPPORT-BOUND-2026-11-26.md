# A late support-conditioned bound on stable pair depletion

Use the 1.206 MeV retained state and inclusive count from the
[32+28 time-step warm history](INDEPENDENT-WARM-COUPLED-B10-2026-11-26.md).
After this anchor assume: only the stated tree-level signed
(chiarchileftrightarrowphiphi) interaction can change the
inclusive stable count; incident stable particles remain within
(p_mathrm{lab}le P/b), (P=17.466) MeV and (b=a/a_mathrm{anchor});
and the photon temperature satisfies (bT_gammage T_mathrm{anchor}).
Previously excluded particles do not re-enter the reaction. Bound
the scalar occupation at forward final momenta by a time-uniform
(F). The snapshot's peak alone does not establish this last assumption.
The fixed-count (Omega h^2) normalization also assumes no change in
comoving visible entropy.

For equal incoming masses and equal outgoing masses, (t) and (u)
are nonpositive, so the product of the two heavy-exchange denominator
magnitudes is at least (M^4). In the spin-averaged amplitude used in
[the pair benchmark](INDEPENDENT-SCALAR-PAIR-ARCHIVED-RATES-2026-11-26.md),
bound (|xsqrt{1-x^2}|le 1/2), (k^2le m^2+P^2), and
(p_mathrm{cm}le P/b). This gives
(overline{|mathcal M|^2}le C p_mathrm{cm}^2),
(C=1.4570	imes10^{-9} {m GeV}^{-2}). This is an analytic
envelope; the separate angle-grid check uses at most 0.732 of it.

The identical final scalars give
[
sigma v_mathrm M
=rac{k}{64pisqrt{s}E_1 E_2}
    int_{-1}^{1} dx,overline{|mathcal M|^2}
le rac{C P^2}{64pi m^2 b^2}.
]
Both Bose factors are at most ((1+F)^2) together. If the total
inclusive comoving stable count is (N) measured per physical volume at
the anchor, then the incident charge densities satisfy
(n_chi n_{archi}le (N/(2b^3))^2). Each forward collision removes
two fermions; omitting inverse production can only weaken a lower bound
on (N). With (Hge H_gamma/b^2) from photons alone, these relations
give the integrated inequality
[
 rac{N(infty)}{N(1)}ge
 left[1+rac{N(1) C P^2}{384pi m^2 H_gamma}(1+F)^2ight]^{-1}.
]
The coefficient before ((1+F)^2) is 0.279727. Even taking (F=1)
everywhere yields conditional (Omega_chi h^2ge 27.11), versus
0.120. For (F=0.001) the lower bound is 44.87. Reaching .120
requires that the allowed occupation supremum exceed roughly 40.3
(or a hypothesis of the bound fail); this condition isn't sufficient
for such depletion. For comparison, the anchor peak in the
resolved cells is (7.06	imes10^{-4}). No estimate of a thermal
or actual late high-occupation tail has been inferred from that peak.

A growing comoving incident support, a new number-loss channel,
strong Bose growth outside the adopted cap, or nonstandard visible
expansion or entropy can invalidate this restricted result. It
does not exclude other reheating histories or establish the
exact asymptotic yield of this one.

```sh
OPENBLAS_NUM_THREADS=1 python -W error warm_scalar_pair_late_support_bound_reference.py \
    results/warm-coupled-pair-b10-28-shape-2026-11-26.json \
    results/warm-coupled-pair-b10-count-audit-2026-11-26.json
```
