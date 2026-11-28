# A shell-shape-independent upper on the added retained dark pair gains

The [tagged dark-pair fold](INDEPENDENT-WARM-COLD-LIGHT-PARENT-DARK-PAIR-2026-11-26.md)
is sensitive to the choice of reconstructed shape. Obtain a
uniform bound on those two **added positive** heavy sources
on the whole parent recoil interval of the innermost light
shell, without fixing either incident shape inside its shells.
Retain the archived support and shell numbers.

Fix heavy energy (E_H) and momentum (p_H). One incident
particle must carry at least (E_*=(E_H+M)/2), from energy
conservation with the other heavy output. Its (N_*) is at
most the sum of counts of shells touching that energy range.
For either incoming mass, integrate the remaining heavy
momentum using the three-momentum delta, then integrate the
orientation of the incident total (mathbf q) using the
energy delta. The latter yields (1/(q p_H)) after cancelling
the untagged heavy energy in the invariant measure. Change
incident opening cosine to (q), with
(dcos	heta=q,dq/(p_1 p_2)). The allowed (q) interval
has width no more than (2min(p_1,p_2)) by the triangle
inequality. Thus the two remaining radial integrals can
be bounded by (N N_*/p_*), where (N=int p^2f,dp)
and (p_*=sqrt{E_*^2-d^2}); count both choices of the
energetic incident leg. The squared matrix element divided by
(E_1E_2) is bounded cell by cell on ((s,cos^2	heta_{m CM}))
for the scalar channel and by the negative (t) endpoint
for the light channel. Four times the spin-averaged heavy
matrix element bounds the spin sum even without using the
identical-initial scalar factor.

To turn the rate into a bound on an additional stationary
heavy occupation, use the positive decay response. Pauli
(f_Lle1) and the represented radial number alone give
(int p f_L,dple(3N_L)^{2/3}/2), from filling the lowest
momenta first. Converting that integral to the light-energy
average bounds the heavy decay blocking without imposing
an additional pointwise half-occupation threshold. Scalar
occupation only increases the decay damping in this linear
response.

On the first q32 endpoint, the added stationary upper is
(2.47	imes10^{-22}) for the two retained pair gains combined.
Adding it to the **frozen reference** with the
[selected visible corrections](INDEPENDENT-WARM-COLD-LIGHT-VISIBLE-MARGIN-2026-11-26.md)
uses less than 19.6% of the remaining sampled negative
margin across the six q32 endpoints. The upper is much
looser than the explicitly reconstructed rates, but it
does not depend on their intra-shell profiles. The q16
sources remain absent by energy support.

Only this extra positive pair source bound is shape independent.
The baseline decay and inverse terms are still frozen to the
archive. This does not prove the full rate for re-shaped
profiles or establish the effects of unrepresented tails,
elastic scattering or an evolving heavy distribution.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_cold_parent_retained_pair_envelope_reference.py
~~~

[Diagnostics](results/warm-cold-parent-retained-pair-envelope-2026-11-26.json).
