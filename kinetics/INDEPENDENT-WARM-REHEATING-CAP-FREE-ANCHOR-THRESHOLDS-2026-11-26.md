# Locate the cap-free obstruction on the sampled late history

The [first-exit survival bound](INDEPENDENT-WARM-REHEATING-CAP-FREE-LATE-BARRIER-2026-11-26.md)
retains at least .826 of the stable fermions after the
modeled visible 80 keV anchor without a pointwise scalar
occupation cap. Here we reuse the same rate inequalities
but compute a better passage-time lower bound, and check
the warmer anchors of the existing thermal diagnostic.
This stage does not evolve a new spectrum.

For $X=b^3n/n_0$ and $R=\rho_{\rm cap}/(m_\chi n_0)-X$,
the previous dark energy budget and photon expansion floor
give $dX/du\ge -G(X,u)$, where $u=1-1/b$ and
$$
 G(X,u)=A X^2+\min(BX,CX^2)+D R^2
         +S R X^2(1-u)^3+J R X^{4/3}(1-u).
$$
Here $A$ is the stable vacuum loss, $B,C$ bound the
same three-lepton process two ways, $D$ is the scalar
quadratic term by reversal, and $S,J$ are the soft and
hard linear terms. In the softer term the scalar
density brings one more $b^{-3}$, while in the harder
term the two Pauli integrals have combined $b^{-4}$
instead of $b^{-3}$. Both additional dilution factors
can be dropped to form an autonomous upper loss rate.

To prove a bound at a test value $x$, split $[x,1]$
into slabs; on each use the highest $X$ for all
positive powers of $X$ and the lowest $X$ for the
largest $R$. Sum width divided by this upper rate.
If the sum reaches 1, the entire photon-normalized
future is too short to pass below $x$, even if the
real trajectory is nonmonotone. This certifies a
conditional floor using positive rectangles rather
than interpolating the thermal history or trusting
a fast loss integration. The stored 128/2048 slab
arms, direct quadrature and a continuous comparator
with additional dilution give independent numerical
controls. A tiny outward numerical pad is applied
to the rectangle result.

| Sampled visible temperature | Stable retained fraction lower | Stable \(\Omega h^2\) lower |
|---:|---:|---:|
| 187.49 keV | No positive bound | No positive bound |
| 166.00 keV | .05188 | 8.69 |
| 147.23 keV | .45719 | 76.62 |
| 130.71 keV | .63332 | 106.13 |
| 102.84 keV | .79584 | 133.36 |
| 80.00 keV | .87049 | 145.86 |

166 keV is the warmest **sampled** anchor at which
this test certifies more than .12. There is no
claim that the physical branch has a transition
at that temperature or that a warmer inconclusive
comparison permits efficient depletion. At the
166 keV sampled anchor the inequality still clears
.12 if the energy allowance is at most about
1.0068 of nominal, keeping the initial stable
count fixed. At 80 keV this multiplier is 1.6764;
the stronger claim that at least 80% of the count
survives holds up to 1.1677 of the nominal budget.
Increasing a comparison budget does not solve the
production history or alter the copied input count.

The same qualifications as the original bound
apply: ordinary expansion/entropy, the selected
vacuum reactions, homogeneous isotropic and
factorized resolved distributions, and no coherent
condensate or added energy source outside the budget.
The numbers are conditions on modeled moments, not
evidence that those moments are reached.

~~~sh
python -W error warm_prompt_reheating_cap_free_thresholds_reference.py \
 EOS.late.dat \
 results/warm-prompt-reheating-post-2mev-elastic-2026-11-26.json.gz \
 results/warm-prompt-reheating-cap-free-late-barrier-2026-11-26.json
~~~

[Sampled anchor and energy robustness ledger](results/warm-prompt-reheating-cap-free-anchor-thresholds-2026-11-26.json.gz).
