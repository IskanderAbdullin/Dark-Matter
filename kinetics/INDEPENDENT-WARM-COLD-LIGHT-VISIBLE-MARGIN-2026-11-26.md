# Selected visible pair injection does not flip the frozen q32 cold-light sign

The [fixed-light recoil check](INDEPENDENT-WARM-NEAR-REST-LIGHT-INVERSE-DECAY-2026-11-26.md)
found a small negative decay/inverse difference at the sampled
points of the innermost q32 light shell. Add a separate upper on
specified **positive** vacuum visible pair injection, both directly
into light modes and indirectly through a stationary heavy parent.
Drop the visible reverse losses, which can only lower this comparison.
The sign result remains a sampled frozen check.

For a tagged dark energy (E), detailed balance bounds the
incoming thermal gain by (e^{-(E+E_2)/T}) times the vacuum
cross section for annihilation into the visible states, with
the specified pion statistical factor included in its cap.
Use one Maxwell partner charge (two spin states). Below
(sqrt s=1) GeV the dimensionless angle average is exact for
the bound proportional to (s=2d^2+2(EE_2-pp_2cos	heta));
above it keep half of the Boltzmann exponent to pay the invariant
threshold and use the global Breit–Wigner cap for the rest.
For charged pions the condition (sge4m_pi^2) bounds
the partner energy by ((4m_pi^2-2d^2)/(2(E+p))) from below;
use the earlier stated supported-channel cap for the cross section.
No new hadronic spectrum is assumed.

For the secondary path divide the heavy source bound by
((Gamma M/E_H)(1-max f_L)), a lower bound on the frozen
decay damping, and propagate that excess through the positive
light gain kernel. Use the maximum scalar occupation over
the entire represented entropy reconstruction as an upper,
not just a sampled maximum in the recoil window.

At the first saved point the direct selected visible upper
is (2.018	imes10^{-6}) per Hubble time and the secondary
positive upper is (2.12	imes10^{-13}). The largest
ratio of the combined imposed source to the smallest sampled
dark sink is 19.54% over all six frames. The augmented upper
remains negative at each tested point. A separate direct
rest-frame charged-lepton fold with Maxwell weights and
without final blocking gives (4.38	imes10^{-7}) per Hubble
time in the first frame, below the analytic upper.
The Maxwell partner moments independently agree with direct
radial quadrature to (8	imes10^{-16}).

This bound includes only the chosen thermal charged-lepton and
supported charged-pion vacuum channels. It does not cover
elastic gain, other dark pair production or any heavy lag,
initial heavy transients, unknown hadrons or exterior re-entry.
In particular the sign of the fully coupled evolution is not proved.

~~~sh
D=results
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_cold_light_visible_source_margin_reference.py \
 "$D"/warm-near-rest-light-inverse-decay-2026-11-26.json.gz \
 "$D"/warm-signed-fermion-with-charged-pions-2026-11-26.json \
 "$D"/warm-light-occupation-next-q32-tape-2026-11-26.json.gz \
 "$D"/warm-observed-suboperator-q32-tape-2026-11-26.json.gz
~~~

[Diagnostics](results/warm-cold-light-visible-source-margin-2026-11-26.json).
