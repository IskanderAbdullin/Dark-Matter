# Expanding resolved photon populations retain the conditional gate

Evolve hard photon occupations alongside the early scalar
spectrum and joint fermion count/energy moments. Use the same
pair matrix integrals as the [nonexpanding photon box](INDEPENDENT-EARLY-POLE-SPECTRAL-PHOTON-BOX-2026-11-26.md).
The photon modes outside the modeled hard set remain spectators
in comoving variables. Each pair event conserves number plus
energy between photon and fermion sectors; the common FD closure
for the fermions is still imposed. The expanding continuity
ledger includes the scalar and both fermion traces.

At g=.005, the reversible resolved-photon model crosses the
soft-mode stimulated gate near x=.050928. With 20 instead of 12
photon nodes it crosses near .050932, and with the other saved
scalar/fermion grid near .050922. The matched number-budget BE
photon closure crosses around .050885. At g=.004 and .003 it
does not cross by .052 in these same conditional hard-cut
runs. The final hard photon population differs from the full
moment-matched BE target by about 6.5–6.9% of total photon
count in the sampled active-bin absolute measure.
Imposing extra count/energy-conserving hard-bin redistribution
at rate 100 per reference expansion time barely changes the
.005 tested gate in this particular time interval.

The sign change survives removing whole-spectrum fast thermal
photon redistribution for this one model family. It is not a
calibrated threshold in the real gauge theory. A single-line
damping scale at g near .004–.005 can be comparable to the
kinematic angular resolution of the narrow pole. Establishing
the actual two-line spectral/vertex kernel, and including
soft/collinear processes consistently, remain separate work.
No final abundance or transfer spectrum is inferred.

~~~sh
python -W error early_pole_spectral_photon_decay_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --order 20
~~~

[Raw coupled shape/energy/count checks](results/early-pole-spectral-photon-decay-2026-11-26.json.gz).
