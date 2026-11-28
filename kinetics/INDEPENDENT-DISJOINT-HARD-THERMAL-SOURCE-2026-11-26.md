# Disjoint hard thermal source sensitivity

The [free-quark Born comparison](INDEPENDENT-HARD-THERMAL-QUARK-BORN-2026-11-26.md)
overlaps the BESIII vacuum mass windows. Adding both in full counts
some invariant masses twice. This test avoids that overlap by using
the free thermal contribution **only outside** the measured windows
2.2324–3.08 and 3.40–3.671 GeV. The KLOE window ends below the
free-quark cutoff for all included temperatures. Keep the existing
FD lepton source, KLOE pion band, and separate BESIII vacuum intervals
with their original classical thermal weights. No J/ψ resonance
or near-pole CLEO proxy is included in these particular cases.

In the unmeasured mass regions, use the restricted thermal
Born shape (massless u/d/s, 4 ≤ invariant mass/T ≤ 12;
same boost support as the separate Born test). Compare two
choices for its lower visible-temperature cutoff.

| Positive prescription | T_R for 99% of the .120 optimistic scalar count |
| --- | ---: |
| Leptons and the specified measured sub-bands | 348.62 MeV |
| Also free-quark complement, acting from 300 MeV | 337.30 MeV |
| Also free-quark complement, acting from 250 MeV | 325.18 MeV |

At the measured-only root, 32.8% of the standalone hard Born
count from 250 MeV overlaps the two BESIII windows and is removed.
The corresponding hard count contributions are 0.0269144
before excluding this overlap and 0.0180849 afterward.
Changing quadrature order does not move the roots at the displayed
precision. Choosing the arbitrary lower-temperature cutoff
does, so the variation is a source-model sensitivity, not a confidence
interval or a measured reheating threshold.

This patchwork mixes a prescribed vacuum source with a finite-T
free-parton model on disjoint invariant-mass support. It does not
smoothly match QCD across the crossover, include a fully interacting
spectral function, or prove that all produced dark particles can
turn into scalars. As before, the optimistic count and heating
proxies are not a kinetic or cosmological solution.

Reproduce with the same EOS, KLOE data and BESIII PDF extraction as in
[the measured-source test](INDEPENDENT-PARTIAL-HADRON-REHEAT-SOURCE-2026-11-26.md):

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error disjoint_hard_thermal_measured_source_reference.py \
    EOS.dat KLOExSec.data BESIII-2112.11728.txt
~~~
