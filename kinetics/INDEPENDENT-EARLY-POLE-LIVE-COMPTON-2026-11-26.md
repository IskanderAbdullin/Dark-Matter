# Live reciprocal hard Compton packets on the common energy budget

The [frozen hard Compton diagnostic](INDEPENDENT-EARLY-POLE-HARD-COMPTON-HEAT-2026-11-26.md)
is now placed on live photon and charged moments alongside
the cold scalar mode and optional photon pairs. The photon
window is [.5, 4], the fermion window [.45, 5]; both bounds are
imposed on both directions. The hard invariant/angle cuts and
vacuum amplitude are retained. The packet matrix elements
are massless and the exact charged mass enters the moment basis
and the conservative dressed energy remap.

One photon and one charged particle appear on either side of
each signed event. The final fermion radial coordinate is
determined from the variational energy and is accepted only
if an outgoing momentum triangle can close. Photon count and
both fermion counts are individually unchanged by Compton
collisions. A common-temperature fixture with independent
chemical potentials balances to numerical precision. The
full Bose and Fermi affinity makes the instantaneous entropy
production nonnegative. The finite same-state tests do not
supply a physical photon spectrum.

At finer packet resolution the initial heat into parent/daughter
is −9.275 / +7.607 (fine) and −9.492 / +7.496 (third).
The independent frozen integral with exactly these separate
photon and fermion bounds gives −9.343 / +7.660 and
−10.356 / +7.508, respectively. The hard reconstruction
difference is most visible for third parent, about 8%; the
live packet quadrature itself changes the four instantaneous
entries by less than 1.7% from the coarse to the refined rule.

On the fixed-scale exposure `.0003`, without the screened
unlike elastic channel, this hard Compton subset changes the
cold gain by:

| Archived state | Photon pairs off | Photon pairs on |
| --- | ---: | ---: |
| Fine | +5.80% | +5.48% |
| Third | +8.18% | +8.15% |

These are matched switches on the same archive and box, not
a calibrated physical rate. All sampled scalar increments
remain positive in these trials.

At the separate very short interval `2e-6`, adding Compton
changes cold gain by +.039% (fine) and +.050% (third) without
screened unlike elastic scattering. With the same screened
elastic calculation, it changes the gain by +.052% and +.039%,
respectively. Photon pairs are off in this short control.
This smaller numerical effect at the early interval does not
justify dropping the process everywhere, especially with an
evolving photon spectrum. The tested live channels share one
energy and entropy ledger; no final relic or late-universe
transfer function follows.

From `kinetics/` with one BLAS thread:

```sh
python early_pole_hard_compton_heat_reference.py --level 1 --bounds .45 5 --photon-bounds .5 4 --without-fixtures
python early_pole_archived_live_compton_reference.py fine --level 1
python early_pole_archived_live_compton_reference.py fine --evolve
python early_pole_archived_live_compton_reference.py fine --evolve --strength 0
python early_pole_cold_screened_split_reference.py third --duration 2e-6 --count 12 --compton-strength 1
```

Raw: `results/early-pole-live-compton-2026-11-26.json.gz`.
