# Separate gross visible stable-pair loss from net injection

Use the existing finite-width vector portal for direct χ₂χ̄₂ ↔
charged-lepton pairs, with the same 10 GeV mediator and visible
FD bath. The invariant and boost quadrature folds lepton blocking
into a positive angular weight while factoring out exp(−q/T).
For production, restore exp(−q/T) and the two dark Pauli terms.
For depletion, use the two **occupied** dark legs, without a
standalone exp(q/T). This obtains separate gross rates, not just
positive and negative parts of a net shell rate.

Split at both incoming dark shell boundaries and set occupation
outside the represented support to vacuum. The boost excess caps
are larger than necessary for all represented inverse pairs.
Double event rates for the combined stable fermion and antifermion
number. Results are normalized to its represented density and
the visible radiation Hubble rate.

| Archive | Gross direct visible production / (nχ H) | Gross direct depletion / (nχ H) |
| --- | ---: | ---: |
| 50 MeV reheat, 12 MeV quantum archive | 0.03797 | 3.331e-9 |
| 50 MeV reheat, 9.3 MeV quantum archive | 0.01470 | 2.164e-9 |
| 9.019 MeV reheat, paired quantum archive at 3.157 MeV | 0.007299 | 8.208e-13 |
| 50 MeV reheat, 12 MeV q32 populated control | 0.03801 | 3.328e-9 |
| Same q32 inner state, cropped at q16 | 0.03801 | 3.327e-9 |

The q32 archive can reach only 0.700 GeV in invariant mass,
far below this benchmark's 10 GeV mediator mass. Its added
occupied exterior changes the fixed-state direct visible inverse
number rate by about 0.02%, not by a resonant factor.
This does not ensure negligibility for a harder unknown tail
or a different mediator parameter. The results compare direct
e, μ and τ processes only, not the entire thermal portal or
a sustained background annihilation history.

```sh
OPENBLAS_NUM_THREADS=1 python -W error archived_visible_gross_depletion_reference.py EOS.late.dat results/{12-quantum-48-shape,early-9p3-quantum-shape,low-reheat-early-pair-quantum-shape,wide32-embedded-vacuum-shape}-2026-11-26.json --crop-last 16
```
