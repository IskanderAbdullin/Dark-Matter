# How high must reheating be even to supply enough scalar counts?

Extend the prescribed-radiation **vacuum-final lepton source** through
the QCD-temperature region. The visible central radiation table is used
for its entropy clock. The collision source itself still has only
electrons, muons and taus; it does not include hadron production.
Count both original stable daughters and the one stable plus one scalar
from each promptly decayed heavy daughter. To test whether the count
is enough, assign **every** dark particle from this source to an 86 eV
scalar, with no number multiplication or change in visible entropy.

| Target fraction of observed rest-mass Ωh²≈.120 | Prescribed T_R | Stable no-loss Ωh² proxy at that T_R |
| --- | ---: | ---: |
| 50% | 310.85 MeV | 9,319 |
| 99% | 421.43 MeV | 18,436 |
| 100% | 422.67 MeV | 18,622 |

The two momentum and time quadratures give consistent thresholds.
These are **necessary source-count thresholds only within this lepton
history**, not viable reheating temperatures or a full source bound.
Changing the source, subsequent dark number growth, visible entropy,
or mass changes the threshold. At the 99% threshold even the
optimistic all-conversion allocation says nothing about whether
the scalar particles are cold.

A separate positive energy fold redshifts every deposited energy
element as radiation until photon temperature 1 MeV. If the
prescribed source energy stays in the dark sector and real dark
particles have nonnegative mass, their subsequent trace contribution
can only increase this radiation-redshifted value. At the 99% count
threshold the resulting conditional energy floor is 1.34% of visible
radiation. If stable counts instead remain fixed, their rest energy
alone is about 8.91% of visible radiation at 1 MeV. These are
different conditional histories, not a bound on scenarios where
dark energy returns to the visible sector or affects the source
or expansion. A matter or nucleosynthesis likelihood is not computed.

```sh
OPENBLAS_NUM_THREADS=1 python -W error scalar_count_reheat_target_reference.py EOS.late.dat
```
