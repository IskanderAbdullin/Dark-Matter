# Evolve the resolved scalar-pair channel on expanding archived states

Start from the represented φ and stable fermion momentum cells. Evolve
only the signed χ₂χ̄₂ ↔ φφ channel, physical momentum redshift, and a
prescribed visible-radiation expansion. Use an explicit midpoint
collision update between two conservative redshift steps. Each
collision changes species number but conserves their sum and energy.
Count daughters outside the represented shells in a separate
collisionless energy ledger; these cohorts do not feed back into
the rates. No source from the visible bath, additional dark
channels, or pre-checkpoint evolution is added here.

The positive and negative rates are computed separately from the
actual within-cell occupations. The original scalar first-bin
momentum cut is not inside the low-reheat input support, so no
cosmological cold fraction can be inferred from the low-reheat runs.

## Low-reheat archive after 3.157 MeV

- At fixed background for 0.05 initial Hubble time, the six-step test
  raises the retained scalar count 3.885-fold and the inclusive count
  3.891-fold. The retained mean scalar kinetic energy reaches
  23.97 MeV and the stable count decreases only 1.30e-5.
- With expansion to a/a_archive=10, the 12/24-step runs give inclusive
  scalar count increases of 24.60/24.74 times the initial population,
  and retained-only increases of 21.86/22.00. The finer inclusive
  scalar rest-mass proxy is 5.638e-11. Stable loss is 1.064e-4 of the
  archived count, leaving stable Ωh²≈0.1182255. About 11% of the final
  scalar count was produced outside this finite support. An expanded
  grid to q/T_R=128 starts with the identical inner state and vacant
  outer cells. It retains essentially all of these particles and
  reproduces the inclusive scalar count to 1.4e-9 relative at
  matching temporal resolution. The maximum inverse-to-forward
  ratio in the wide low run was below 5e-9.
- A matched radiation-only redshift control conserves both species'
  counts. The independent vacuum-only fold on the collisionlessly
  redshifted incoming fermion shape gives 5.649e-11 at scale factor
  10 and 5.657e-11 at scale factor 60. The small low-reheat scalar
  yield is thus sensitive to this additional source, but its
  enhancement is not evidence for a dominant scalar abundance.
  The fold discards inverse, stimulation, and all other channels,
  and uniformly scales the fermion count with the one-rate ODE
  solution; it is not a second time-evolved shell solution.
- A separate deliberately optimistic bound places the entire
  low-reheat fixed stable source budget Ωh²=0.120 at reheating,
  admits maximal arbitrary-momentum vacuum annihilation, and
  uses only photon energy density as the radiation floor to the
  Hubble rate. It permits a scalar rest-mass increment at most
  4.88e-9, for **vacuum annihilation only**. This upper bound is
  not a time-extended bound on stimulated or different channels.

## High-reheat archive after 12 MeV

At a/a_archive=1.2 with inverse reactions enabled, the 5/10-step
runs lose 4.65%/4.70% of the archived stable count and gain roughly
12.2% scalar count. A control artificially disabling inverse
reactions loses 18.76%/18.78% of stable count instead. The comparison
shows that its hard initial scalar population affects the evolution,
not just the instantaneous slope. The absolute comoving scalar
count in the original q/T_R < 0.005 band is unchanged at roundoff in
both controls, even though the inclusive scalar count grows.
The retained exterior is not
allowed to re-enter, and omitted prompt, elastic, and visible
operators can change both histories. It does not close the
high-reheat abundance or provide a full dark-sector spectrum.

## Reproduce

Use the [low-temperature radiation table](INDEPENDENT-BENCHMARK-SETUP.md)
from `kinetics`:

```sh
export OPENBLAS_NUM_THREADS=1
low=(python -W error scalar_pair_expanding_continuation_reference.py EOS.late.dat results/low-reheat-matched-vacuum-shape-2026-11-26.json)
"${low[@]}" --stop .05 --fixed --steps 3 6
"${low[@]}" --stop 10 --steps 12 24
"${low[@]}" --stop 10 --pad 128 --steps 12 24
"${low[@]}" --stop 10 --pad 128 --steps 12 24 --no-pair
high=(python -W error scalar_pair_expanding_continuation_reference.py EOS.late.dat results/12-quantum-48-shape-2026-11-26.json --stop 1.2 --steps 5 10)
"${high[@]}"
"${high[@]}" --no-inverse
python -W error scalar_pair_vacuum_redshift_fold_reference.py EOS.late.dat results/low-reheat-matched-vacuum-shape-2026-11-26.json --optimistic-total-stable .12
```
