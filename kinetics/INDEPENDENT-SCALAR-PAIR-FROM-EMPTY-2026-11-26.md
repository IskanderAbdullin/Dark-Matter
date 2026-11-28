# Resolve the low-reheat pair source from the initially empty stage

At T_R=9.0190699 MeV start with the same empty physical pulse, finite
support q/T_R in [0.8,16], and energy-coupled visible radiation as
in the [matched source benchmark](INDEPENDENT-LOW-REHEAT-FRW-2026-11-26.md).
Add a signed χ₂χ̄₂ ↔ φφ midpoint on the fitted distributions at each
source split. Opposite species number changes are tracked separately;
the sum reconciles with exterior daughters. The initial tiny source
pulse does not apply annihilation to a nonexistent population.
Exterior cohorts contribute to the expansion and ledger but do not
return to the collision integrals. This is not a chemical network
including reactions that multiply total particle number.

| At 3.157 MeV | Steps | Inclusive scalar fixed-count Ωh² | Inclusive stable fixed-count Ωh² |
| --- | ---: | ---: | ---: |
| Vacuum lepton source, no dark collision | 24 | 2.2801e-12 | 0.11977062 |
| Same lepton source plus signed scalar-pair conversion | 12 | 2.4102e-10 | 0.11915628 |
| Same lepton source plus signed scalar-pair conversion | 24 | 2.4507e-10 | 0.11971416 |
| Same lepton source plus signed scalar-pair conversion | 48 | 2.4609e-10 | 0.11985385 |
| Signed visible portal, prompt return and elasticity, plus signed scalar-pair conversion | 24 | 2.4506e-10 | 0.11971416 |

Compare at the **same** mesh because source time quadrature is larger
than several physical corrections. At 24 steps the scalar proxy is
about 107 times the matched source-only case, while the stable change
is less than 0.05% relative. The richer 24-step control gives almost
the same integrated count for this low reheating; this does not
generally justify neglecting the moment exchange or reverse portal.
The dual visible-dark-exterior energy ledger and total number
reconciliation are independently checked.

Restart the 48-step earlier vacuum-visible-plus-pair endpoint without
any thermal replacement. Through a/a_archive=10 retain only signed
pair collisions and prescribed radiation, and carry the earlier
exterior counts collisionlessly. The chained inclusive proxy is
Ωh² ≈ 3.0031e-10 for scalars and 0.1198412 for stable fermions.
Compared with the already pair-populated checkpoint this late
contribution is about 22%, not a factor of 25. The previously
reported factor of ~25 used an initial archive with the earlier
pair source switched off. The new numbers are rest-mass count
proxies, not cold fractions or a complete relic prediction.

A separate early segment to a/a_archive=4 gives a padded q/T_R=64
cross-check. Neither grid includes the original very soft q/T_R<.005
band. No low-reheat condensate claim follows.

## Reproduce

From `kinetics`, prepare the [EOS table](INDEPENDENT-BENCHMARK-SETUP.md).
The same defaults without `--pair-conversion` keep the former
operator selection. The quantum run and retained shell checkpoint have
been reproduced with the final logging enabled; see also the
[quantum bridge and streaming control](INDEPENDENT-LOW-REHEAT-QUANTUM-BRIDGE-2026-11-26.md).

```sh
export OPENBLAS_NUM_THREADS=1
run=(python -W error joint_resonant_frw_split_reference.py EOS.late.dat --eos-lower-knot-limit .001 --physical-amplitude --empty-pulse .00002 --empty-reheat .0090190699 --visible-source --visible-backreaction --support-low .8 --empty-cells 18 --empty-stop .35 --trace-exterior --pair-conversion)
"${run[@]}" --source-only --matched-vacuum-order --steps 12 24 48 --save-snapshots early_shape.json > early_vacuum.json
"${run[@]}" --quantum-visible-source --with-elastic --vacuum-scaled-residual --allow-empty-shells --steps 12 24 > early_quantum.json
python -W error scalar_pair_expanding_continuation_reference.py EOS.late.dat early_shape.json --stop 10 --steps 16 32 > bridge.json
python -W error scalar_pair_expanding_continuation_reference.py EOS.late.dat early_shape.json --stop 4 --pad 64 --steps 8 16
python -W error early_scalar_pair_budget_reference.py EOS.late.dat results/low-reheat-matched-vacuum-2026-11-26.json early_vacuum.json early_quantum.json --late bridge.json
```
