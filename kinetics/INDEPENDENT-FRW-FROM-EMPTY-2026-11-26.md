# Finite-support interacting FRW evolution starting from an empty dark sector

This independent numerical control replaces the collisionless 45 MeV seeded snapshot in the [continuing-source calculation](INDEPENDENT-VISIBLE-SOURCE-CONTINUING-FRW-2026-11-26.md) with an actual empty dark sector at visible reheating temperature 50 MeV. The visible background follows the given entropy-conserving equation of state. Evolve to 45 MeV at the physical benchmark rates.

## Empty start

A finite entropy-variable parametrization cannot represent exactly zero occupation. Start without any dark particles and integrate the visible production rates over one short first interval in log(a), freezing rates and the visible-only Hubble rate at its midpoint, and depositing number and kinetic-energy moments at the interval endpoint. Ignore redshift and collisions in this first pulse. Reconstruct separate FD/BE states using only the resulting positive retained moments. There is no imposed occupation floor. A retained shell with unresolvable first-pulse source is rejected rather than silently initialized. Both tails are excluded from dynamics and are reported as separate channel ledgers for the first pulse and later production.

## Result

| q/T_R interval and physical channels | Midpoint intervals after first pulse | Final scalar share of retained dark energy | Scalar collisional heat / total retained injection in the a⁴ energy ledger |
| --- | ---: | ---: | ---: |
| [0.15,16], prompt and below-pole | 12 | 0.300103 | 0.077741 |
| [0.15,16], prompt and below-pole | 24 | 0.301519 | 0.079158 |
| [0.05,16], prompt and below-pole | 12 | 0.299950 | 0.077604 |
| [0.05,16], prompt and below-pole | 24 | 0.301363 | 0.079018 |
| [0.05,16], only source and redshift | 24 | 0.222438 | 0 |
| [0.05,16], with prompt, without below-pole | 24 | 0.301366 | 0.079021 |

At [0.15,16], changing the first-pulse interval from 1e-10 to 1e-6 changes the 24-step scalar energy share by 1.23e-8. At [0.05,16], both projected species counts and the energy balance including redshift, production and equal-opposite collisional heat close to 4.3e-14 relative at 24 steps. The change in the final scalar share between 12 and 24 steps exceeds the change from lowering the cut and rearranging the fixed 18-cell mesh. No separate momentum-resolution bound is inferred from that comparison.

On the same [0.05,16] mesh, switching off elastic scattering changes the 24-step scalar energy fraction by only 2.8e-6. This is a bulk statistic on a finite support and cannot say how fast unrepresented soft modes fill.

These are finite-support moment results, not a relic density or a demonstration of thermal equilibrium. The same vacuum dark final-state approximation for visible injection and no inverse visible reactions are used. The heavy distribution is not evolved and decays are taken as immediate and isotropic. The visible EOS is imposed and not evolved in response to dark energy withdrawal. The very cold bins are excluded: even though elastic scattering can populate them, direct production of a very soft scalar from a boosted parent is strongly suppressed. For scalar lab momentum p, the initial parent requires E greater than or equal to M(Es²+p²)/(Es sqrt(p²+mphi²)+P p). At 50 MeV that lower bound is roughly 96 T_R at p/T_R=0.01, and 192 T_R at 0.005. Runs with those smaller momentum cuts have not passed the centered implicit collision solve and are not reported as viable evolutions.

## Reproduce

Run from `kinetics` with the same external central EOS as earlier checks; the table is not redistributed. The empty-start runs do not read the old 45 MeV snapshot.

```sh
export OPENBLAS_NUM_THREADS=1
python -W error joint_resonant_frw_split_reference.py EOS.dat --physical-amplitude --visible-source --with-elastic --empty-pulse 1e-10 --support-low .15 --steps 12 24
python -W error joint_resonant_frw_split_reference.py EOS.dat --physical-amplitude --visible-source --with-elastic --empty-pulse 1e-6 --support-low .15 --steps 12 24
python -W error joint_resonant_frw_split_reference.py EOS.dat --physical-amplitude --visible-source --with-elastic --empty-pulse 1e-9 --support-low .05 --steps 12 24
python -W error joint_resonant_frw_split_reference.py EOS.dat --physical-amplitude --visible-source --source-only --empty-pulse 1e-9 --support-low .05 --steps 12 24
python -W error joint_resonant_frw_split_reference.py EOS.dat --physical-amplitude --visible-source --empty-pulse 1e-9 --support-low .05 --steps 24
```
