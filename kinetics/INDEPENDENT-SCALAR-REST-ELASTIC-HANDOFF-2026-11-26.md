# Can elastic scattering reconnect the scalar rest mode to fast decay?

The [preceding decay-only gate](INDEPENDENT-THREE-SPECIES-DECAY-IR-GATE-2026-11-26.md) cannot by itself rule out equilibration of the soft tail. With the same Yukawa coupling y=.001, evaluate φχ2 → φχ2 and its antiparticle copy when the incoming scalar is at rest. Use the exact below-pole tree-level s/u trace and final-state FD/BE factors with all other modes held in neutral thermal equilibrium. This is the local diagonal scalar relaxation coefficient, not the full projected elastic operator: an elastic event conserves scalar number and moves occupation to other bins.

At rest the incoming fermion energy determines s=m2²+mφ²+2mφ Eχ. On the truncated thermal-support quadrature s remains far below the χ1 pole, so these events do not duplicate real inverse decay. The handoff gate counts only events placing the scalar above the first local decay damping/H crossing at the same temperature. The outgoing scalar energy runs from mφ to mφ+2mφ pχ²/s. Include the detailed-balance difference factor (1−exp(−mφ/T)) rather than treating the vacuum loss frequency as the net relaxation rate. The kinematic opening of the hard outgoing gate is placed at a quadrature boundary; the scalar s/u trace also agrees with independent complex gamma-matrix multiplication to 1.2e-15.

| Common temperature (GeV) | Elastic diagonal rest damping/H | Direct handoff to decay-fast momenta/H |
|---:|---:|---:|
| .68825 | 9.45e3 | 1.34e3 |
| .30 | 1.87e4 | 2.80e2 |
| .10 | 5.50e4 | 9.60e-3 |
| .020 | 9.07e4 | Not reached on included thermal support |
| .010 | 5.28e4 | Not reached on included thermal support |

Thus direct elastic handoff is available at the hot common-bath benchmark, but can become slower than expansion while local loss or diagonal damping still appears fast. This does not prove decoupling at lower temperature: many elastic scatterings can traverse intermediate scalar momenta, and the diagonal rate cannot capture that transport. The calculation also omits other dark processes, medium corrections, and actual evolving nonthermal fermion populations. It neither detects nor excludes a condensate.

Reproduce using the same authors' central SM EOS data:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error kinetics/scalar_rest_yukawa_elastic_ir_reference.py EOS_TABLE
```
