# A passive three-mode approximation for the finite-chemistry crossover

Use the [static finite-chemistry controls](INDEPENDENT-CHARGED-LEPTON-FINITE-CHEMICAL-CROSSOVER-2026-11-27.md)
to test when a compact acoustic closure preserves both frequency
and attenuation. At vanishing pair rate, the elastic collision
operator has three conserved directions (energy, momentum and pair
number). In the instantaneous fast-elastic limit, their leading
streaming and chemical decay generate the characteristic cubic

    λ³ + γ_c λ² + c_f² k² λ + γ_c c_e² k² = 0.

Here c_f and c_e are the number-frozen and chemically equilibrated
sound speeds, respectively, and γ_c=εα is the leading chemical decay
rate. For any positive ε, the small-k damping pole of the cubic has
residue (c_f²−c_e²)/(2α) in εD, exactly the residue projected from the
five-mode model. The cubic is useful for phase speed, but it omits
finite-elastic attenuation. Over the tested domain that omission can
misestimate the decay by up to 68.2%, even with a small phase error.

## Exact peak of the isolated chemical relaxation

Let r=(c_e/c_f)² and q=c_f k/γ_c, and write the real chemical
root as −a in units of γ_c. Then q²=a²(1−a)/(a−r), with r<a<1.
The squared acoustic attenuation-to-frequency ratio is N/(4r−N),
where N=(a−r)(1−a)/a. Since N=1+r−a−r/a, its maximum is at
a=√r, giving q_peak=r^(1/4). For r>1/9, N<4r so the acoustic
pair stays oscillatory for every nonzero k in this isolated cubic.

At .015 and .050 GeV the peak is at k/γ_c=1.7372 and 1.7327.
The amplitude lost in one local acoustic cycle of the *isolated
chemical cubic* there is .221925% and .0031634%, respectively.
Direct polynomial roots at and on either side of each maximum
reproduce this check. The singular long-wave coefficient as
γ_c→0 is not a singular finite loss per cycle. The numerical
peak does not bound total damping with the omitted elastic losses.

## Retaining the stationary fast-elastic response

In entropy-normalized coordinates, split the exact five-mode
generator M=−C−ikB into the three elastic null coordinates S and
their two fast complements R. At fixed retained amplitudes,
reconstruct fast amplitudes with

    fast = −(M_RR)⁻¹ M_RS slow.

The resulting three-coordinate generator is the zero-frequency Schur
complement. It keeps the full wave-number dependence of the fast
block but drops its frequency-dependent memory. It is passive:
for a reconstructed vector with vanishing fast residual, the real
quadratic form of the reduced generator equals the negative entropy
loss of that vector under the full collision matrix. Compute
attenuation from that positive loss, rather than subtracting
nearly imaginary eigenvalues.

| Pair kernel multiplier, two T controls | Max relative attenuation error, bare cubic | Max relative attenuation error, stationary fast elimination |
|---:|---:|---:|
| 10⁻³ | .682 | .0113 |
| 10⁻⁴ | .180 | 5.33e−5 |
| 10⁻⁵ | .0215 | 3.43e−7 |
| 10⁻⁶ | .00219 | 1.49e−8 |

The tests use T=.015 and .050 GeV, h_V/h_D=10, and seven wave numbers
from .01 to 100 times the slow chemical rate. The quoted errors
are maxima **within these static controls**, not error bounds for
larger k, other thermal states, or a time-dependent cosmological
closure. In particular a phase-speed match alone is insufficient
to validate attenuation.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  charged_lepton_chemical_three_mode_reference.py \
  results/charged-lepton-portal-combined-momentum-2026-11-27.json \
  results/charged-lepton-finite-chemical-sound-crossover-2026-11-27.json
```
[Full receipt](results/charged-lepton-chemical-three-mode-2026-11-27.json).
