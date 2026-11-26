# Retain outgoing pair scalars from 3.43 MeV onward

The [previous widened test](INDEPENDENT-WARM-PAIR-SUPPORT-CLOSURE-2026-11-26.md)
began at 1.206 MeV. Start earlier, from the retained 32-step
3.433865 MeV snapshot of the 50 MeV history, and continue the signed
pair operator alone until 0.413763 MeV. This time range includes
the period when even unpadded retained scalars can make a stable pair.
Compare q/T_R ≤ 16 with an otherwise identical initially empty
extension to q/T_R ≤ 96.

| Steps | Support | Inclusive stable number / same initial retained stable number | Peak reverse/forward | New exterior scalar count / initial retained scalar count |
| ---: | ---: | ---: | ---: | ---: |
| 36 | 16 | 0.935684134 | 0.021399 | 0.0082446 |
| 36 | 96 | 0.935706550 | 0.021524 | 1.17e-27 |
| 52 | 16 | 0.935660909 | 0.021918 | 0.0082428 |
| 52 | 96 | 0.935683169 | 0.022006 | 3.94e-26 |

Retaining the new tails raises the final stable count by about
2.4e-5 fractionally in either paired comparison. Changing step
count shifts either standalone result by a similar amount.
The energy ledger is within 3.1e-14. The extra cells strongly
change how the newly produced scalars are partitioned between
retained and exterior populations, but have little effect on
the inclusive stable count in this test.

This is a pair-only test on a prescribed radiation clock.
The prompt and elastic dynamics, dark energy contribution to
the clock, and exterior cohorts produced before 3.43 MeV are
not included here. The warm coupled benchmark is not obtained
by substituting this pair-only result. Finer angular and
momentum quadrature are independent from the time-step control.

~~~sh
export OPENBLAS_NUM_THREADS=1
r=results/warm-coupled-pair-b35-32-shape-2026-11-26.json
for n in 36 52; do
  python -W error scalar_pair_expanding_continuation_reference.py EOS.dat "$r" \
    --stop 8.571428571428571 --steps "$n" --quadrature 3 7 6 2
  python -W error scalar_pair_expanding_continuation_reference.py EOS.dat "$r" \
    --stop 8.571428571428571 --steps "$n" --quadrature 3 7 6 2 --pad 96
done
~~~
