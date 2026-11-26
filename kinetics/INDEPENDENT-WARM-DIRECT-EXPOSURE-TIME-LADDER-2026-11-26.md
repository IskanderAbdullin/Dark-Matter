# Separate path and sample error in isolated high-inverse exposure

For the two exact moment-matched .52/.58 GeV packet seeds in the
[isolated direct stress](INDEPENDENT-WARM-DIRECT-PROMPT-STRESS-2026-11-26.md),
evolve only the direct-occupation prompt operator from Hdt=0 to 1e-5.
Hold the visible clock, volumes, and archived background fixed;
the initial H contains the specified visible and represented dark
energy. Tighten the geometric time ladder from 4 to 8 to 16 steps
per decade. Use the signed state integration but accumulate the
positive **bare** inverse energy rate in the band starting at
1 GeV. That accumulated term is not the signed net energy change,
a real expanding comoving dose, or an upper bound.

| Hard packet (GeV) | 4 / decade | 8 / decade | 16 / decade | Fraction of finest result by Hdt = 1e-7 |
| ---: | ---: | ---: | ---: | ---: |
| .52 | 3.3537e-28 | 2.9984e-28 | 2.8591e-28 | 99.09% |
| .58 | 3.1182e-28 | 2.7898e-28 | 2.6620e-28 | 92.88% |

The table integrates each reconstructed rate with trapezoidal
quadrature and gives the resulting fraction of the fixed dark energy.
For the 8→16 change, evaluate the finer trajectory on the same
coarse time nodes as well. In both cases roughly −0.99% of the
old result comes from additional integration nodes and −3.6%
from evolving the path with smaller steps. Quadrature alone
cannot remove the latter. A Simpson check in log time on each
already solved n=16 path changes its bare total by only about
0.34%; the unsampled first interval remains trapezoidal. No
asymptotic trajectory order or percent-level dose accuracy
is claimed. The .58 coarse terminal rate was still over twice
the fine one, but the two finest terminal rates differ by
less than 7e-5 relative in either case.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p dose
for h in 052 058; do
  for n in 4 8 16; do
    python -W error warm_direct_stress_exposure_reference.py EOS.dat \
      --hard ".${h#0}" --steps-per-decade "$n" >"dose/${h}_$n.json"
  done
  python -W error warm_direct_stress_exposure_ladder_audit_reference.py \
    "dose/${h}_4.json" "dose/${h}_8.json" "dose/${h}_16.json"
done
~~~

[.52 comparison](results/warm-direct-high-dose-hard052-audit-2026-11-26.json)
and [.58 comparison](results/warm-direct-high-dose-hard058-audit-2026-11-26.json).
Raw .52 trajectories: [4](results/warm-direct-high-dose-hard052-n4-2026-11-26.json.gz),
[8](results/warm-direct-high-dose-hard052-n8-2026-11-26.json.gz),
[16](results/warm-direct-high-dose-hard052-n16-2026-11-26.json.gz).
Raw .58 trajectories: [4](results/warm-direct-high-dose-hard058-n4-2026-11-26.json.gz),
[8](results/warm-direct-high-dose-hard058-n8-2026-11-26.json.gz),
[16](results/warm-direct-high-dose-hard058-n16-2026-11-26.json.gz).
