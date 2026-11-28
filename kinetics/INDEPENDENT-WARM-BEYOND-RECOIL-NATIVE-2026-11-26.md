# Complete the native support baseline beyond the recoil cutoff

Independently continue the exact q32 and q64 refined recoil
endpoints from 8.133100 to 7.394049 MeV, with 24 equal intervals
in ln(a) and relative expansion 1.10. This baseline uses the
original pair quadrature; it deliberately leaves the subsequently
diagnosed reverse-projection defect in place. Both trajectories
use their own accepted past for numerical seeds. The ordinary
fallback is accepted in one q32 and five q64 blocks.

The q32 first-cell native prompt rate is geometrically zero
throughout. The q64 rate is nonzero and negative, ending at
−6.03e-9 per Hubble time in first-cell occupation units.
Both **native** trajectories exhibit their only cold-cell rise
at step 12. The preceding common-chart and cold-lift
replays independently identify it as a pair projection artifact,
not a q32-specific production or recoil feature.

The shared-shell number L1 at the endpoint is 3.49e-6 for
scalar and 1.16e-6 for light, relative to the original q32 species
inventory. At the common start those were already 3.45e-6
and 8.89e-7. The endpoint full-energy L1 values normalized to
initial total dark energy are 6.65e-7 and 8.38e-7. These comparisons
include inherited initial differences and original-quadrature
effects, and are not measurements of support error alone.
Inherited positive exterior cohorts enter gravity and do not
supply old incoming phase occupation. The maximum cumulative
total number defect including new exterior births is below 2e-13
on each branch.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p original
for q in 32 64; do
  zcat "results/warm-recoil-crossing-q$q-terminal-2026-11-26.json.gz" > "original/seed$q.json"
  s=$(python -c 'print(1.10**(1/24))')
  python -W error warm_q64_history_serial_reference.py EOS.dat "original/seed$q.json" "original/q$q" --step-expansion "$s" --blocks 24
done
python -W error warm_beyond_recoil_matched_reference.py EOS.dat original 24
~~~

[Complete native support baseline](results/warm-beyond-recoil-native-matched-2026-11-26.json.gz);
[q32 endpoint](results/warm-multishell-trajectory-native24-q32-2026-11-26.json.gz);
[q64 step 23](results/warm-beyond-recoil-native-q64-step23-2026-11-26.json.gz)
and [q64 endpoint](results/warm-beyond-recoil-native-q64-step24-2026-11-26.json.gz).
