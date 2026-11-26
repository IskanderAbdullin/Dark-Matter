# Compare explicit q64 support over three already archived q32 intervals

Extend the [earlier coupled comparison](INDEPENDENT-WARM-MATCHED-EXPLICIT-SUPPORT-COUPLED-2026-11-26.md)
from the actual first q32 snapshot at 10.391 MeV to the already
archived q32 endpoint at 10.238 MeV. Replay the same three
expansion intervals on q64 with initialized vacuum outside q32,
the same inherited collisionless gravity, and the original
prompt, selected elastic and signed scalar pair operators.

The q64 run accepts all three blocks, retrying the nonlinear
initialization once in each of the last two. The terminal
relative expansion is about (1.014963). Relative to initial
species numbers, the common-shell (L^1) count differences
are (4.89	imes10^{-7}) for scalars and (8.36	imes10^{-8})
for the light fermion. Newly represented q64 shell count
fractions outside q32 are (1.13	imes10^{-11}) and
(2.19	imes10^{-11}). Compare traced positive inherited
cohorts with their original values under the exact expected
rebase, and match the two temperatures and scale factors.
Combined number accounts for new collisionless receipts.

As a separate frozen rate diagnostic at the three actual q64
checkpoints, measure the bare high-parent inverse source into
(q_Hge1) GeV in the comoving anchor frame. The local
energy rate divided by current retained density times Hubble
decreases from (3.79	imes10^{-32}) to (3.67	imes10^{-32}).
This is only the represented (L+phi	o H) source, not
an upper on other dark or visible high-heavy channels.

The comparison accumulates a measured support difference
on the archived time grid. It is not a full time-convergence
claim at this longer interval and does not recover older
unknown phase distributions or predict the late dark abundance.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_extended_coupled_support_archive_audit_reference.py \
 EOS.dat matched-q64 --execute
~~~

[Diagnostics](results/warm-extended-coupled-support-2026-11-26.json).
Exact q64 checkpoints: [block 1](results/warm-extended-coupled-q64-1-blocks-2026-11-26.json.gz),
[block 2](results/warm-extended-coupled-q64-2-blocks-2026-11-26.json.gz),
[block 3](results/warm-extended-coupled-q64-3-blocks-2026-11-26.json.gz).
