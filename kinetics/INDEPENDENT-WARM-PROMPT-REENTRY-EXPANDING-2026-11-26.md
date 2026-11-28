# Short expanding prompt control with an explicit outer reservoir

The previous [tail comparison](INDEPENDENT-WARM-PROMPT-REENTRY-SHAPE-AUDIT-2026-11-26.md)
held the clock and physical momenta fixed. Start now from
the actual retained q32 checkpoint and compare a q32 closed
mesh with q64 explicit, initially empty extra shells. Both
systems undergo the same prescribed visible-entropy cooling
from the same radiation table. Retained dark energy contributes
to the collision Hubble normalization. The cold and outer shell
edges redshift on fixed comoving grids.

Split each log-expansion block into two conservative redshifts
around an implicit prompt step. Compare four and eight blocks
from (b=1) to (1.0001); final visible temperature is
10.389474 MeV. The comoving number of each species is conserved
to (5.6	imes10^{-14}) relative. Energy after including the
integrated expansion work closes to (5	imes10^{-15}).
The q64 additional shells at the fine endpoint contain
(1.144	imes10^{-11}) and (2.215	imes10^{-11}) of the
scalar and light counts. Their respective species energy
fractions are (1.238	imes10^{-10}) and (1.732	imes10^{-10}).

The absolute shared-shell number (L^1) change, normalized to
each species' initial count, is (4.44	imes10^{-9}) (scalar)
and (4.10	imes10^{-10}) (light) between q32 and q64 on the
fine time grid. The estimated time-discretization changes
from four to eight blocks are larger, about (6.0	imes10^{-6})
and (4.0	imes10^{-6}); do not read the smaller support
difference as a continuum precision claim.

Only the prompt operator acts in this short expanding check.
There is no new injection, elastic scattering, signed scalar
pair conversion, older inherited exterior phase distribution,
or physical heavy evolution. Past omitted exterior is not assigned
a distribution in a newly added shell.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_prompt_reentry_expanding_reference.py EOS.dat
~~~

The [radiation table](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat)
digest is retained in the
[diagnostics](results/warm-prompt-reentry-expanding-2026-11-26.json.gz).

## Longer prompt-only control

Run the same matched four/eight-block setup to (b=1.01)
((T=10.287946) MeV). The eight-block q64 population
beyond the q32 edge is (1.0984	imes10^{-11}) scalar and
(2.1493	imes10^{-11}) light, in species number fractions.
The common-shell support differences at fine time resolution
are (4.50	imes10^{-7}) and (7.41	imes10^{-8}), smaller
than the step refinement differences (1.02	imes10^{-4})
and (9.39	imes10^{-5}). Each q32 block accepts on its
first nonlinear initialization, whereas the fine q64 run
records seven rejected initializations with all eight blocks
ultimately accepted. This measures solver effort for these
controls, not a verified resource-limit failure. Comoving
counts and integrated expansion work still close.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_prompt_reentry_expanding_reference.py EOS.dat --b 1.01
~~~

[Longer diagnostics](results/warm-prompt-reentry-expanding-b101-2026-11-26.json.gz).
