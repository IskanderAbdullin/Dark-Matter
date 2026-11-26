# Independently continue q32 and q64 from 9.748 to 9.027 MeV

Start both supports from the respective accepted
[matched 9.748 MeV checkpoints](INDEPENDENT-WARM-Q64-HISTORY-SERIAL-2026-11-26.md).
Each support uses only its own past shapes to initialize the
next nonlinear solve. Evolve the same 16 intervals to
relative expansion 1.08, selected elastic, signed pair and
prompt interactions, and the inherited collisionless gravity
of each branch. No narrower solution is used in q64.

Both branches accept all 16 steps using their own history.
The endpoint temperature is 9.026983 MeV. The common-shell
q32/q64 count (L^1) differences are (1.787	imes10^{-6})
and (4.333	imes10^{-7}), normalized to the two initial
q32 counts. Full energy changes relative to initial total
dark energy are (2.72	imes10^{-7}) and (3.43	imes10^{-7}),
and the corresponding pressure changes (9.05	imes10^{-8})
and (8.09	imes10^{-8}). The added q64 region contains
(9.06	imes10^{-12}/2.07	imes10^{-11}) of its terminal
scalar/light counts. The latest represented light peak is
about (5.123	imes10^{-5}) per spin, a reconstruction
diagnostic, not a proof about unrepresented modes.

The independent audit checks each checkpoint parent link,
EOS and seed identity, rebases stored exterior gravity and
verifies the combined retained and new-collisionless number.
The terminal q64 bare high inverse rate, choosing the original
10.391 MeV comoving reference, is (2.69	imes10^{-32}) of
current retained density times Hubble. This is a frozen
(L+phi	o H) term only. Neither other high-heavy inputs
nor a physical evolved heavy distribution is included.

This is a same-time-grid support comparison; both grids can
have larger time-discretization error. Missing visible and
untracked exterior physics still preclude a final abundance.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p matched
STEP=$(python -c 'print(1.08**(1./16))')
for q in 32 64; do
 zcat results/warm-history-lower-seed-q$q-2026-11-26.json.gz > matched/seed$q.json
 python -W error warm_q64_history_serial_reference.py EOS.dat matched/seed$q.json matched/q$q \
  --step-expansion "$STEP" --blocks 16
done
python -W error warm_history_lower_matched_audit_reference.py EOS.dat matched
~~~

[Matched audit](results/warm-history-lower-matched-16-2026-11-26.json.gz).
Exact terminal pairs: [q32](results/warm-history-lower-q32-16-terminal-2026-11-26.json.gz)
and [q64](results/warm-history-lower-q64-16-terminal-2026-11-26.json.gz).
