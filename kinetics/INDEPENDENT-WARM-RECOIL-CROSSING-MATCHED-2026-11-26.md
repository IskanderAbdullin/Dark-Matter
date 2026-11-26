# Through the cold-light recoil gate on matched q32 and q64 trajectories

Earlier [closed-prompt geometry](INDEPENDENT-WARM-COLD-PROMPT-SUPPORT-GRID-2026-11-26.md)
distinguished a lost q16 event from a physical absence.
Continue both supports toward the corresponding q32 boundary.
Use sixteen matched half-size intervals from accepted q32/q64
states at 8.568461 MeV. The q64 branch uses stored and
interpolated narrow entropy profiles only as numerical
initializations; the physical targets and collision
quadratures are wider and independently solved.

A light particle at rest requires a scalar recoil energy near
.240 GeV. The redshifting q32 edge crosses that rest-point
requirement at predicted visible temperature 8.253704 MeV.
The first finite light shell remains partly connected until
8.238281 MeV because its nonzero momenta admit a range of
recoil energies. Thus a single rest-mode calculation and the
projected first cell are not interchangeable.

At 8.239831 MeV the native q32 first-cell prompt rate is
already only (-1.44	imes10^{-9}) in shell-averaged
occupation per Hubble time. At 8.213018 MeV there are zero
native prompt events in that cell. On q64 the corresponding
rates are (-1.72	imes10^{-7}) and (-1.54	imes10^{-7}).
Cropping the **same** q64 distribution to q32 reproduces
the disappearing term. Independently projecting the open
tagged-light kernel against (p^2dp) agrees with the native
q64 rate. The q32 frozen open surrogate still permits gain
from a stationary heavy parent into an unrepresented scalar;
that gain is not in the q32 update.

Both accepted trajectories still show declining first-cell
mean occupation across the next intervals, at order
(6	imes10^{-5}) per unit log expansion, when the q32
prompt term is zero. The quoted prompt difference is small
relative to this measured combined update. The result
identifies a local support artifact, not shutoff of the
full light population or a measured late abundance shift.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p crossing
zcat results/warm-recoil-crossing-q32-seed-2026-11-26.json.gz > crossing/q32-seed.json
zcat results/warm-recoil-crossing-q64-start-2026-11-26.json.gz > crossing/q64-start.json
zcat results/warm-recoil-crossing-guide-{1,2,3,4}-2026-11-26.part.gz > crossing/guides.json
STEP=$(python -c 'print(1.11**(1./32))')
python -W error warm_q64_history_serial_reference.py EOS.dat crossing/q32-seed.json crossing/q32 \
 --step-expansion "$STEP" --blocks 16
python -W error warm_q64_guided_serial_reference.py EOS.dat crossing/q64-start.json crossing/guides.json crossing/q64
python -W error warm_recoil_crossing_matched_reference.py EOS.dat crossing
~~~

[All frame diagnostics](results/warm-recoil-crossing-matched-2026-11-26.json.gz);
raw pre/post bracket [q32](results/warm-recoil-crossing-q32-bracket-2026-11-26.json.gz),
[q64](results/warm-recoil-crossing-q64-bracket-2026-11-26.json.gz);
exact final two [q32](results/warm-recoil-crossing-q32-terminal-2026-11-26.json.gz)
and [q64](results/warm-recoil-crossing-q64-terminal-2026-11-26.json.gz) checkpoints.

## First-cell channel budget after crossing

At the 8.133100 MeV endpoint, freeze each accepted state and
evaluate all three **native selected** first-cell rates.
With the evolution quadrature on q64, the prompt, selected
elastic and scalar-pair pieces are respectively
(-1.11	imes10^{-7}), (-3.75	imes10^{-5}) and
(-2.53	imes10^{-5}) per unit log expansion. Thus the
missing q32 prompt component is small in the coupled cold
cell. The instant total is close to, but different from,
the full-step change measured between accepted states.

Refining the endpoint pair and elastic angles/radial rules
leaves the component hierarchy unchanged. It shifts the
q64 signed pair estimate to (-1.79	imes10^{-5}) and the
elastic estimate to (-3.87	imes10^{-5}). The first-cell
pair rate therefore has appreciable numerical sensitivity;
these two quadratures do not supply a converged value.
This is a frozen-budget test, not a recomputation of the
trajectory with every refined operator.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_recoil_endpoint_channel_budget_reference.py EOS.dat
~~~

[Native component diagnostics](results/warm-recoil-endpoint-channel-budget-2026-11-26.json).
