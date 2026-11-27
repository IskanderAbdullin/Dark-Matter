# Admit heavy regeneration in the cold weighted-count ceiling

The [initial reheating count bound](INDEPENDENT-WARM-DIRECT-RESTRICTED-COUNT-UPPER-2026-11-26.md)
omitted heavy regeneration while comparing unmodeled
earlier rates. This is a different conditional test:
start at the [modeled 2 MeV or 80 keV anchor](INDEPENDENT-WARM-REHEATING-POST-2MEV-ELASTIC-2026-11-26.md)
and allow the selected vacuum regeneration kernels
after that point. Do not continue a two-temperature
shape by assumption.

The count is Q=N_scalar+N_light+2N_heavy. Prompt
decays and their inverse, elastic redistribution, and
stable-pair conversion to scalars or back preserve Q.
Dark heavy regeneration increases it by two per
pair. Use the [independent all-momentum energy-weighted
kernel](INDEPENDENT-WARM-FROM-EMPTY-RECYCLING-ENVELOPE-2026-11-26.md)
to upper-bound the simultaneous positive count source
by K ρ_dark²/M_heavy². This coefficient already counts
both parents per event and omits final heavy
Pauli blocking in the conservative direction.
Drops from reverse-heavy events are omitted.
Do not insert a favored scalar momentum shape.

If the sector receives no more energy or particles,
nonnegative pressure implies ρ_dark(b)b³≤ρ_dark(1),
with b normalized to the anchor. The photon floor
H≥H_gamma/b² then gives a comoving count increase
at most K ρ_dark(1)²/(M_heavy² H_gamma), or that
quantity divided by the initial visible entropy
for the extra yield. Number-conserving momentum
redistribution need not be weak for this inequality.

| Modeled visible anchor | Starting weighted yield | Extra yield upper, even allowing heavy regeneration | Scalar Ωh² upper if all remaining Q becomes scalars |
|---:|---:|---:|---:|
| 2 MeV | 6.59290e−4 | 5.37306e−5 | .016920 |
| 80 keV | 6.59290e−4 | 6.96620e−8 | .015647 |

The lower-temperature result is about 26.1% of the
.060 value used for the half-dark target. This is
an upper within the isolated cold future and the
stated weighted-count network, not a final yield
or a claim about unmodeled sources. In particular
the bound does not assume the actual anchor was
reached in the full kinetic theory. No external
lepton injection/heating or extra weighted-count
process is included in this upper, unlike the
separate lower bound granting an optional vacuum
leptonic annihilation loss.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_post_anchor_scalar_ceiling_reference.py EOS.late.dat late.json \
 results/warm-energy-weighted-recycling-kernels-2026-11-26.json
~~~

[Energy, count and normalization ledgers](results/warm-prompt-reheating-post-anchor-scalar-ceiling-2026-11-26.json).

## Allow the specified vacuum charged-lepton energy input

Use the [separate vacuum late lepton comparison](INDEPENDENT-WARM-LATE-LEPTON-RECYCLING-2026-11-26.md)
to enclose pair input by C_lepton T_visible^9
and fractional elastic heating by κ T_visible^5.
With entropy degrees at least two, take
T_visible≤A T_anchor/b where A=(g_s,anchor/2)^(1/3).
The photon expansion floor then bounds the
energy integral by C_lepton(A T_anchor)^9/(4H_gamma)
and the heating exponent by κ(A T_anchor)^5/(3H_gamma).
Grant the entire positive input from the start
of the subsequent regeneration test. Each produced
count unit costs at least m_light in the retained
processes, including two fast daughters per heavy
parent. Include both the direct count and the extra
regeneration due to this energy.

This increases the .016920 isolated result to
.016938 at 2 MeV. At 80 keV the upper remains
.01564651 at the shown precision. The extension is
still restricted to the stated vacuum lepton and
dark channels, not a bound for arbitrary external heating.

~~~sh
python -W error warm_prompt_reheating_cold_lepton_energy_allowance_reference.py EOS.late.dat \
 results/warm-prompt-reheating-post-anchor-scalar-ceiling-2026-11-26.json
~~~

[Lepton energy allowance and updated count upper](results/warm-prompt-reheating-cold-lepton-energy-allowance-2026-11-26.json).
