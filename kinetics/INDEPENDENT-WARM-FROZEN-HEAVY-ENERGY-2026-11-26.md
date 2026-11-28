# How much heavy energy is in the frozen decay balance?

The [direct light source test](INDEPENDENT-WARM-LIGHT-DIRECT-PAIR-DOSE-2026-11-26.md)
needs a time-uniform heavy energy constraint much tighter than the
total dark budget. As an independent diagnostic, freeze the
represented light/φ shapes and integrate the heavy stationary
distribution from decay plus inverse decay alone.
This freezes the daughters and omits every other heavy source.
It is **not** an upper on the coupled heavy population.

The occupied daughter energies cap the heavy support by
(E_H≤E_{light,max}+E_{φ,max}). Integrate the exact
angle-split stationary distribution as an unpolarized
two-charge gas with four states. Check three radial
orders on 69 panels. Relative to the actual retained
light/φ energy of each input, the q16/q32 10.442 MeV
fractions are (1.42517×10⁻⁶) and (1.42900×10⁻⁶).
At 3.434 MeV the number is (4.01143×10⁻¹⁴).
The 1.206 and .414 MeV resolved daughter supports
cannot reach the heavy energy threshold, so this
stationary source is identically zero in those
representations. There is no claim about omitted
daughter modes, transient histories or external sources.

This checks a useful numerical scale, not the 4.83%
uniform-budget condition for the partial direct light
comparison. The dynamical heavy fraction and the
remaining decay/elastic light-occupation gains still
need independent checks.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_frozen_heavy_energy_reference.py \
 "$D"/warm-early-q{16,32}-time32-checkpoint-2026-11-26.json.gz \
 "$D"/warm-coupled-pair-b{35-48,10-42,30-32}-shape-2026-11-26.json
~~~
