# Live-parent energy ledger separates centroid error from heat exchange

This audit replays the [resolved-cell momentum redistribution tests](INDEPENDENT-EARLY-POLE-PARENT-MOMENTUM-RELAXATION-2026-11-26.md)
with a separate integral for both particle traces, the collisional
bin-centroid energy defect, and the energy transferred by the
parent equilibration operator. It uses the earlier count-conservative
cell overlap, not the later [exact projected-energy stencil](INDEPENDENT-EARLY-POLE-EXACT-ENERGY-PARENT-STENCIL-2026-11-26.md).

| Parent treatment | Coarse centroid defect | Refined centroid defect | Refined heat transfer |
|---|---:|---:|---:|
| No mixing | +0.5376% | −0.0535% | 0 |
| Fixed temperature, rate 100 | +1.8825% | +0.4824% | −0.03954 |
| Isolated mixing, rate 100 | +1.6831% | +0.4732% | < 4 × 10⁻¹⁴ absolute |

Centroid defects are signed fractions of the final daughter energy.
Heat is energy into tracked species in units of $T^4$, with
$T\propto1/a$. The negative sign is heat *leaving* the tracked
parent sector for the imposed thermostat. In the refined fixed-T run
the heat removed is 1.78 times the final daughter energy (coarse 1.84).
Thus the temperature-maintained and internally isolated controls
cannot be compared as if they had the same energy budget. The heat
removal can be a small fraction of a much larger radiation bath;
this diagnostic does not evolve that bath.

Separate parent, daughter and total projected continuity residuals
all stay below $10^{-6}$ relative to the final daughter energy.
This only verifies the ledger of each numerical scheme: the centroid
terms would be absent in the exact physical reaction. Isolated
mixing alone matches both target moments and transfers no heat.

~~~sh
# mode: none, isolated, or fixed-T
python -W error early_pole_live_parent_energy_reference.py --mode fixed-T
python -W error early_pole_live_parent_energy_reference.py --mode fixed-T --scalar 65 --fermion 90 --parent 72
~~~

[Six raw continuity ledgers](results/early-pole-live-parent-energy-ledger-2026-11-26.json).
