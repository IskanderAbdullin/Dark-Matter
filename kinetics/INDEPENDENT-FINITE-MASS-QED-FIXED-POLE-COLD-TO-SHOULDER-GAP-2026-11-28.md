# Closing the small scalar-momentum gap at the frozen transverse pole

The [cold and shoulder comparison](INDEPENDENT-FINITE-MASS-QED-FIXED-POLE-SCALAR-SHOULDER-2026-11-28.md) retained two disjoint bands. Here integrate exactly the missing q interval [1.5e−5,4e−5] at the same fixed k=.0015. The union is now contiguous from q=1e−6 to 3e−4, but still omits both tails. Its boundaries are taken from the independent [algebraic knot audit](INDEPENDENT-FINITE-MASS-QED-SCALAR-Q-KNOT-EXACT-2026-11-28.md), which agrees with the existing bracket construction.

The two gap meshes agree on all retained coefficients within 1.74e−13 and 7.57e−14 (fine and third); their common support contains 23 and 25 panels. At the trial occupation n=1e−8 the gap number source is positive but small, +4.42e−21 and +1.60e−21 before the common rate factor. On the same trial, the shoulder loses number and the cold band gains it. In particular the small omitted gap cannot explain the measured two-archive difference.

| Archive | Gap number-neutral n | Joined number-neutral n | Joined scalar-energy-neutral n |
| --- | ---: | ---: | ---: |
| Fine | 4.9504753e−9 | 1.5791222e−8 | 1.6233979e−8 |
| Third | 4.4026414e−9 | 4.0805173e−8 | 4.1382184e−8 |

Adding the gap shifts either joined neutral occupation by at most 1.20e−6 relative to the two-band result. Number and scalar-energy balances still differ, and each cancels only after summing signed fluxes across the retained support. The joined number-neutral occupations still differ by about 2.58 between archived inputs. The fixed-pole test does not determine the actual n, supply an uncertainty bound for q outside the stated union, or include other collision terms.

Replay from kinetics:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error finite_mass_qed_gap_band_merge_pilot.py
# To regenerate the four new raw gap meshes:
mkdir -p gap-replay
for g in fine third; do
  for m in 1 2; do
    python -W error finite_mass_qed_fixed_pole_scalar_q_band_reference.py \
      "$g" "$m" gap "gap-replay/$g-gap-center-$m.jsonl"
  done
done
python -W error finite_mass_qed_gap_band_merge_pilot.py --gap gap-replay --raw-gap
```
The [receipt](results/finite-mass-qed-fixed-pole-cold-to-shoulder-gap-2026-11-28.json) reads the earlier cold and shoulder panels together with four new gap panel logs. All rates omit the common g²y² k² dk/(64π⁵).
