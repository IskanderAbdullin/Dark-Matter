# Scalar number and energy in the first resolved future-transfer band

Continue the [positive-energy spacelike transfer test](INDEPENDENT-FINITE-MASS-QED-FUTURE-SPACELIKE-POLE-2026-11-28.md) by integrating q in [.004,.015] at the same frozen transverse pole k=.0015. Use actual Pauli factors, the scalar occupation on each archived cell, and exact angular and radial partitions. The scalar number, scalar energy, and energy above rest mass are accumulated separately. The omitted scalar momenta between the cold band and this band are not interpolated or set to zero.

There are 22 fine and 30 third exact scalar panels. The first mesh pair failed the tight coefficient resolution test (drift .00560 and .00444). The accepted intermediate-to-high maximum changes are 3.56e−6 and 1.67e−6. The raw forward coefficient for this band is larger than its archive's separate [cold band](INDEPENDENT-FINITE-MASS-QED-FIXED-POLE-COLD-SCALAR-BAND-2026-11-28.md) by 9.40e6 and 7.09e7, for fine and third. These ratios compare positive coefficients at this **one pole**, not total dark-sector production, physical photon populations, or the net contribution when reverse and forward terms balance.

| Archive | Conditional number-neutral n | Conditional scalar-energy-neutral n | Mean forward / inverse scalar energy |
| --- | ---: | ---: | ---: |
| Fine | 6.7835248e−5 | 4.8430344e−5 | .01174030 / .00838204 |
| Third | 6.5807624e−5 | 4.6748933e−5 | .01173862 / .00833913 |

Between these two conditional neutral populations on either grid, scalar number decreases while scalar energy increases. Photon energy increases because it is minus the scalar net number times the fixed pole energy; the charged sector loses the sum of the two energy gains. The receipt independently reconstructs that charged transfer from the extensive forward/inverse ledger. This is a possible redistribution at a trial n, not a determination of n or a cosmological cooling prediction.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error finite_mass_qed_future_scalar_band_receipt_pilot.py --meshes 1 2
# Regenerate on explicit work paths:
for g in fine third; do
  for m in 0 1 2; do
    python -W error finite_mass_qed_future_scalar_band_pilot.py "$g" "$m" first_future "$g-first-future-$m.jsonl"
  done
done
```
Compare the [result](results/finite-mass-qed-first-future-scalar-band-2026-11-28.json) to the archived panel logs. No medium continuum, matched virtual term, or other momenta are included; the common g²y² k² dk/(64π⁵) is omitted.
