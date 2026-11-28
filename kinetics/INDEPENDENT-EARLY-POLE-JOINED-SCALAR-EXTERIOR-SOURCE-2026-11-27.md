# Frozen exterior scalar emission carries more energy than its number fraction suggests

Evaluate the reciprocal decay/inverse source outside the live scalar window, without changing the [completed five-channel hard trajectories](INDEPENDENT-EARLY-POLE-FINE-EXTENDED-MEMORY-2026-11-27.md). Replay their instantaneous fitted hard distributions, fixed hard-leg support and actual scalar Hartree pole. Use live occupations in q=[1e−6,1e−5] and the original frozen archived occupations outside it. The outgoing momentum test imposes both triangle inequalities, including q ≤ p_parent+p_daughter at higher scalar momentum. No supported packet on this hard box can reach q=5. The lower numerical cut is 1e−8. The independently integrated kernel reproduces the saved instantaneous cold source mode sums at the original eight nodes within 8.1e−13 per mode.

| Fine archive, initial hard reconstruction | Initial exterior/cold net count source | Final exterior/cold net count source | Final exterior/cold scalar energy source |
| --- | ---: | ---: | ---: |
| Native steps | 4.494e−6 | 3.604e−6 | 0.2633% |
| Count-conserving continuous | 6.137e−6 | 3.619e−6 | 0.2659% |

The ratios are to the cold-window source, not to the number or energy already present. Most of the exterior energy source is in q=[.1,5], and spontaneous creation dominates that band's gross balance. At the final step in the native case, the q≥.1 band gives a net number rate 2.614e−7 and scalar energy rate 2.690e−7 in the reference normalization. Doubling the outer quadrature order and increasing the inner order moves the quoted endpoint exterior ratios by less than 0.05% relative.

Integrating all nine native snapshots with two time interpolants gives exterior/cold number-source integrals 3.78–3.79e−6 and exterior/cold energy-source integrals 0.287–0.288%. This is a frozen-source diagnostic: the hard legs do not receive the missing reciprocal work, and the exterior occupations never evolve. Even the interpolated cold number gain differs from the true solver's live cold count gain by about 0.14% for the higher-order interpolant. None of these ratios is a physical hot-tail abundance, a comparison to an observational constraint, or a bound including the unmodeled hard exterior.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
b=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
python -W error early_pole_joined_scalar_exterior_source_reference.py "$b" --q-order 4 --radial-order 8
python -W error early_pole_joined_scalar_exterior_source_reference.py "$b" --q-order 8 --radial-order 12
python -W error early_pole_joined_scalar_exterior_source_reference.py "$b" --q-order 8 --radial-order 12 --projection cell-conservative
python -W error early_pole_joined_scalar_exterior_source_reference.py "$b" \
    --q-order 6 --radial-order 10 --snapshots all > exterior.json
python -W error early_pole_joined_scalar_exterior_history_reference.py exterior.json
```
Raw: `results/early-pole-joined-fine-216-scalar-exterior-2026-11-27.json.gz`.
