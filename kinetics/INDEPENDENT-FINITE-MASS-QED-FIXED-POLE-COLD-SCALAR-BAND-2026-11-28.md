# Cold scalar number and scalar energy at one fixed transverse photon momentum

The [two-scalar-mode control](INDEPENDENT-FINITE-MASS-QED-TRANSVERSE-TWO-SCALAR-MODES-2026-11-28.md) shows that a single isotropic transverse photon occupation can drive opposite changes at two resolved q. Here integrate the actual piecewise-constant scalar occupancy in the declared cold interval [1e−6,1.5e−5] at k=.0015. This is a finite retained-band calculation, not the entire scalar spectrum or a momentum-integrated photon solution.

All three outer quadratures split at the actual scalar shell edges, the extreme-cosine openings of either charged histogram, and the boundary between timelike and spacelike charged transfer. Two independent bracketing grids check the scalar q knots. At each q, the full cosine integral includes the positive residue in ω and separate inverse and forward Pauli inputs. Source coefficients are multiplied by q²/E_s for number and q² for scalar energy. A common normalized factor k²/(64π⁵) times g² y² and dk is not folded into either set of coefficients. This uses a single charge orientation and the isotropic scalar phase-space measure.

The scalar number and scalar energy each have the form n D−B, with both coefficients integrated independently and n the common transverse pole occupation at this fixed k. In addition to the full scalar energy, accumulate a separate centered moment using the exact excess
`E_s−m_s = q²/(E_s+m_s)`.
Then the relative difference between energy-neutral and number-neutral n is computed without subtracting two almost equal thresholds:
```text
n_energy/n_number − 1 = (B_excess/B_num − D_excess/D_num) /
                         (m_s + D_excess/D_num).
```
The ordinary full energy moment is separately reconstructed from the number and centered moments. Compare the change of the **paired threshold difference** under refinement; the ratios are much more correlated than the two extensive coefficients.

| Archive | Number-neutral n | Scalar-energy-neutral n | Relative shift | Last two-mesh change in shift |
| --- | ---: | ---: | ---: | ---: |
| fine | 2.9272790170228e-09 | 2.927279195386e-09 | 6.09314083e-08 | 6.95e-20 |
| third | 2.3756478716318e-09 | 2.3756480023326e-09 | 5.5016916e-08 | 4.25e-20 |

The centered coefficient changed by up to 5.36e−5 between the coarse and intermediate meshes, motivating the third quadrature. The last two meshes agree to 1.78e−13 on all retained coefficients.

In this isolated band, the trial number-neutral occupation leaves a small net loss of scalar energy, because the scalar energy-neutral occupation is slightly higher. Photon heat is neutral at the scalar number-neutral occupation here, since every event has the same fixed photon ω. This differs from the [prior varying-k longitudinal comparison](INDEPENDENT-FINITE-MASS-QED-TIMELIKE-POLE-MOMENTUM-BAND-2026-11-28.md), where varying photon energy rather than varying scalar energy separated the balances. Neither comparison fixes the actual n or a matched energy flow.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p cold-work
for g in fine third; do
  for m in 0 1 2; do
    python -W error finite_mass_qed_fixed_pole_scalar_q_band_reference.py \
      "$g" "$m" cold "cold-work/$g-cold-center-$m.jsonl"
  done
done
python -W error finite_mass_qed_fixed_pole_scalar_q_band_receipt.py cold-work
python -W error finite_mass_qed_fixed_pole_scalar_q_band_receipt.py results --gzip-input
```
The [summary](results/finite-mass-qed-fixed-pole-cold-scalar-band-2026-11-28.json) is reproduced from six compressed panel receipts. No scalar shoulders outside the named interval, real-time mode population or full reaction rate is inferred.
