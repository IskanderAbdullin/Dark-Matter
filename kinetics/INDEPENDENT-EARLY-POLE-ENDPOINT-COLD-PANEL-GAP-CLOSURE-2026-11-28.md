# Check the scalar gap under the same frozen finite-panel reweights

The [cold source envelope](INDEPENDENT-EARLY-POLE-ENDPOINT-COLD-PANEL-REWEIGHT-ENVELOPE-2026-11-28.md) holds the fitted scalar gap fixed. Here use the identical S2 and S4 q-panel occupation degrees of freedom and count/optional mean excess constraints, but solve the frozen Hartree equation for the gap itself. The excess constraints are evaluated at the original nodal reference mass and remain fixed. The prescribed, positive out-of-band step shape is always retained.

For each evaluation at a trial mass (M>0), a positive shape gives (I(M)=int dN/sqrt{q^2+M^2}). Independently optimize that inverse-energy integral over the finite occupied panels while preserving their declared moments. Every allowed (I) decreases with (M), and the maximized and minimized envelopes are continuous and non-increasing as well. Accordingly, their roots of (M^2-m_mathrm{bare}^2-lambda I(M)/2) delimit the mass of the finite-panel feasible family. Compare the baseline fixed shape and stored nodal gap independently; the root-procedure transports no reaction coefficient through that change.

| Chart | Step-reconstructed baseline gap minus stored nodal, relative | Width at fixed interval counts | Width at fixed counts plus mean excess |
| --- | ---: | ---: | ---: |
| S2 | −6.67e−8 | 6.11e−7 | 2.12e−13 |
| S4 | −1.70e−8 | 1.68e−7 | about 9e−15 |

The last S4 interval is near ordinary floating-point resolution; 24/48 positive rules gave the same values in this run, but this is not directed-rounding certification. The strong reduction also does not show that the moments themselves are accurate. Most importantly, this is an alternate *frozen finite closure* at fixed charge populations, not actual cosmological mass uncertainty, a self-consistent time step, or a bound on how varying that mass changes the cold collision kernel.

Reproduction from `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28
for s in 2 4; do
 python -W error early_pole_endpoint_cold_panel_gap_closure.py \
  results/fine-charge2-scalar$s-endpoint-continuum-cold-{1,2}.jsonl.gz \
  results/early-pole-fine-charge2-scalar$s-x05-state-$d.json.gz \
  results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json --scalar-div "$s"
done
```
