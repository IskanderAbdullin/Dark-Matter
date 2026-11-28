# Numerical forward-only windows in retained cold D8 spectra

Fix each evolved D8 numerical state at x=.05. Refitted transverse pole data are checked against the independent [parent cap scan](INDEPENDENT-EARLY-POLE-ENDPOINT-PAST-K-CAP-APERTURE-2026-11-28.md) at k=.00140, .00145 and .00150. Then evaluate the **weighted** gross forward, inverse, and affine slope of the actual retained charged populations at q=3.5e−6 and q=1.5e−5, including their scalar Bose/occupation factors. These are 12 point diagnostics, not a k- or q-integrated source. At k=.00150 the baseline also reproduces the [independently partitioned parent-energy integral](INDEPENDENT-EARLY-POLE-ENDPOINT-PEER-D8-PARENT-TAIL-2026-11-28.md).

The inverse channel is limited by the daughter occupation support: its largest parent-equivalent energy is the daughter cap minus the positive transfer (omega-sqrt{q^2+m_s^2}). That may be far below the forward parent cap even when many angles are allowed for forward scattering. Independently classify both angular apertures directly from the two caps. Where paired quadratures and finite support give zero inverse, state it as zero **on that retained numerical support**, not a physical absence. At a fully closed point both rates are zero, so no unique conditional balancing occupation can be defined. At a retained forward-only point (B=0,D>0), the root of (n_gamma D-B=0) is zero, not a prediction that the photon population is zero.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-cell-endpoint; d=2026-11-28
cap=results/early-pole-endpoint-past-k-cap-aperture-$d.json
tail=results/early-pole-endpoint-peer-d8-cold-parent-tail-$d.json
python -W error early_pole_endpoint_peer_d8_k_rate_points.py \
  "$p-fine-d8-x05-state-$d.json.gz" "$p-third-d8-x05-state-$d.json.gz" "$cap" /tmp/d8-rate-points.jsonl
python -W error early_pole_endpoint_peer_d8_k_rate_receipt.py \
  /tmp/d8-rate-points.jsonl "$cap" "$p-fine-d8-x05-state-$d.json.gz" \
  "$p-third-d8-x05-state-$d.json.gz" --tail "$tail"
python -W error early_pole_endpoint_d8_forward_inverse_cap_window.py \
  "$p-fine-d8-x05-state-$d.json.gz" "$p-third-d8-x05-state-$d.json.gz"
```

At q=3.5e−6 and k=.00145 the forward coefficient is positive but only (2.5244	imes10^{-4}) (fine D8) and (2.0890	imes10^{-4}) (third D8) of its own k=.00150 value, even though the forward cosine aperture is 0.885 and 0.628. The inverse is exactly zero on both retained supports. At q=1.5e−5 and k=.00145 the forward ratios instead are .0029843 and .0021653 and the inverse is nonzero (ratios .0017270 and .0011167). At k=.00140 and q=3.5e−6 both channels vanish; at q=1.5e−5 a tiny forward-only contribution survives. The raw mesh discrepancy is at most 1.10e−12 relative; the independent fixed-k split baseline agrees within 7.20e−14. Thus allowed angle is not a proxy for the observed numerical weight.
