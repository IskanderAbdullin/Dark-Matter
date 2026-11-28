# Use a measured positive moment to bound a missing pion slice

The [new BABAR publication, Table 1](https://arxiv.org/html/2611.04596v1)
reports the 0.30–0.50 GeV integral for `a_μ^ππ`, rather
than a numerical spectrum usable over the narrow low gap.
Use the bare, final-state-radiation-inclusive convention
and the [standard positive dispersive kernel](https://arxiv.org/pdf/1006.5313),
`H(s) = K_μ(s)/(4π³)`, with cross sections in GeV⁻².
No actual pointwise spectral shape is reconstructed.

For the **nonoverlap only**, .30 GeV≤√s≤√.10 GeV, let
`W(s)` be the expansion-integrated positive 50 MeV warm
vacuum energy weight, including the previously bounded
massive two-pion Bose factors. Positivity gives

`∫_gap W σ ds ≤ sup_gap(W/H) × a_μ^ππ(0.30–0.50 GeV)`.

Assigning the whole experimental moment to the gap
overestimates the input there. Do **not** add a BABAR source
above √.10: the [existing KLOE fold](INDEPENDENT-WARM-PION-BAND-AND-SOURCE-HEADROOM-2026-11-26.md)
already supplies the overlapping region. The referenced
baseline has Bose-weighted KLOE but no subthreshold
form-factor cap or additional pion heat, so the gap is
not counted twice. Use the experimentally reported
moment plus the stated multiple of the two quoted
uncertainties *added linearly* as a conditional cap.
This is not a statistical confidence limit.

The finer positive cells give `sup(W/H) ≤ .01821833 GeV⁴`;
the coarser version gives .01880438. The maximum controlling
cell is next to the upper gap edge. A separate thermal
normalization agrees with the existing hadronic kernel,
and the conservative muon-kernel lower cells differ by
less than .67% from an independent Gauss calculation.

| Assumed whole moment | Extra gap energy upper, GeV⁴ | Share of remaining positive-source budget |
| --- | ---: | ---: |
| New BABAR, central plus one linear error sum | 1.1623×10⁻¹⁰ | 10.50% |
| New BABAR, central plus two linear error sums | 1.2188×10⁻¹⁰ | 11.01% |
| BABAR combined, central plus two linear error sums | 1.1201×10⁻¹⁰ | 10.12% |

With the wider new-BABAR allowance, the previously
conservative early-deposit comparison still has
9.852×10⁻¹⁰ GeV⁴ of headroom. Convex interpolation of
the quadratic energy bound gives an upper scalar
Ωh²≤.03174 when only this gap is added. The result is
specific to the effective two-pion vacuum input and its
population bound. The region below .30 GeV, other
hadronic and radiative channels, and any in-medium
departure from the vacuum kernel remain separate.
No BABAR/KLOE spectral combination or statistical
consistency assumption is made.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_babar_dispersive_gap_reference.py EOS.dat \
  "$D"/warm-pion-bose-and-threshold-2026-11-26.json
~~~
