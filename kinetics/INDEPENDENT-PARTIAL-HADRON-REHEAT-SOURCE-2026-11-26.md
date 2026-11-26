# Measured vacuum hadronic sub-bands in the high-reheat count benchmark

The earlier all-convert scalar count threshold used only charged
lepton injection. Here add disjoint vacuum spectral contributions
to the same 10 GeV vector portal coupled to the electromagnetic
current:

- The [combined bare KLOE π⁺π⁻ data](https://arxiv.org/abs/1711.03085)
  cover s=0.10–0.95 GeV². Convert each bin's cross section to R(s)
  with the corresponding point-muon denominator.
- [BESIII Table II](https://arxiv.org/html/2112.11728) gives inclusive
  R at 14 points. Integrate only the sub-bands 2.2324–3.0800 and
  3.4000–3.6710 GeV in √s. The gap between them is deliberately not
  assigned the continuum interpolation across the J/ψ pole.
- As separate sensitivities, add an isolated narrow J/ψ using the
  [PDG mass, width and branching fractions](https://pdg.lbl.gov/2025/listings/rpp2025-list-J-psi-1S.pdf),
  and a local constant R proxy of 3.497 anchored to the single
  [CLEO point](https://doi.org/10.1103/PhysRevD.76.072008) at
  √s=9.996 GeV. That last constant is used only on 9.986–10.006 GeV
  as a sensitivity trial; it is not a measured continuous band.

Fold each R with the same finite vector Breit–Wigner and two dark
pair thresholds, `2 mχ` and `2 Mχ`, with classical weights for the
vacuum hadronic parts and the existing FD lepton source. Use the
same prescribed visible-radiation expansion and entropy. An
independent reproduction of the earlier KLOE instantaneous kernel
checks the pair-count and deposited-energy normalization.

| Included vacuum spectral prescription | T_R for optimistic 99% of .120 |
| --- | ---: |
| Charged leptons only | 421.43 MeV |
| Add measured KLOE ππ | 373.97 MeV |
| Add measured BESIII continuum sub-bands | 348.62 MeV |
| Add separate narrow J/ψ sensitivity | 342.64 MeV |
| Also add illustrative CLEO point window | 342.57 MeV |

All rows make the optimistic count allocation that **every** dark
particle becomes an 86 eV scalar. They do not compute how to do it.
At the measured KLOE+BESIII threshold, a coherent ± quoted-input-error
shift gives 347.54–349.72 MeV. That is a sensitivity test, not a
covariance fit or confidence interval; the unmodeled spectral
regions and in-medium changes can be more important. The charged
lepton *local* equilibrium FD response at the lepton-only threshold
has a fast mode around 4.1 H, another reason not to mistake the
prescribed vacuum calculation for the actual thermal history.
No complete hadronic spectrum, QCD plasma spectral function,
inverse or backreacting cosmological solution is inferred.

## Input provenance and reproduction

Download KLOE's [combined cross-section data](https://www.lnf.infn.it/kloe/ppg/ppg_2017/KLOEppg2017/CrossSection/Combination/KLOExSec.data)
from the collaboration's cross-section combination directory and
[arXiv:2112.11728v2 PDF](https://arxiv.org/pdf/2112.11728v2).
For the latter, run `pdftotext -layout BESIII.pdf BESIII.txt`.
The script parses Table II. Both checked input hashes are stored
with the numerical output; the measurements are not redistributed.

```sh
curl -fLO https://www.lnf.infn.it/kloe/ppg/ppg_2017/KLOEppg2017/CrossSection/Combination/KLOExSec.data
curl -fL https://arxiv.org/pdf/2112.11728v2 -o BESIII.pdf
pdftotext -layout BESIII.pdf BESIII.txt
OPENBLAS_NUM_THREADS=1 python -W error measured_partial_hadron_reheat_fold_reference.py EOS.late.dat KLOExSec.data BESIII.txt
```
