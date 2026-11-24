# Conditional 86-eV Yukawa thermal-stability quartic versus halo self-scattering (2026-11-24)

## Physical test

The preceding independent thermal one-loop two-Dirac off-diagonal Yukawa analysis in PHI-BOSE-YUKAWA-FERMION-THERMAL-MASS-AND-Z2-GATE-2026-11-24.md found a negative STATIC real-scalar curvature around T=3.245 to 35.482 MeV, conditional on both Dirac species being in a common FD bath, with m1=0.1, m2=0.02 GeV, y=0.001 and vacuum real scalar mass mphi=86 eV. With no other positive plasma contribution, a scalar quartic potential lambda*phi^4/4! and its leading thermal self-energy lambda*T^2/24 would need lambda at least 1.08515714e-6 to keep local phi=0 curvature nonnegative throughout the examined interval. This is NOT an established cosmic phase transition or a real-time scalar pole mass.

But the SAME quartic yields tree-level identical REAL scalar contact scattering. Accounting for final 1/2!, the differential Born cross section is lambda^2/(128*pi^2*s) per solid angle. In the nonrelativistic limit s->4*mphi^2, total sigma=lambda^2/(128*pi*mphi^2). With 1 GeV^-2=3.89379e-28 cm^2 and 1 GeV=1.78266192e-24 g, this implies

- lambda_min = 1.08515714e-6.
- sigma_NR = 0.3959392 GeV^-2.
- sigma_NR/mphi = **1005.62027 cm^2/g**.
- Isotropic viscosity-weighted sigma_V/mphi = (2/3)*sigma/mphi = **670.41351 cm^2/g**.

The simple Born contact elastic amplitude is derived independently of the thermal determinant. The derivation assumes a dominant dilute nonrelativistic real-scalar dark halo, not the unestablished actual f(q) or condensate.

## Cluster comparison as a scoped consistency obstruction, NOT a new likelihood

Published cluster-merger SIDM constraints are model and hydrodynamics dependent: [hydrodynamic Bullet Cluster analyses](https://academic.oup.com/mnras/article/530/3/3160/7649378) discuss bounds around 2–5 cm^2/g, and [2026 JWST Bullet Cluster analysis](https://doi.org/10.1051/0004-6361/202659214) discusses stringent and relaxed interpretations. We therefore compare to **illustrative** reference maximum sigma/m, not a universal statistical 95% exclusion.

| Reference maximum (cm^2/g) | Corresponding lambda maximum at 86eV | Required stabilizing lambda divided by this maximum |
|---:|---:|---:|
| 0.2 | 1.53035e-8 | 70.9 |
| 1 | 3.42197e-8 | 31.7 |
| 2 | 4.83939e-8 | 22.4 |
| 5 | 7.65175e-8 | 14.2 |

Even under a generous reference of 5 cm^2/g the simple stabilizing lambda predicts roughly 201 times more Born sigma/m than the reference. Conversely at lambda allowed by that reference, the same CONDITIONAL equilibrium-fermion static curvature remains negative at T=20 MeV. Holding lambda_min fixed purely algebraically, sigma/m=1 would require mphi about 862 eV, but this is NOT a physical mass scan since the thermal-stability condition also changes with mass.

This is a genuinely useful **model-closure tension**, not a blanket exclusion of 86eV Bose dark matter. Out-of-equilibrium fermions, other positive portal thermal mass, full effective potential and nonperturbative phase history, partial phi abundance, quantum occupation and opacity can change its physical relevance. No physical freezeout, BEC or CLASS/halo likelihood has been established. Cross sections may require appropriate viscosity/identical-boson transport factors, but factors of order unity do not resolve O(200) difference.

## Independent reproduction

Source: kinetics/phi_bose_quartic_halo_tension.py, tests: tests/test_phi_bose_quartic_halo_tension.py, numerical JSON: kinetics/results/phi-bose-quartic-halo-tension-2026-11-24.json, permanent archive: archive/stages/2026/11/2026-11-24-yukawa-scalar-quartic-cluster-tension.md.

Commands: python -m kinetics.phi_bose_quartic_halo_tension ; python -m pytest -q tests/test_phi_bose_quartic_halo_tension.py.

NEXT_ACTION: derive self-consistent thermal and nonequilibrium fermion occupancies, scalar real-time pole and quartic/SM portal; compute true 1<->2 plus pair 2<->2 quantum collision rates and freezeout f(q), then compare predicted halo scattering including astrophysical model uncertainties. CLASS-ready remains false.
