# Independent full relativistic χ–lepton 5D Pauli quadrature error gate (2026-11-22)

## Question, falsification criterion and model scope

Earlier finite-mass t-channel χℓ→χℓ work showed positive microscopic Pauli entropy, exact pairwise Qχ=−Qℓ and an entropy-compatible energy-bin source. Its 10-bin nonlinear time pilot is **not physical rate evidence** if the 5D continuous-energy integration itself is under-resolved. Test two allowed but NON-evolved isotropic Fermi spectra at fixed masses/couplings by varying **one numerical quadrature axis at a time**, not by selecting favorable grid points. The main falsification criterion is whether successive incoming radial, incoming polar, outgoing CM polar, outgoing azimuth and energy-tail truncations stabilize Qχ and entropy. Energy and H-theorem checks are independent controls, NOT rate convergence criteria.

On-shell χ(M=.1)+ℓ(m=0)↔χ+ℓ with finite vector mV=.1 and G=.001, inherited gχgℓ=4×4 *numerical convention* (not physical ν flavor accounting); independent Gauss p,k,cosμ and CM outgoing cosθ with periodic trapezoid φ. Eχ+Eℓ=Eχ′+Eℓ′, spin trace

    |M|² = 2G²[(s−M²−m²)²+(u−M²−m²)²+2(M²+m²)t] [mV²/(mV²−t)]².

J=F−R with full Fermi Pauli blocking; H=∫W J log(F/R)≥0. This standalone implementation independently reproduces the prior FD χ temperature .9, ℓ temperature 1: Qχ=6.935904667187969×10⁻¹⁴, entropy=7.706560741319965×10⁻¹⁵ at nrad=nμ=4, ncos=nφ=12 and kinetic cutoff12. Separate common-temperature FD with μχ=.2, μℓ=−.1 gives |Q|<3×10⁻²⁸ and H≈2×10⁻⁴² at the same grid.

## Independent experiment

Energy units arbitrary with reference T=1. Spectra: smooth fχ(E)=expit[−E/.9+.15exp(−E/.9)], fℓ(ω)=expit[−ω+.05exp(−ω)]; oscillatory fχ(E)=expit[−E+.3sin(4E)], fℓ(ω)=expit[−ω+.3cos(4ω)]. Both strictly obey 0<f<1. One-axis scans hold the other axes fixed at nrad8,nμ8,ncos24,nφ16, cutoff12. The following Qχ values are in units 10⁻¹⁴ from a synthetic integral; sign changes are **numerical instability**, not model disproof.

| Axis / order | Smooth Qχ (×10⁻¹⁴) | Oscillatory Qχ (×10⁻¹⁴) |
|---|---|---|
| incoming nrad=4,6,8,10,12,16 | 8.2080, 9.0759, 9.0116, 9.0120, 9.0122, 9.0134 | +30.9727, −6.4605, +2.0409, −41.1258, +14.2008, +0.0091 |
| incoming nμ=4,6,8,10,12,16 | 9.0823, 9.0290, 9.0116, 9.0042, 9.0008, 8.9981 | 2.0198, 2.0344, 2.0409, 2.0419, 2.0432, 2.0437 |
| outgoing ncos=8,12,16,24,32,48 | 7.6719, 8.4024, 8.7519, 9.0116, 9.0674, 9.0534 | 2.9528, 3.0730, 2.8634, 2.0409, 1.4583, 1.0691 |
| outgoing nφ=8,12,16,24,32 | 9.0117, 9.0116, 9.0116, 9.0116, 9.0116 | 2.2162, 2.0681, 2.0409, 2.0287, 2.0275 |
| kinetic cutoff/T=8,12,16,20 | 8.7081, 9.0116, 9.0330, 9.0852 | **−35.1256, +2.0409, −25.1233, −27.6789** |

Radial and cutoff oscillatory heat is **not converged even in sign**, while every tested microevent has nonnegative entropy and Qχ+Qℓ exactly zero. In the smooth case outgoing ncos8→32 changes Q by ~18.2%, despite apparently stable incoming nr8→10. A green conservation/entropy test is insufficient.

As a stronger *separate multi-axis* diagnostic (not a continuum certificate), the same oscillatory state gives Qχ:
- coarse nrad4,nμ4,ncos12,nφ12: +2.55856732856×10⁻¹³;
- nrad40,nμ8,ncos32,nφ16: +4.535657994×10⁻¹⁵;
- nrad48,nμ16,ncos64,nφ32: +4.619773235×10⁻¹⁵;
- nrad64,nμ16,ncos64,nφ32: **+4.614781070×10⁻¹⁵**.

The old coarse Q exceeds the latest multi-axis value by **55.44×**. This refined value is CONDITIONAL: at nrad40,nμ12,nφ24, ncos24/32/48/64/96 yields Q≈4.5780/4.5204/4.5683/4.6004/4.6115 ×10⁻¹⁵, with H ~7.40–7.54×10⁻¹³. UV, even narrower physical mediators and IR features require further checks. No real physical production history, cosmological Γ or viability follows.

## Reproduction, negative decision and next gate

Independent [5D kernel](independent_chi_lepton_5d_quadrature.py), [single-axis scan](chi_lepton_5d_quadrature_gate.py), [four regression tests](../tests/test_chi_lepton_quadrature_gate.py), [numeric evidence](results/chi-lepton-5d-quadrature-gate-2026-11-22.json), [dated stage](../archive/stages/2026/11/2026-11-22-chi-lepton-5d-quadrature-negative.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_chi_lepton_quadrature_gate.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.chi_lepton_5d_quadrature_gate --full

All four tests passed locally; this is NOT yet the separately confirmed GitHub Scientific CI on the new commit. The independent kernel agrees with the earlier thermal reference; eventwise |ΔEχ+ΔEℓ| <5×10⁻¹⁴. No RNG/seed/external file is required. The committed script regenerates the complete numerical scan.

**Scientific decision:** previous nrad4 synthetic χℓ Q, entropy and time results MUST be labeled *qualitative/instrument validation only*, not a numerical physical relaxation parameter. Original data remain in history. P0 next: physically adaptive/segmented incoming momentum integration, 3+ resolutions and cutoff controls, resolved outgoing forward angles, independent 5D comparisons, then paired strong χ1/χ2 and finite SM/φ Bose+FRW. 86-eV relic, Γel(T), ΩDM, BEC and CLASS-ready ALL remain unproven; status in_progress.
