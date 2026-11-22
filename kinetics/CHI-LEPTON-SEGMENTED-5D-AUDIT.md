# Segmented incoming momentum quadrature for relativistic nonthermal χℓ 2→2 Pauli weak moments (2026-11-22)

## Falsification target and independent method

Earlier [independent 5D negative](CHI-LEPTON-5D-QUADRATURE-NEGATIVE-AUDIT.md) showed an oscillatory but fully physical Fermi occupation whose Gaussian incoming radial heat changes sign under nrad/cutoff changes, and whose inherited coarse Q overestimates a multi-axis refined estimate ~55-fold. A correct quantum H theorem and exact energy do not prevent false heat. Test whether *segmenting p,k by fixed energy intervals* provides a controllable incoming radial quadrature and cleanly separates UV cutoff from loss of resolution. Compare two mappings (direct Gauss in p between energy-guided p-bounds, and Gauss in E with exact dp/dE=E/p Jacobian); keep the original unsegmented method as independent comparison.

For each species mass m, choose **fixed maximum energy-segment width** ΔE≤1 in arbitrary T=1 units. Set energy edges m+jΔE, physical p_j=√(E_j²−m²), then apply nlocal-point Gauss directly to each [p_j,p_{j+1}]. This **exactly reproduces** the free radial 0th measure ∫p²dp=pmax³/3 for any mass with nlocal≥2 because the measure is polynomial in p. Direct Gauss in energy instead uses p(E) and dp/dE=E/p; near massive threshold its square-root cusp leaves finite quadrature error (m=.5, ΔE1, nlocal5: relative radial-volume error +4.52×10⁻⁷ vs ≈2×10⁻¹⁶ using segmented p). This is a measure-normalization gate, not full collision convergence.

The same independent full on-shell χℓ↔χℓ Pauli 5D amplitude, F−R and H as the previous stage, synthetic parameters Mχ=.1,mℓ=0,mV=.1,G=.001, gχgℓ=4×4 numeric convention, nμ8,ncos32,nφ16, Q in arbitrary units. Inputs are prescribed smooth and oscillatory nonthermal Fermi distributions, not cosmologically generated χ1 or χ2.

## Reproduced radial/cutoff/angular results

At cutoff12 (max ΔE=1), oscillatory Qχ for segmented-p nsegments12 and nlocal3,4,5,6 is

| local Gauss nodes / total radial nodes per species | Qχ | H |
|---|---:|---:|
| 3 / 36 | +4.645350308×10⁻¹⁵ | +7.532886595×10⁻¹³ |
| 4 / 48 | +4.499933299×10⁻¹⁵ | +7.541675008×10⁻¹³ |
| 5 / 60 | +4.568755689×10⁻¹⁵ | +7.540452035×10⁻¹³ |
| 6 / 72 | +4.556593552×10⁻¹⁵ | +7.540669156×10⁻¹³ |

Unsegmented uniform-p Gauss nr24 on **identical** nμ8,ncos32,nφ16 gives **−1.398867936×10⁻¹⁴** (incorrect sign vs segmented and higher-resolution), whereas unsegmented nr64 gives +4.55619107×10⁻¹⁵. Segmented 12×6 and unsegmented nr64 therefore agree to ~0.009% on this *fixed angular cutoff/occupation scenario*. Segmented 8×3 is notably poor Q+2.27825×10⁻¹⁴ and 8×4 gives +1.29928×10⁻¹⁵: **segmentation is not a universal cure** without its own nsegments×nlocal convergence. Smooth f is much less oscillatory; 8×3/8×4/8×6 yield Qχ≈(9.06892,9.06836,9.06981)×10⁻¹⁴ with positive entropy. Do not advertise 0.009% as total physical error.

**Corrected cutoff test:** keep ΔE≤1 and nlocal5 when increasing the cutoff, so radial node count *increases* and tail tests are not contaminated by shrinking resolution:

| kinetic cutoff/T | nsegments × nlocal | oscillatory Qχ | oscillatory H | smooth Qχ |
|---:|---:|---:|---:|---:|
| 8 | 8×5 | +5.365357751×10⁻¹⁵ | +7.308681284×10⁻¹³ | +8.759008311×10⁻¹⁴ |
| 12 | 12×5 | +4.568755689×10⁻¹⁵ | +7.540452035×10⁻¹³ | +9.069869394×10⁻¹⁴ |
| 16 | 16×5 | +4.594669672×10⁻¹⁵ | +7.548936845×10⁻¹³ | +9.082503394×10⁻¹⁴ |
| 20 | 20×5 | +4.594452896×10⁻¹⁵ | +7.549190132×10⁻¹³ | +9.082914679×10⁻¹⁴ |

The cutoff16→20 relative Q change is ~0.0047% **at these fixed angular resolutions**, and entropy ~0.0034%; this is *encouraging* on tail truncation but not global continuum accuracy. Fixed 12 segments across changing cutoff, as in an earlier probe, shifted Q20 to 4.301×10⁻¹⁵ simply because each segment widened; this is precisely why UV cutoff and radial resolution must not be conflated.

**Unresolved outgoing angle:** keeping the segmented 12×5 radial rule and cutoff12, oscillatory ncos16/24/32/48/64 produces Qχ=(4.87480,4.60634,4.56876,4.62164,4.65110)×10⁻¹⁵ and H=(6.8560,7.3998,7.5405,7.5429,7.4984)×10⁻¹³. Changes in H and Q remain several percent at coarse angular grids and are NONMONOTONIC across 32→64. Thus neither the continuous nor the projected strong collision is yet physically certified at claimed subpercent precision. High angular and cutoff 3+ independent rates remain the next scientific gate.

## Tests, reproducibility and conclusion

[Radial scheme extended without altering the previous uniform-p default](independent_chi_lepton_5d_quadrature.py), [segmented deterministic gate](chi_lepton_segmented_5d_gate.py), [4 new regressions](../tests/test_chi_lepton_segmented_5d_gate.py), [compact results](results/chi-lepton-segmented-5d-2026-11-22.json), [archive](../archive/stages/2026/11/2026-11-22-chi-lepton-segmented-radial.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_chi_lepton_segmented_5d_gate.py tests/test_chi_lepton_quadrature_gate.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.chi_lepton_segmented_5d_gate --full

Eight relevant local pytest checks passed including the old negative event constraints; no GitHub CI success claimed until verified on this stage. **Scoped positive:** structured energy-guided p segments recover stable exact radial number measure and substantially reduce an identified oscillatory weak-source alias. **Negative:** adaptive outgoing CM-angle physics, true ν/e± masses/degeneracies, Fermi 0/1 limits, production of χ2/real φ and full SM/Zprime/FRW are absent. No physical γ_el, mφ=86 eV ΩDM, BEC, frozen f(q) or CLASS-ready conclusion.
