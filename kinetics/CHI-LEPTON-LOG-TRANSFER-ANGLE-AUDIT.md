# Relativistic finite-vector χℓ 2→2: analytically compensated outgoing CM polar-angle integration (2026-11-22)

## Question, physics and prespecified independent checks

Energy-guided incoming p Gauss now avoids gross oscillatory fχ,fℓ aliases, but the full Pauli 5D heat/entropy remain nonmonotone under outgoing CM polar Gauss refinement. For a finite *t-channel vector* the squared propagator has strong forward enhancement, so uniform cosθ quadrature can miss essential small |t|. A partial analytic cancellation that makes only the zeroth angular factor exact need NOT resolve the **energy-transfer-weighted** moment. Test an exact angular *logarithmic transfer variable* against two independent closed analytic vector-propagator moments and then against fully on-shell χℓ↔χℓ Pauli Q/S at 3+ outgoing angular resolutions. Retain all other numerical uncertainties and compare two physical mediator masses.

For fixed incoming invariant s, center-of-mass outgoing momentum p*, define c=cosθ*, z=1−c, A=2p*²/mV² >0 and t=−2p*²z. The propagator squared factor is [mV²/(mV²−t)]²=(1+A z)^−2. Its two **independent analytic** CM angular controls are

    I0(A) = ∫_{−1}^1 dc/(1+A(1−c))² = 2/(1+2A),
    I1(A) = ∫_{−1}^1 (1−c)dc/(1+A(1−c))²
          = [log(1+2A)+1/(1+2A)−1]/A².

I1 is the leading momentum/energy-transfer angular **proxy**, not the complete physical 5D χ drag. For very small A it has the stable series 2−16A/3+12A²−128A³/5+…; for all tested positive A the exact formulas are independent of the quadrature under audit.

Three integration representations:
1. Uniform Gauss in c with explicit (1+A z)^−2.
2. Uniform u=1/(1+A z): cancels propagator² exactly for I0, but leaves a difficult logarithmic I1 integrand.
3. **Log-transfer** y=log(1+A z), y∈[0,log(1+2A)], c=1−expm1(y)/A. Then dc/(1+A z)² = exp(−y)dy/A. With Gauss x∈[−1,1], y=(x+1)log1p(2A)/2, use nonnegative effective weights w_x log1p(2A)exp(−y)/(2A) and the amplitude's numerator **without** the already absorbed propagator. No energy/Pauli/H clipping, no rate post-normalization. The transformation is exact at continuous integration level, including for A→0; it does not approximate a physical BGK Γelastic.

## Independent analytic measurement, including a falsification of I0-only schemes

At ncos=8, normalized absolute errors for the two exact angular controls:

| A | Uniform c I0 rel.error | Uniform c I1 rel.error | Propagator-cancel u I1 rel.error | Log-transfer I0 rel.error | Log-transfer I1 rel.error |
|---:|---:|---:|---:|---:|---:|
| 1 | −1.69×10⁻⁸ | +2.34×10⁻⁸ | −2.66×10⁻⁹ | ~2×10⁻¹⁶ | ~0 |
| 10 | −1.196% | +0.396% | −0.148% | ~0 | ~0 |
| 100 | **−51.56%** | +1.334% | −10.57% | −2.75×10⁻¹² | +6.35×10⁻¹³ |
| 1000 | **−93.11%** | −19.77% | **−33.87%** | −4.37×10⁻¹⁰ | +6.62×10⁻¹¹ |
| 10000 | **−99.28%** | −39.11% | **−50.26%** | −1.42×10⁻⁸ | +1.60×10⁻⁹ |

Thus making only I0 exact is insufficient: the u scheme still misses I1 by 50.26% at A=10⁴. Eight-node log-transfer gives ~1e−8 relative or better for both on-shell analytic controls over A=.01–10000. This is a scoped mathematical validation of angular integration.

## Actual 5D nonthermal Pauli weak moments: positive and negative

The **same** independent incoming energy-guided p/k segmented rule and relativistic outgoing CM boost are used; only the angular rule changes, all full Pauli F/R factors and eventwise Qχ+Qℓ, Sdot≥0 remain. Synthetic fixed Mχ=.1,mℓ=0,G=.001,Tref=1, energy cutoff12, historical numeric gχgℓ=4×4. Two mediator cases:
- mV=.01, smooth nonthermal f, incoming segments8×4, nμ6, nφ12.
- mV=.1, oscillatory nonthermal f, incoming segments12×5, nμ8, nφ16.

**Narrow mediator mV=.01**, Qχ in arbitrary units:

| CM polar ncos | Uniform c Qχ | Log-transfer Qχ |
|---:|---:|---:|
| 8 | 8.845413038×10⁻¹⁸ | 1.799941528×10⁻¹⁷ |
| 16 | 1.140027874×10⁻¹⁷ | 1.799939450×10⁻¹⁷ |
| 32 | 1.387397025×10⁻¹⁷ | 1.799939448×10⁻¹⁷ |
| 64 | 1.600559905×10⁻¹⁷ | 1.799939448×10⁻¹⁷ |

At ncos32 the uniform source **underestimates log-transfer ncos64 by 22.92%**, while log-transfer ncos8 differs from its ncos64 by only **0.000116%**. Entropy has the same controlled behavior (log-transfer ncos8 2.2362502×10⁻¹⁸ vs ncos64 2.2362441×10⁻¹⁸). This is an unusually strong, independently testable numerical improvement without changing physical Pauli distributions.

**Larger mediator mV=.1**, oscillatory synthetic Qχ for log-transfer ncos8/16/32/64 = (4.6749992,4.6625450,4.6587064,4.6586688)×10⁻¹⁵, entropy = (7.4321991,7.4323698,7.4323836,7.4323860)×10⁻¹³. Log-transfer ncos32→64 heat shifts about **0.00081%** for these *fixed incoming/UV parameters*, while the earlier uniform-c ncos16/24/32/48/64 series was nonmonotone (Q 4.8748/4.6063/4.5688/4.6216/4.6511×10⁻¹⁵). A different numerator or f(E) could require a different ncos: these are demonstrable **specific cases**, not universal convergence.

## Scope and scientific decision

The outgoing angle transformation reduces the dominant spurious forward peak error at **two mediator-to-energy ratios** while preserving exact four-momentum/Pauli gain-loss and microscopic H. It does not itself prove the *complete* χ1/χ2→SM collision: time/momentum-space Pauli strong Cχ and Cℓ, physical spin/flavors, actual mediator population, χ1→χ2+real φ Bose, finite e±γν heat/entropy, thermal masses, H/FRW and realistic coupling constraints are not validated. In particular, the χ+ℓ weak source uses prescribed synthetic Fermi occupation functions rather than the cosmologically produced daughter fχ2. Full cutoff, incoming nμ, outgoing φ and angular mapping need **joint** independent tests across both shapes and mℓ=.5 before physical Γel/H claims.

Reproduce [source with all three angular schemes](independent_chi_lepton_5d_quadrature.py), [analytic gate](chi_lepton_outgoing_log_angle.py), [10 tests](../tests/test_chi_lepton_outgoing_log_angle.py), [compact measured data](results/chi-lepton-outgoing-log-angle-2026-11-22.json), [dated archival stage](../archive/stages/2026/11/2026-11-22-chi-lepton-log-transfer-angle.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_chi_lepton_outgoing_log_angle.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.chi_lepton_outgoing_log_angle --full

All 10 tests passed locally; actual GH Scientific verification status must be independently read after commit. Negative results from ordinary Gaussian and from the I0-only compensation remain archived; no spurious physical Γel/86eV ΩDM/BEC/CLASS-ready assertion. The next meaningful stage is a **joint** incoming radial+μ and mapped outgoing (c,φ) cutoff/finite-lepton mass error budget, followed by coupled two-χ and real-φ Bose quantum kinetics.
