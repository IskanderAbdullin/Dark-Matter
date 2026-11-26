# Independent EXACT FD ultrarelativistic direct Yukawa cold-daughter fraction and finite-mass controls (2026-11-24)

## Physical question and strict boundary

For the existing on-shell massive Dirac χ1(.1GeV)→Dirac χ2(.02GeV)+real φ(86eV) illustrative Yukawa benchmark, a fully coupled 1↔2 Bose-stimulated spectrum can only be interpreted after actual nonthermal evolution. Yet a useful **independent fixed-temperature vacuum direct-decay reference** asks:

> What proportion of daughter scalars born in **ordinary, unstimulated direct decays** of a thermal Fermi-Dirac parent have laboratory energy kφ/T below a prescribed cold threshold α?

This is a *microscopic one-way weak diagnostic*, NOT the fraction of eventual cosmological DM, a Bose enhancement calculation, a zero-mode, SM/FRW, number-conserving finite-time history, or a transfer function. We deliberately switch off Pauli final-state blocking (f2=0), Bose final-state occupation (bφ=0), and reverse fusion to isolate the kinematic parent-origin cold tail. These processes must be restored in the real physical relic model.

## Analytic ultrarelativistic distribution: actual FD parent and time dilation

Use a fermionic parent thermal FD distribution f1(p)=1/[exp(p/T)+1] when T≫M1, and the correct vacuum time-dilated event weight

    dNdecay ∝ p² dp (M1/E1) Γrest f1(E1)
             ∝ [q/(exp(q)+1)] dq,  q=p/T.

For an isotropic two-body rest-frame decay, the observed massless-in-this-limit daughter scalar energy is uniform over

    0 < kφ < A p,
    A=(Eφ*+k*)/M1≈0.96

(the small positive 86-eV mass and the exact boosted kmin are retained in the separate finite-T integral below). Thus at fixed parent q, the fraction emitted with scalar x=kφ/T<α is min(1,α/(Aq)).

The FD-weighted numerator and denominator can be **independently integrated analytically**:

    D=∫_0^∞ q/(e^q+1)dq = π²/12,
    P(x)= [12/(A π²)] ln[1+exp(−x/A)], x>=0,
    F_cold(<α)=1+[12/π²] Li₂(−exp[−α/A]).

The normalized direct-daughter energy density P(x) is NOT an equilibrium Bose occupation, and the answer is independent of y because the same *one-way Γrest* cancels from the fraction. Its mean scalar energy is

    <kφ/T>_direct = 9 A ζ(3)/π² = **1.05229867594**.

For comparison, a *separate* exactly thermal relativistic BE scalar occupation has <k/T>=π⁴/[30 ζ(3)]≈2.70118; neither expression is the actual late cosmological f(q) of the Yukawa model.

| Cold threshold α | Exact ultrarel FD direct-decay fraction F(k/T<α) |
|---:|---:|
| 0.001 | **0.000877551395** (0.087755%) |
| 0.01 | **0.008745886701** (0.874589%) |
| 0.1 | **0.084547145515** (8.454715%) |

These are strictly **conditional instantaneous direct emission fractions**. In particular, a very large stimulated low-k condensate is NOT implied by the small direct-cold fraction; inverse/scattering/thermal history can qualitatively change the spectrum.

## Independent MASSIVE scalar-parent-kinematics FD numerical reference

Retain the exact massive boosted daughter shell from the already validated physical Yukawa on-shell event measure

    E1=sqrt(p²+M1²),
    kmax=(E1 Eφ*+p k*)/M1,
    kmin=[p² mφ² + M1² Eφ*²]/[M1 (E1 Eφ* +p k*)].

At fixed p, probability of kφ<αT is the *exact* clipped uniform rest-angle fraction

    f_cold(p,K)=max(0,min(1, (K−kmin)/(kmax−kmin))), K=α T.

Normalize the parent FD decay weight p²(M1/E1)dp/[exp(E1/T)+1]. We integrate the full parent tail **starting at the exact inverse-rapidity p_low(K)**; not a sparse parent p grid that accidentally gives zero for rare but physically allowed states. A physically positive high-parent cutoff is chosen only after the exact threshold (parent q_min+70), while denominator is integrated over q≤75; the ignored FD tail is exponentially tiny in the stated comparisons.

The independent massive, numerical FD angle/parent integrals give:

| Parent T (GeV) | Cold F(k/T<.001) | F(k/T<.01) | F(k/T<.1) | Mean k/T of direct daughter |
|---:|---:|---:|---:|---:|
| 0.1 | **4.65676e−110** | **2.76587e−14** | 0.00386116 | 1.339497 |
| 1 | 2.79493e−5 | **0.00521753** | 0.07842246 | **1.061156** |
| 100 | 0.000876154 | 0.008744159 | 0.084545259 | 1.052301 |
| 1000 | **0.000877530** | **0.008745862** | **0.084547120** | **1.052299** |
| T/M1→∞ (analytic) | 0.000877551 | 0.008745887 | 0.084547146 | 1.052299 |

The T=1000GeV numerical fractions agree with exact ultrarel formulas within ≈3e−5 relative across α=.001/.01/.1; energy mean agrees within ≈2e−8 relative.

**New negative control:** at T=.1GeV and α=.001 the exact physical parent threshold is **p_min≈23.99990 GeV** or q_min≈239.999. A simplistic parent quadrature ending at q=75 would report an **exact zero** cold production rate, but the real FD tail gives small yet strictly positive conditional fraction **4.65676e−110**. This is a mathematically different issue from an energy bin that is genuinely forbidden under an explicitly hard parent cutoff P; accidental numerical truncation must not be interpreted as exact kinematic decoupling or a condensate.

At T=.1GeV, α=.01 the fraction is 2.77e−14, while at T=1GeV it is .00522 and in the high-T limit ≈.00875. Such strong conditional temperature dependence is a reminder that the original fixed-GeV-temperature thermal near-86eV source bound **cannot** be applied unchanged to much hotter early-Universe epochs. There is still no calculation of the actual time-integrated cosmological cold fraction; that requires the full strongly coupled F/R and FRW dynamics.

## Reproduction and independent checks

- [Physical massive FD parent / on-shell exact angle numerical source + closed-form Li₂](phi_bose_direct_cold_fraction.py)
- [Five independent source/analytic cross-check tests](../tests/test_phi_bose_direct_cold_fraction.py)
- Commands: OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_direct_cold_fraction.py; python -m kinetics.phi_bose_direct_cold_fraction.
- Tests cover direct integration of the normalized probability density and mean, analytic dilog CDF at five α thresholds, exact FD massive-to-ultrarel three-α convergence, the ultracold nonzero q≈240 thermal tail, and the strict difference between direct fraction and full Bose-stimulated dark relic.

New public Scientific CI must be checked for the source/test commits; do not equate software checks with physics. The primary source of the effect (Bose-enhanced cold production) is discussed by [Yin, JHEP05(2023)180](https://doi.org/10.1007/JHEP05(2023)180) and related 2→2 stimulated scattering in [Sakurai–Yin, JHEP03(2025)202](https://doi.org/10.1007/JHEP03(2025)202), but these works do NOT independently validate our particular physical y=.001 finite-T fraction; the scalar FD kernel above is independently derived and tested.

**NEXT_ACTION:** use this normalized *direct only* FD p-angle cold-fraction curve as a source benchmark for the actual physically coupled three-field quantum implicit collision solver at thermal source conditions; add Fermi blocking, inverse rates, Bose-stimulated occupation, SM finite bath and FRW, track comoving q and distinguish cold/hot number and energy fractions with parent UV/qmin/time convergence. Do not use the analytic fraction as late-time Ω, CLASS-transfer or evidence of BEC.
