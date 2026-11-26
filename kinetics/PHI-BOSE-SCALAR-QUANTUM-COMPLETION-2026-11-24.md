# Scalar quantum completion: fermionic Coleman-Weinberg matching, self-interaction and cold field diagnostics

## Missing independent scalar parameter
Benchmark real scalar with off-diagonal yφ(barχ1χ2+h.c.), m1=M=.1GeV, m2=m=.02GeV, y=.001, physical desired mφ=86e-9GeV. Combined Z2 (φ→−φ, χ2→−χ2) allows φ² and φ⁴, forbids odd powers if exact; no shift symmetry protects mφ. Renormalizability requires scalar mass and quartic counterterms. A specified mφ and y alone do not determine the physical self interaction. Do not treat renormalized λ as zero without a matching condition. No claim that this is the unique UV completion.

For constant field φ, fermion eigenmasses z±=(M+m±sqrt((M−m)²+4y²φ²))/2. For two Dirac species MSbar one-loop fermion potential is −sum z^4[ln(z²/Q²)−3/2]/(16π²). General convention verified in [Vevacious, EPJC 73 2588 (2013)](https://doi.org/10.1140/epjc/s10052-013-2588-2). At φ≈0, Δ=M−m, x=y²φ², z+ = M+x/Δ−x²/Δ³+2x³/Δ⁵+..., z−=m−x/Δ+x²/Δ³−2x³/Δ⁵+... . Set F(z)=z⁴(ln(z²/Q²)−3/2), a=F'(M)−F'(m), b=F''(M)+F''(m), c=F'''(M)−F'''(m). Expansion coefficients in x: a/Δ; −a/Δ³+b/(2Δ²); 2a/Δ⁵−b/Δ⁴+c/(6Δ³). Convert to δm²φ²/2 + δλφ⁴/4! + c6φ⁶/6!. These are potential (zero-momentum) matching contributions, not full momentum-dependent pole mass/coupling. The renormalized m² and λ require independent matching. c6 finite one-loop fermionic contribution only, independent of Q at fixed masses/y; general UV EFT may contain independent higher operators.

| Q GeV | δm² GeV² | δλ | c6 GeV^-2 | |δm²|/mφ² |
|---|---|---|---|---|
| 0.01 | -2.281044e-9 | -1.869692e-12 | 2.122576e-16 | 308416 |
| 0.02 | -1.410186e-9 | -1.448309e-12 | 2.122576e-16 | 190669 |
| 0.05 | -2.589738e-10 | -8.912714e-13 | 2.122576e-16 | 35015.4 |
| 0.1 | 6.118843e-10 | -4.698884e-13 | 2.122576e-16 | 82731.8 |
| 0.2 | 1.482742e-9 | -4.850546e-14 | 2.122576e-16 | 200479 |
| 1 | 3.504812e-9 | 9.299155e-13 | 2.122576e-16 | 473879 |

Independent 70-digit direct differentiation agrees with analytic derivatives to <1e−11; six scales; dδm²/dlnQ=1.256382677e−9GeV², dδλ/dlnQ=6.079271019e−13 match analytic Tr(M_f⁴) log derivative. c6 invariant to 1e−10 relative. At Q=M, δm²=6.118843e−10 vs mφ²=7.396e−15: ratio8.27318e4, curvature correction magnitude equivalent to24.736keV. In this scheme retaining 86eV requires a compensating parameter at roughly 1.2e−5 relative unless UV mechanism protects it; scheme-dependent naturalness diagnostic, not observed pole mass or exclusion. δλ at Q=M negative. Do NOT identify sign or magnitude of full physical λ with this correction. Positive illustrative total λ=|δλ|=4.698884e−13 is only a scale diagnostic.

## Connection to halo and number change
With V⊃λφ⁴/4!, identical scalar vacuum threshold σ=λ²/(128π mφ²), scattering length a_s=λ/(32πmφ). Repulsive dilute NR coherent TF n=1 polytrope R=π sqrt(λ/(32π G mφ⁴)). For positive illustrative λ=4.698884e−13, σ/m=1.88555e−10 cm²/g and TF R=.467672 AU. Independent prior halo benchmark σ/m=1.25 cm²/g requires λ=3.82587e−8, TF R133.447AU (same prior scale). To support 1kpc pure TF object with m86eV requires λ≈91403, outside perturbative validity; this does not exclude noncondensed galaxy profiles or compact-star populations.

For NR incoherent 4→2, ten tree quartic diagrams partition six external lines into two triples. Four internal q²=9mφ², six q²=−3mφ²; sum denominators −1/mφ². Thus signed threshold amplitude λ²/mφ²−c6 if effective local sextic included (derivative higher operators omitted). With 1/(2!4!) phase-space symmetry, event R=√3 n⁴ |M|²/(12288π mφ⁴). Quartic-only exactly matches [Arcadi et al. JHEP 08 (2019) 050](https://arxiv.org/abs/1906.07659), eq55 and appendix C; their formula includes both factorials. This is the unoccupied-final incoherent NR approximation; do not use it as bound for relativistic 2→4, coherent field or Bose-enhanced gas.

If all observed DM number were already present in a cold harmonic zero-mode at T_SM, naive amplitude²=2n/mφ and quartic/quadratic ratio λ n/(6mφ³). For positive illustrative loop-size λ, gS=10.75, Tγ0=2.34865e−4eV and gS0=3.91, this ratio reaches 1 at T_SM≈6.971MeV. At TR20MeV, n≈1.918e−7GeV³, amplitude≈2.112GeV, ratio23.615: quadratic extrapolation inconsistent. At 5MeV ratio.36899; at1MeV .0029519. Formal dilute 4→2 depletion rate/H at 1MeV4.52e−5 for that λ, but coherence/occupation and expansion may change it; at 20MeV formal value5.78e4 is invalid as cold quadratic/dilute approximation. Use exact coherent mixed-potential expansion or full interacting kinetic theory, not this naive extrapolation. No physical origin of such coherent field claimed.

Next: solve adiabatic homogeneous oscillations with actual renormalized λ as parameter, assess self-resonance and scalar number-changing/missing visible channels. The energy recycling no-go was explicitly conditional on a selected tree particle network, not a closed scalar completion. GOAL and CLASS-ready=false unchanged. Reproduce python -m kinetics.phi_bose_scalar_quantum_completion (optional --high-precision requires mpmath).
