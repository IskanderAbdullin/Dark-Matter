# Ход работы: теоретическое продолжение III

2026-11-21 — no GitHub Actions, no Computer Use, no changes to main; отдельная ветка для защиты общего info.

1. **Source provenance:** непосредственно распакован исходный архив \`rv2_JCAP_041P_0721(1).zip\`, прочитан TeX (≈70k characters), в том числе BECdens, Rc, RbsSIDM, distribBE, dispV0, BE equations, source background and conclusions.
2. **Algebraic radius audit:** subbing \(\sigma=8\pi a^2\) into \(R_c=\pi\sqrt{\hbar^2a/(Gm^3)}\) gives \(\pi^3/8\) in fourth power; later RbsSIDM gives \(1/(8\pi)\). Ratio \(\pi\), reproduced physically at m50/86; labelled *definition inconsistency to verify*, not blanket error of DM hypothesis.
3. **Microscopic scales:** with source F583-4 R0.6kpc,V0 44.1km/s, hypothetical m86eV, \(\lambda_{\rm th}=39.10\) µm, separation=28.30 µm, \(n\lambda^3=2.637\); tension with \(\nu_c=2.612\) is source rounding (m=86.225 needed at source ν2.61).
4. **Bose occupations:** normalized exact critical 3D momentum distribution. Independently integrated fractions \(f_{\rm mode}>1,10,100,1000\), highlighting limits of classical field kinetic transfer.
5. **New momentum moments:** source V0²=kT/m is not actual σ1D² in degenerate case; exact Li5/2/Li3/2, critical σ1D/V0=0.716598, excess kurtosis≈1.906876 from Li7/2. Implications for source interpretation, not a claimed violation of potential equation.
6. **Local isothermal HF scale:** formal \(c_T^2=V_0^2 g_{3/2}/g_{1/2}+2gn/m\), critical \(\lambda_{J,T}=2\pi R\sqrt b\) and ratio to full TF first-zero radius √8; **local diagnostic**, not dynamical stability.
7. **Cosmological split:** derived conditional matter-era growth exponent p=[√(25−24 f_w)−1]/4 for a smooth warm fraction on given k; emphasized \(f_w\) primordial != late noncondensed fraction, not a numerical Lyman-alpha bound.
8. **Literature:** searched/verified editor records and abstracts across cosmology, Bose thermal/hydrodynamics, wave kinetic, King stability, 2026 SPARC covariance and stimulated eV dark-matter genesis; 44-item curated crosswalk with provenance cautions.
9. **Theoretical scenarios:** standard thermal, Bose-stimulated cold generation, early condensate / late gas, stationary Vlasov DF. Each scenario has independent measurable go/no-go conditions and failure modes.
10. **Stability/kinetic:** derived sign \(df/dE<0\) for a lowered Bose DF and distinguished it from a complete stability theorem. Wrote schematic quantum Boltzmann–Poisson requirements, no pretence solver executed.
11. **Reproducibility:** purely local mathematical consistency code checked numerical pi, moments, HF/TF scales and conditional growth. No data fit/CLASS run in this continuation.
12. **Next:** confirm radius convention against original coauthor prior papers; construct full momentum-dependent production and kinetic mapping; mathematically rigorous stability of lowered Bose; revisit JCAP phase-space probability and common-mass inference with systematics.

### Negative-result discipline
- No CLASS/CAMB transfer function, observational exclusion, physical particle Lagrangian, halo microphysical formation or stability theorem established.
- 40 thematic references are deliberately *annotated screening*, not 40 complete PDF full-text reviews.
- No claim that a Bose-gas momentum-moment identity or classic King lowering is mathematically new.

13. New close prior art: 2023/2024/2025 Proukakis–Rigopoulos–Soto GP/Boltzmann/Poisson + stochastic/linear perturbations (full framework prior art); 2022 Chavanis Casimir constrained collisionless relaxation. 44 curated studies, not full text reading of all.
14. Added stand-alone Python numerical reproduction of exact Bose moments, occupancy, TF radius pi issue and conditional mixed growth, with explicit tolerances and an archived numerical JSON summary. No GitHub Actions.
