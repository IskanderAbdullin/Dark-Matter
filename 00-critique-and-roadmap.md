# Independent critical audit, research questions and roadmap

## Main physical model reconstructed from the supplied TeX
For a noncondensed ideal Bose component in spherical gravitational equilibrium, \`f(p,r)=[(2πℏ)^3(exp[(p²/(2m)+mΦ−μ)/(kT)]−1)]^−1\`, \`n=λ_dB^−3 Li_{3/2}(z)\`, \`z=exp[(μ−mΦ)/(kT)]\`, \`λ_dB=√(2π)ℏ/(mV0)\` and \`V0²=kT/m\`. Source eqs. (distribBE), (nonconddens), (deBroglie), (z), (Bose-gas) give \`ν=n0 λ_dB³∈(0,ζ(3/2)]\`, \`R²=V0²/(4πGmn0)\`. With \`y=m[Φ(r)−Φ(0)]/(kT)\`, \`ξ=r/R\`, \`z=z0 exp(−y)\`, the stable equivalent spherical system is **\`(ξ²y')'=ξ² Li_{3/2}(z0 exp(−y))/ν\`**, \`y(0)=y'(0)=0\`. From Gauss's law, dimensionless circular speed satisfies **\`u²=ξ y'\`**.

With a heuristic, isotropic, isothermal population of compact condensates, the source adds a Bose-star density \`x_BS=α exp(−γ y)\` with \`γ=V0²/v_BS²\`, giving \`(ξ²y')'=ξ²[x_n+α exp(−γy)]\`. Do not impose \`α=γ\` as a physical identity: source uses it as a simplifying fitting ansatz motivated by an order-of-magnitude statement.

## Results verified in a new calculation
1. Exact degeneracy limit \`ν_c=ζ(3/2)=2.612375348685488\`.
2. At \`ν=ν_c\`, \`Li_{3/2}(e^(−y))=ζ(3/2)−2√π√y+O(y)\` and regular Poisson potential gives \`y=ξ²/6+O(ξ³)\`. Consequently, **\`ρ_n(r)/ρ_n(0)=1−[2√π/(√6 ζ(3/2))](r/R)+O((r/R)²)\`** for \`r≥0\`. The one-sided coefficient is **0.5539795458**. Finite central density does *not* imply a differentiable three-dimensional spatial density at exact criticality. For \`ν<ν_c\`, the density instead begins with a quadratic term at the mathematical center; the limiting operations \`r→0\` and \`ν→ν_c\` do not commute.
3. The published table's galaxy F583-4 (\`R=0.6 kpc\`, \`V0=44.1 km/s\`, \`ν=2.61\`) reproduces **86.2246 eV** using the article's analytic mass formula; rounding agrees with its 86 eV entry. This is a formula cross-check **not** a reproduction of an observational fit.

## Critical vulnerabilities
**A. Cosmological free-streaming:** an 86 eV thermal relic is not equivalent to a keV thermal WDM relic. Modern Lyman-alpha work gives \`m_thermal ≳5.7 keV\` (95% CL) for a specified thermal relic and IGM model, while alternate priors yield 4.1 keV: [Rogers et al., PRD 109 043511 (2024)](https://doi.org/10.1103/PhysRevD.109.043511). These values cannot simply be transferred to Bose statistics; nevertheless they make a self-consistent Boltzmann/transfer-function computation the first go/no-go condition. If a much colder nonthermal primordial distribution or an early cold condensed fraction is proposed, its production, abundance, and timing must be derived, not assumed.

**B. Abundance / degrees of freedom:** for *one real scalar*, a relativistically decoupled thermal population with \`μ≈0\`, entropy conservation, and no cold condensate yields approx. \`Ω_b h²=(2/3)(m/93.14 eV)(10.75/g_*S,dec)\`. With \`m=86 eV\`: \`g_*S,dec=25 → Ωh²≈0.265\`; \`g_*S,dec≈55 → Ωh²≈0.120\`. This simplified estimate has explicit assumptions and is not the original article's interacting/mixed-phase calculation.

**C. Infinite-mass halo:** dilute asymptotic Bose gas tends to an isothermal envelope, \`ρ∝r^−2\`, so an untruncated model has \`M(<r)∝r\` and must be matched to a finite virial/tidal environment for cosmological and lensing tests.

**D. Equilibration vs static fit:** equilibrium Bose–Einstein occupation needs an actual physical relaxation mechanism. A rough contact-collision time \`τ=1/[ρ(σ/m)v]\` equals ~2.54 Gyr at \`ρ=10^−24 g cm^−3\`, \`σ/m=1.25 cm² g^−1\`, \`v=100 km/s\`; it scales inversely with cross-section and local density. Collisionless phase mixing does not automatically establish a BE chemical potential.

**E. Statistical inference:** the source ignores baryons in 22 LSB galaxies, excludes two poor-fit systems from the common \`m=86 eV\` scan, and assigns both exceptional systems an additional condensate population; most reported \`ν\` hit the \`ζ(3/2)\` boundary. A robust population constraint needs all galaxies, baryonic mass-to-light nuisance parameters, distance/inclination uncertainties, physical priors and out-of-sample model checks.

**F. Microphysics:** the Bose-star radius bound \`R_c≲83.7 AU\` at \`m=50 eV\`, \`σ/m=1.25 cm²/g\` is a particular self-interacting TF estimate. Neither typical masses, abundance, formation times nor microlensing detectability follow from the radius alone.

## Prioritized go/no-go work
P0. Derive the primordial phase-space distribution, solve linear perturbations, constrain \`P(k)\`/Lyman-alpha and abundance for the **same** mass/couplings; quantify if a viable window exists.

P1. Refit all 22 original galaxies with transparent likelihoods; then all quality-controlled SPARC galaxies with gas, disk and bulge terms; compare cross-validated models and infer a *shared* global mass.

P1. Quantify true kinetic times and condensate-star mass function, including instability, evaporation and host halo response.

P2. Derive finite-mass distribution with escape energy; study phase-space positivity and stability.

P2. Test lensing and wide binaries with extended-lens cross-sections, not point-lens assumptions.

P3. Develop Boltzmann-to-Vlasov/hydrodynamics solver and eventually phase-space / multi-fluid N-body validation.

## Sources worth reading first
2025 condensate review: https://doi.org/10.3389/fspas.2025.1538434
2024 Lyman-alpha: https://doi.org/10.1103/PhysRevD.109.043511
2025 mixed WDM: https://arxiv.org/abs/2504.06367
2023 kinetic Bose-star formation: https://doi.org/10.1103/PhysRevD.108.043535
2024 self-interacting nucleation: https://doi.org/10.1103/PhysRevD.109.016002
2016 SPARC: https://arxiv.org/abs/1606.09251
2026 SPARC degeneracies: https://doi.org/10.1093/mnras/stag2148
2026 Jeans/SIDM: https://doi.org/10.1093/mnras/stag969
2026 astrometric lensing: https://doi.org/10.1088/1475-7516/2026/09/072
