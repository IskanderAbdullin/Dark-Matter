# Coherent positive-quartic alternative: adiabatic background and exact self-resonance
2026-11-24 UTC. Independent conditional initial history, not a production solution.

The earlier high-entropy thermal and low-reheat particle-injection gates do not apply to an independently prepared low-entropy coherent field. Take a real scalar, V=m²φ²/2+λφ⁴/24, m=86 eV, independently matched positive λ. Ignore interactions with fermions and higher operators for this specific stage. Normalize today's homogeneous energy to Ωh²=.120, not the actual observed abundance from a solved source.

Set x=φ sqrt(λ/(12m²)), τ=mt, amplitude A, e=A²+A⁴ and energy ρ=6m⁴e/λ. Then x''+x+2x³=0. Define I0=∫dθ/[1+A²(1+sin²θ)]^1/2 and Ic=∫cos²θ[1+A²(1+sin²θ)]^1/2 dθ on 0..π/2. Exact oscillator action J=(6m³/λ)(4A²Ic/π), conserved J a³ to adiabatic accuracy. ω=mπ/(2I0), averaged P=ωJ−ρ; virial P=<V4>. Thus w from weighted sin⁴ integral and c_a²=dP/dρ=dlnω/dlnJ=−dlnI0/dlnA/(2+dlnIc/dlnA). Limits w→0 and 1/3, c_a²→0 and 1/3. This is the homogeneous adiabatic branch; its perturbations can invalidate it.

Exact solution x=A cn(z|p), z=τ sqrt(1+2A²), p=A²/(1+2A²). For k_phys/m=κ, the fluctuation equation is δ''_z+[h−6p sn²(z|p)]δ=0, h=(κ²+1+6A²)/(1+2A²). Lamé n=2 eigenvalues: 2(1+p)±2sqrt(1−p+p²), 1+p, 1+4p, 4+p. Physical h≥1+4p, only open positive-momentum gap is 4+p<h<2(1+p)+2sqrt(1−p+p²). Exactly κ_lo²=3(1+A²); width Δκ²=(1+2A²)3p²/[2sqrt(1−p+p²)+2−p]. Compute Floquet monodromy over 2K(p), normalize determinant, μ/m=sqrt(1+2A²) acosh(|tr|/2)/(2K). Analytic gap is not replaced by noisy near-zero numerical μ. Very narrow gaps at late low A can be below double-precision resolution.

Local independent checks: direct integration of original nonlinear oscillator plus perturbation matrix (not elliptic-function coefficients) reproduces center exponent and w for A=.05,.1,.3,1,3,10,100; maximum sampled relative center-exponent error ~2.72e−7; pressure differences <2e−15, |det−1|<3.1e−14. Outside sampled gap at −.2 and 1.2 gap widths trace stable; tiny edge traces are not evidence of extra bands. Quarter period agrees with elliptic identity.

Use SM entropy normalization a=Tγ0/Tγ(gS0/gS)^(1/3), Tγ0=2.34865e−13 GeV, gS0=3.91. Illustrative constant gS=gρ=10.75 above MeV and asymptotic post-annihilation gS=3.91, gρ=3.36, Tν/Tγ=(4/11)^(1/3). Not full BBN calculation. H includes the conditional coherent density plus radiation.

| λ (independent positive examples) | A at 20 MeV | w | H/ω | ρ/ρν,1 | μ/H at gap center | q today (eV) |
|---|---:|---:|---:|---:|---:|---:|
|4.69888e−13|2.72085|.296782|6.07e−16|4.717e−7|5.866e13|3.732e−9|
|3.82587e−8|119.995|.333312|1.44e−17|1.931e−5|2.944e15|1.555e−7|

The first λ is an illustrative magnitude of the previous finite one-loop quartic correction, NOT a prediction of positive physical λ; the second corresponds to the prior illustrative σ/m=1.25 cm²/g. At 20 MeV both backgrounds are rapidly oscillatory and exhibit catastrophic linear instability if seeds exist. A mode initialized at center and redshifted k_phys∝a^-1 exits the band at T≈15.319 and .514 MeV in the two constant-g toy histories; integrated linear exponents ≈1.399e13 and 9.505e16 (65/129-grid relative changes .023% and .168%). These enormous numbers are only formal linear extrapolations; backreaction necessarily invalidates the extrapolation once fluctuations grow. No prediction of final fragmentation, spectrum, abundance, or gravitational observables. The center comoving q is below the prior chosen 1.14e−6 eV rms diagnostic for both examples; fragmentation need not imply hot DM. Redistribution and self-number-changing can alter this conclusion. The ratio ρ/ρν,1 is just an energy-density diagnostic, not a constant ΔNeff when w varies. Low ratio does not certify BBN, and actual stable heavy fermion production still matters.

References: Turner, Phys Rev D 28 1243 (1983), https://doi.org/10.1103/PhysRevD.28.1243 ; Massó et al., action-angle variables https://arxiv.org/abs/astro-ph/0501381 ; Greene et al., self-resonance in quartic preheating https://arxiv.org/abs/hep-ph/9705347 ; Lamé eigenvalues https://dlmf.nist.gov/29.3 . The finite-amplitude reduction and numerical values here are explicitly derived, not observational fits.

Next: physical scalar number-changing and fragmentation/transfer to fermions, or independent halo profile stability. Do not apply linear growth to nonlinear saturation. CLASS-ready=false.
