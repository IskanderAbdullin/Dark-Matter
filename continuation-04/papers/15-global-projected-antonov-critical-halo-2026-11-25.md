# 15 — Global projected Antonov form for bounded and exactly critical lowered Bose halos
2026-11-25. Independent halo direction of GOAL. Newtonian, isolated, isotropic collisionless model, no baryons, external tides, scalar coherence or collisional Bose terms. This is a quadratic-form/conditional linear-energy result, NOT a proof of nonlinear orbital stability for an unbounded critical DF, a time-dependent spectrum or galaxy fit.

## 1. Existing nonlinear theorem applies to bounded models

Lemou, Méhats, Raphaël, *Orbital stability of spherical galactic models*, Invent. Math. 187 (2012), 145–194, [arXiv:1007.4095, Theorem 1.2, Proposition 3.6](https://arxiv.org/abs/1007.4095), [DOI](https://doi.org/10.1007/s00222-011-0332-9), establish orbital stability (mod translations, including nonspherical perturbations) for continuous nonnegative compactly supported isotropic Newtonian Vlasov–Poisson steady states with strictly energy-decreasing DF and stated energy/L¹/L∞ bounds on the perturbed solution. For both linear and quadratic lowered Bose families, at 0<psi0<eta the occupation is bounded, strictly increasing with binding s=E_t−E in units kT, continuous to zero at s=0. The quadratic form's derivative is positive for s>0. The exterior vacuum potential −GM/r fixes E_t<0 in gauge Phi(infty)=0. Below the minimum physical energy the DF can be monotonically C¹ extended (for example an exponential extension matching value and derivative); this does not affect any occupied phase-space point. Thus for a self-consistent finite-radius bounded lowered model these hypotheses hold. No need to solve individual l=0,1,2 growing modes just to infer the existing theorem. This is a conditional theorem application, not proof that such equilibria form. Exact critical psi0=eta has unbounded occupation and is NOT covered by this theorem as stated.

## 2. Dimensionless global form

Let r be xi, psi relative binding potential, rho=I(psi)/Ic, g=−psi'=Q/r², g'=rho−2g/r. I(psi)=1/2 ∫_0^psi F(s) sqrt(psi−s) ds, with F(s) the dimensionless lowered occupation. W=I'(psi)/Ic is the local susceptibility. Potential perturbation h has energy gradient norm A_l=∫_0^∞[r² h'²+l(l+1)h²]dr. Use normalized angular harmonics. For l>=1 the full energy-shell average vanishes. The generalized projected Antonov form is A_l−B_l, B_l=∫_0^rt r² W h²dr. For l=0 the correct B is instead
B_0=(1/(2Ic))∫_0^psi0 ds F'(s)∫_0^r_s dr r² sqrt(psi−s) [h(r)−bar h(s)]²,
where bar h is the average with weight r² sqrt(psi−s) on the shell. This is the energy-shell projection (not the full orbital E,L projection). For fixed radial support harmonic continuation outside rt gives exterior A=(l+1)rt h(rt)² and minimizes the exterior gradient term. The unprojected l=0 local susceptibility alone is not the dynamically admissible response.

## 3. Direct critical radial bound, including singular endpoint

Adapt the energy-shell Hardy proof in the cited Proposition 3.6 to this dimensionless normalization. Put p=psi−s, w=r² sqrt(p), j(r,s)=∫_0^r w(h−bar h)dr, T=w^−1 d/dr and G=r³p^(3/2). j vanishes at both endpoints. Direct differentiation using g'=rho−2g/r yields
−T²G/G = 3(rho+g/r)/(2r⁴ p²).
Integration of the ground-state transform j=Gv gives ∫w(Tj)² >= (3/2)∫(rho+g/r)j²/(r² p^(3/2)) dr. For smooth h and s strictly interior the boundary term vanishes: j=O(r³) at 0 and O((r_s−r)^(3/2)) at the outer turn, same powers as G.

Let c=1/(2Ic); integrate over a compact interior interval of s with weight c F'(s)≥0, defining B_cut. By integration by parts B_cut=−∫h' (∫c F' j ds)dr. Apply Cauchy–Schwarz in r and then s with weights p^(3/2). Since F(0)=0,
∫_0^psi c F'(s) p^(3/2) ds = (3/2)rho(r).
Hence B_cut² <= A_0 C_cut where C_cut=(3/2)∫dr ds c F'(s) rho j²/(r² p^(3/2)). Meanwhile the Hardy bound gives B_cut >= (3/2)∫ c F'(s)(rho+g/r) j²/(r² p^(3/2)). Monotone rho(r) implies g/r=Q/r³>=rho/3. Thus B_cut >= (4/3) C_cut. Combining and increasing the interior cutoffs to the full domain (nonnegative integrands, monotone convergence) gives

**B_0 <= (3/4) A_0, i.e. A_0−B_0 >= A_0/4.**

This derivation does not require bounded F at s=psi0. At exact criticality psi0=eta, prior central expansion gives W~π sqrt(6)/(4 Ic r); compact support and W∈L^(3/2)(R³) allow continuity of the form on homogeneous H¹ by Sobolev, and extension from smooth test functions. The energy-shell projection reduces the local weighted norm; it is bounded on that space. This provides a quantitative radial projected form bound even for the unbounded critical DF. It is not a bound on the unrestricted energy-Casimir Hessian with arbitrary mass-preserving perturbations.

## 4. Nonradial sectors and linear energy implication

Differentiated Poisson gives −(r²g')'+2g−r²Wg=0, with g>0, g~r/3 at center, g=Qt/r² outside, including critical cusp in density. Ground-state identity for l>=1:
A_l−B_l = ∫_0^∞ r²g²[(h/g)']²dr + [l(l+1)−2]∫_0^∞h²dr.
It follows initially for smooth compact support and extends by continuity. l=1 equality is only h proportional to g (translations); l>=2 strict positivity. Because the compactly supported W~1/r lies in L^(3/2), multiplication defines a relatively compact quadratic form on homogeneous H¹ (truncate singularity, local Rellich, control remainder by Sobolev). Thus there is a positive nonradial coercivity margin modulo translations, though no analytic numerical value for the dipole complement is claimed. For l>=2 one may also use C=max r²W<∞: B_l<=A_1 and B_l<=C∫h² yield B_l/A_l<=C/[C+l(l+1)−2].

For a regular linearized perturbation δf supported in the equilibrium's occupied domain, with finite weighted energy C_f=∫|δf|²/|F_E'| in a consistent Poisson normalization, and zero energy-shell integrals, Poisson and weighted Cauchy give A<=sqrt(C_f B), where B is the projected susceptibility of its δPhi. Thus C_f−A>=0, strictly positive modulo translations (radial C_f−A>=C_f/4). Dynamically accessible δf={generator,f0} has zero shell integrals under valid integrations by parts; a nonzero-frequency growing eigenmode also must preserve linear energy-dependent Casimirs. Under a well-defined linear flow preserving this quadratic energy and those constraints, exponential growth in this admissible domain is ruled out (translation null mode is stationary). Exact critical F_E' singularity requires specifying the domain, integrability and well-posedness; no proof of linear evolution existence in full singular norm, nonlinear orbital stability, damping, actual normal-mode frequencies or collisional stability. This closes the sign question for the projected potential quadratic form, not the entire critical dynamical question.

## 5. Independent local numerical form calculation

Code lowered_antonov_hessian.py uses Bose pole subtraction in the radial quadrature (I singular part sqrt(psi)−sqrt(y) atan sqrt(psi/y), I' singular part atan sqrt(psi/y)/(2sqrt(y)), y=eta−psi), independent Poisson integration; piecewise-linear potential basis r_i=rt (i/n)²; exact exterior harmonic energy, 8-point element quadrature, energy projection integrated by Gauss-Legendre s=psi0(1−u²). Report the largest generalized Rayleigh quotient B/A, NOT a physical mode frequency. Linear/quadratic critical eta4 radii 40.2987556 and 130.1782370; masses 26.3570137 and 26.3193883 agree with prior profile to a few ppm. King depth4 independent regular comparator.

| model | n / energy nodes | radial B/A max | unprojected radial B/A max | l1 first | l1 second | l2 first |
|---|---:|---:|---:|---:|---:|---:|
|linear critical|48/128|.557709|5.583020|.999570|.547914|.407888|
|linear critical|96/256|.558805|5.585870|.999892|.548743|.408000|
|linear critical|144/384|.559009|5.586399|.999952|.548898|.408020|
|quadratic critical|48/128|.549932|5.696460|.998528|.555329|.405523|
|quadratic critical|96/256|.553651|5.706716|.999629|.557961|.405909|
|quadratic critical|144/384|.554353|5.708629|.999835|.558458|.405981|
|King|144/384|.579173|4.990367|.999948|.442977|.425952|

The l1 first approaches one (neutral translation); l1 complement and l2 positive margin. Radial bound .75 satisfied. The local unprojected radial test would give a spurious negative form; not a global Vlasov instability. In subcritical linear sequence delta=eta−psi0=.2,.02,.002,2e−5,2e−7 with 96/256, radial max=.566412,.559685,.558887,.558806,.558805; quadratic sequence .561901,.554665,.553746,.553652,.553651. Numerical endpoint convergence consistent with proof, not a rigorous error bound for exact eigenvalues. Source 4dc3b2f, reproduce python continuation-04/code/lowered_antonov_hessian.py --model linear --delta 0 --n 144 --ne 384, repeat quadratic/king, change delta. Local runs only, no remote CI.

## 6. Remaining work

P0 halo mathematical sign gate is advanced: bounded cases have an existing orbital stability theorem; critical projected form positivity with explicit radial estimate; admissible linear-energy implication conditional on well-posedness. Critical nonlinear theorem needs extension of rearrangement differentiability/compactness beyond bounded distribution; cannot take bounded-theorem limit without uniform hypotheses. For actual observational/kinetic goal, prioritize finite-coupling perturbation theory and collisional relaxation or concrete dynamical orbit response, not more local Jeans tests. No claim that early cosmology produces these profiles, nor that 86 eV is measured.
