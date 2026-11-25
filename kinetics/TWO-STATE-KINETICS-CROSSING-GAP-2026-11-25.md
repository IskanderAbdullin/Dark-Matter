# Two-state activation: coherent crossing bounds, kinetic halo conversion, and gap tradeoff

2026-11-25. Continuation of 177b6f3. All results are conditional on the soft mixing EFT, not a complete microscopic model or observational exclusion.

## Coherent flavor transfer

For a positive occupation matrix, |rho_AB| <= sqrt(rho_AA rho_BB). With a diagonal A source/loss and dephasing, no direct B source or off-diagonal bath interaction, total B/N_injected <= (integral |H_AB| dt)^2. This follows by integrating d sqrt(B)/dt <= |H_AB| sqrt(N_injected); use Cauchy-Schwarz for momentum-dependent occupations with redistribution that preserves total trace. Arbitrary diagonal detuning (including crossings) cannot violate this inequality. For the positive energy square root of a two-state mass-squared matrix with B diagonal mass >=m, |H_AB| <= |delta chi/v|/sqrt(p^2+m^2). No assumption about a linear resonance or instantaneous decoherence is required.

With r=1+z, choose r_cut=1e14 (after production at reference parent x well above 80). Require H >= .5 H0 sqrt(Omega_r) r^2, scalar trace bounded by twice the final comoving number, early field seed 1e-8 with zero derivative, and no separate B source. The homogeneous seed envelope, extended above r=1e16 with a conservative absolute-force integral, is verified up to r_cut by the specified roll. Worst q=0 coherent probabilities before r_cut are 7.45e-18, 2.24e-17, 9.32e-17 for .1,.3,1.25 cm2/g; at q_peak they are 9.36e-28, 2.81e-27, 1.17e-26.

After cutoff assume the vacuum diagonal mass difference, no scattering or source, and follow the instantaneous eigenbasis. Before z=1000 the integrated absolute mixing-angle variation is O(1e-8). Thereafter integration by parts bounds the transition amplitude using endpoint |g/W| + total variation of g/W + integral |g^2/W| (g=dtheta/dln a, W=Delta E/(hbar H)). Refine 4k/16k/64k samples; late amplitude estimates <=7.9e-22. Combined heavy eigenstate probability estimates <=1.6e-15. These estimates are for a smooth prescribed homogeneous field history and vacuum postcut Hamiltonian; not interval-arithmetic bounds, not a full finite-temperature QKE or a bound on collision-induced B. A positive quasiparticle Hamiltonian and disappearance of other off-diagonal interactions must be justified independently. Source at early high T and B self scattering remain open.

## All six pair channels in a local Bose halo box

Implement LL,LH,HH pair types with six reversible pair channels, equal inertial masses and retained internal energy d=gap/T. With K_i=q_i^2, the angular phase-space factor is length of overlap of P intervals [|q1-q2|,q1+q2] and [|q3-q4|,q3+q4]. Exact NR relative-speed/final-flux factor cancels the angular Jacobians, leaving dk1 dk2 dk3 dP/(32 pi^4). Include pair symmetry, final distinguishability and reverse orientation factors. Positive nodal interpolation gives exact represented number/energy and detailed balance, and stable geometric Bose net flux. Numba is optional; NumPy fallback included.

For sigma_LL=.3 cm2/g, V0=44.1 km/s, nominal full mixing gap/T=.037284, initial near-critical L alpha=1e-6 and H empty, tau=n_critical sigma_LL vT t=3.9564 for 10 Gyr. Three completed runs: legacy n48, focused n64, focused n88. Initial H gain 16.96,16.89,16.47 per tau versus independent angular integral 15.90. The fine initial discrepancy remains 3.6%, partly initial remesh error. By tau~.055-.058 H fraction reaches .40. Final H fractions .463065,.462771,.461147; independently solved continuum initially critical box .458248. Discrete N/E-matched BE equilibrium H fractions .462771/.461147 for n64/88; final weighted L1 deviations 1.86e-10/1.02e-11. Number and full kinetic+internal energy residuals <2.5e-15; near-empty negative numerical occupations have negligible weighted fraction (n64 final 8.1e-12, n88 <1e-18). Detailed balance <2e-15, Jacobian <1e-9, positive entropy checks. Fixed-volume local chemistry is therefore numerically established within this discretization, rather than inferred solely from a frozen initial source. No gravitational dynamics, escape, condensate, number-changing, spatial transport or time-varying mixing in this box. The rapidly reacting initial state invalidates transferring an unmodified one-species equilibrium halo to this EFT.

## Can the channels instead be kinematically closed?

Optimize over all u in (0,.5), Delta >= lambda_LL*m^2/(32 pi^2*u^2), delta=Delta sqrt(u(1-u))/(1-2u), X>=delta^2/(16 pi^2), and Xcap=min(pi^2/9 N_eta T_eta^4,.04 rho_total). For MB relative velocities, exact threshold factor G(d)=d/2 exp(-d) kve(1,d/2), equivalent to integral sqrt(y(y+d)) exp(-y-d) dy. Ratios LL->LH / LL->LL =2(1-u)/u G(d), LL->HH / LL->LL=((1-u)/u)^2 G(2d). Scan 120 cases including .1,.3,1.25, V0=10,20,44.1,100, z_nom=3,5,8,10,20 and Delta N_eff=.1,.3.

F583-like V0=44.1, sigma=.3, nominal z8, Delta N_eff=.1: maximum permitted gap/T=.11885, minimum total conversion/LL=3.54; to reduce conversion below .1 and two-state critical-density shift below 5%, optimized required gap/T=8.361, u=.012976, X/Xcap=584.16. Even z20, Delta N_eff=.3 needs X/Xcap=6.57, and earlier activation must separately pass the cosmological kinetic test. This strongly disfavors simply closing all inelastic channels within the stipulated illustrative loop and thermal caps. These are model-dependent naturalness screens, not a universal no-go; fine tuning or different UV/thermal origin may change them.

## Sources and next priorities

- Density-matrix production and dephasing: https://doi.org/10.1103/PhysRevD.76.085011 ; https://doi.org/10.1103/PhysRevD.97.115017 .
- Inelastic SIDM halo simulations (different masses/cross sections; cannot transfer constraints): https://doi.org/10.1093/mnras/staa3315 ; https://doi.org/10.1093/mnras/stae1270 .
- Discrete quantum Boltzmann equilibrium: https://doi.org/10.1051/cocv/2021079 ; stiff quantum collision numerics https://doi.org/10.1051/M2AN/2011051 .
- Two pNGB UV example with independent soft breaking: https://arxiv.org/abs/2604.12435 (not a completion of our EFT).
- Independent 2026 nonthermal self-thermalization analysis: https://arxiv.org/abs/2610.21817 (different production and approximations; no direct transfer of mass limits).

Next: two-state halo rotation-curve comparison and finite gravitational halo transport, plus early direct B-collision bootstrap and UV protection. Raw evidence: results/two-state-kinetics-crossing-gap-2026-11-25.json.
