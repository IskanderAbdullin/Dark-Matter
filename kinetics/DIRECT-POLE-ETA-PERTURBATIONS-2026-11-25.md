# Eta and trigger perturbations in the direct cubic branch

## Frozen collisionless response

The post-roll background includes repeated intervals of negative effective trigger curvature. On two choices of redshift, momentum and wavenumber sampling, min(m_eff²/H²)≈−1.74e4 and peak local γ/H≈130. At low k the integrated positive local root is 13.3–13.9, almost entirely reproduced by the unforced scalar curvature. Restricting to k/(aH)≥10, γ/H≥5 and γ≥3|d ln m_eta/dt| gives an illustrative local integral ≈6.0 at k≈.003/Mpc. The maximum change in the integrated root from collisionless η feedback is 1.14–1.19 over the sampled k values. These integrals add positive instantaneous roots, even though the stable intervals can undo phase growth. They are diagnostics rather than true accumulated mode amplification.

For η with E(p)=sqrt(p²+g²χ²), perturbations have δE=(gm/E)δχ in the positive-χ patch. At a frozen snapshot the monopole Vlasov response contains 1−atan(y)/y with y=k_phys v/γ. The scalar has a positive direct mass derivative D=Ng²/(2π²)·∫p⁴f/E³ dp and a negative collisionless feedback S=Ng²m²/(2π²)·∫(p/E)(−df/dp)[1−atan(y)/y]dp. The local dispersion is γ²+3Hγ+k_phys²+Vχχ+D−S=0. The finite γ floor is 1e−5 H. This diagnostic avoids treating the formal infrared-sensitive static susceptibility of the decoupled Bose distribution as a measured instantaneous rate. Direct angular integration and the pressureless limit independently check the kernel normalization and signs.

The local framework is related to mass-varying-species calculations, with the actual η Bose spectrum and the oscillating χ potential substituted here. See [Afshordi, Zaldarriaga & Kohri](https://doi.org/10.1103/PhysRevD.72.065024), [code](direct_pole_eta_stability.py), and [two-resolution summaries and trajectories](results/direct-pole-eta-snapshots-2026-11-25.json). Metric feedback, time-dependent streaming and nonlinear evolution remain outside this check.

## Expanding metric-free phase evolution

The linear calculation with time-dependent χ curvature and collisionless η response finds strong amplification around k=.02–.03/Mpc for the benchmark trigger. The fixed homogeneous metric is used for the Liouville and scalar equations. The η line responds to perturbations through δE=(gm/E)δχ and follows the full momentum and angular transport; no pressure closure is imposed. At k=.03/Mpc, the direct χ growth exists also with directional η feedback switched off. This control helps identify the time-dependent trigger as a source of amplification. It does not isolate the fraction due to conventional parametric resonance versus tachyonic intervals.

The first numerical experiment used finite-angle quadrature for the freely streaming η seed. At .03/Mpc it returned η gain 3.09e12 with 12 angular nodes and 8.89e10 with 20, along with opposite phases. That result fails angular convergence. The revised code transports the initial contribution analytically in angle as q·sin(Ψ)/Ψ in the monopole, where Ψ=k∫v dτ includes the changing mass. Only the induced response uses angular quadrature. For an η temperature seed imposed at z=8, the revised gain is 8.90e10 and 8.89e10 on two choices of background, radial, angular, and ODE resolution.

For a metric-free η temperature seed imposed at 1+z=1e6, tested also with an independently propagated χ seed, the predicted k=.03/Mpc η number contrast per unit η temperature seed is about 4.24e10. The analytic freely streaming monopole by the coupling switch is only −0.00687 per unit seed. The factor persists after increasing radial and angular resolution and after moving the coupling switch earlier. The eta mass divided by its momentum-temperature scale is below 2e-9 before this switch on the reference path, and the early scalar equation directly includes forcing from the ballistic seed. An illustrative 1e-5 initial temperature shift reaches |δn_eta/n_eta|=1 near z=2.03 in the metric-free linear pilot. Shifting the background seed χ_initial by a factor of .5 or 2 still leads to this threshold at z≈2.02–2.10 for k=.03/Mpc. Extrapolated late linear gains above this threshold provide a sensitivity measure only.

The instability depends on actual early perturbations. The imposed temperature mode is a diagnostic initial condition without metric forcing, and the calculation does not establish its primordial amplitude, a GR adiabatic spectrum, the nonlinear η distribution, or an observational exclusion. Treat the available smooth-component BAO/growth analysis as conditional. Prioritize a metric-coupled Einstein–Boltzmann calculation with physically sourced η and χ perturbations, followed by nonlinear evolution if modes become large. Revisit trigger designs without deep negative curvature before refitting cosmological parameters.

Reproduce the refined early-seed case from the repository root with `python kinetics/direct_pole_eta_evolution.py --mesh 320 --radial 10 --angle 24 --rtol 7e-8 --background-samples 42000 --primordial-eta --primordial-chi --k .03 --output /tmp/direct-pole-eta.json`. For the data including the rejected finite-angle calculation see [the time-dependent mode results](results/direct-pole-eta-evolution-2026-11-25.json).

## Stress-energy source for the metric calculation

Use conformal Newtonian gauge with time potential ψ and spatial potential φ. Let ε²=q²+a²m², α=∂χ ln m and J=α(ρ_eta−3P_eta). The background exchange is ρ_eta′+3ℋ(ρ_eta+P_eta)=χ′J. The homogeneous scalar has the opposite exchange, with χ″+2ℋχ′+a²(Vχ+J)=0.

With primes on f0 indicating momentum derivatives, the linear η equation is

$$
\partial_\tau\delta f + i k \mu (q/\epsilon)\delta f
 = -q f_0' \partial_\tau\phi
 + i k\mu f_0'\left[\epsilon\psi +
 \frac{a^2m m_{,\chi}}{\epsilon}\delta\chi\right].
$$

The scalar perturbation has a source χ′(ψ′+3φ′)−2a²(Vχ+J)ψ−a²δJ. The mass change also contributes directly to the stress integrands, with δε=a²m m_,χ δχ/ε in the energy weight and δ(1/ε)=−δε/ε² in the pressure weight. These pieces are required in the Einstein sources even at fixed δf. Evaluate m²=(m0χ/v)² and m m_,χ=m0²χ/v² directly so the mass-gradient source stays finite as χ crosses zero.

The [stress audit](direct_pole_eta_stress.py) verifies the trace, exchange, J and its curvature at ten points covering the roll. Maximum relative discrepancies were 2.1e−16 for the trace, 1.1e−9 for the exchange, 5.4e−9 for J and 1.6e−8 for its curvature. These are moment checks before implementing the metric hierarchy. The public CLASS ncdm hierarchy supplies a constant-mass control, and the scalar-gradient term must be added with the same stress-energy convention. In the CLASS multipole convention this adds −k a² m m_,χ (d ln f0/d ln q) δχ/(3qε) to the η dipole. The existing code uses a fixed mass in ε and in its stress integrals. See [Ichiki and Keum's derivation](https://arxiv.org/abs/0705.2134), [the existing CLASS hierarchy](https://github.com/lesgourg/class_public/blob/1542f2c9bf5b7d3df598ab62ec7b3945e72f8525/source/perturbations.c#L11072), and [the stress results](results/direct-pole-eta-stress-2026-11-25.json).
