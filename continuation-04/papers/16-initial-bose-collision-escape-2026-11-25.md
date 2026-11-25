# Initial Bose C22 escape creation — 2026-11-25

The equilibrium profiles studied in paper 15 are collisionless. Here evaluate their initial local contact-scattering operator without treating the distribution as a thermal bath. Use nonrelativistic identical bosons, isotropic s-wave total cross section σ, molecular chaos and the Uehling–Uhlenbeck factors. No condensate, GP feedback or gravitational scattering. The initial occupation above the imposed energy cutoff is zero.

Let q=v/(sqrt(2)V0), ψ=(E_t-Φ)/V0², F(s) the occupation at s=ψ-q², zero for s≤0. Incoming q1²+q2²=S; A=S/2; c is their directional cosine; outgoing q±²=A±Bu, B²=A²-q1²q2²c² and u uniform in [-1,1]. For threshold P=ψ+ε, ε≥0, at most one outgoing particle exceeds P because S≤2ψ. This requires S>P and |c|<sqrt[P(S-P)/(q1²q2²)]. The outgoing average of the count times stimulation is (1/B) integral_{ψ+P-S}^{ψ-A+B} [1+F(s)] ds, zero if the upper limit is not greater. For linear lowered BE, F(s)=1/expm1(η-s)-c0, c0=1/expm1(η); primitive uses ln(1-exp(s-η)). For quadratic lowering subtract d0*s, d0=c0(1+c0), in F. The incoming angular average uses dc/2, and relative velocity sqrt(S-2q1q2c). Identical incoming pair counting is 1/2; σ is already the identical-particle total cross section.

Write I(ψ)=integral q² F(ψ-q²)dq, n=D I with D=(sqrt(2)mV0)^3/(2π²ℏ³). The dimensionless collision integral J(ψ,ε) integrates q1²dq1 q2²dq2, incoming angle, relative speed in sqrt(2)V0 units, and the outgoing stimulated count. Integrated birth rate per total halo particle is
  Γ = ρ0 (σ/m) sqrt(2)V0 / (2 Ic² Qt) integral r² J dr.
Here r in core radius units; Ic=I(ψ0), Qt=integral r²ρ/ρ0 dr. Adopt illustrative R=.6 kpc, V0=44.1 km/s, σ/m=1.25 cm²/g, ρ0=V0²/(4πGR²)=6.76485e-21 kg/m³. These are benchmarks, not measured properties of a fitted galaxy. Quantum normalization determines m separately for each profile. For the physical isolated unbound energy E=0, ε=Qt/rt; imposed cutoff ε=0 is a different threshold. Particles crossing the latter can remain bound.

## Refined local numerical results (η=4, critical center, n=48 angular/momentum; nr=160 radial)
| model | threshold | integral J no final stimulation | integral J Bose | rate/Gyr no final stimulation | rate/Gyr Bose |
|---|---|---:|---:|---:|---:|
| linear | cutoff | .0724180814 | .0997269911 | .0019341117 | .0026634666 |
| linear | E=0 | .0127376721 | .0193388316 | .0003401924 | .0005164934 |
| quadratic | cutoff | .0137953565 | .0200287047 | .0004308361 | .0006255068 |
| quadratic | E=0 | .0074655424 | .0112132336 | .0002331527 | .0003501951 |

m=87.5663 and 89.2796 eV respectively (not fixed to 86 eV). n=24/nr=96 agree with refined totals within 3e-6 relative; integrated mass crosscheck. Subcritical linear δ=.2, m=105.4768 eV, gives .0051246533 and .0008614141/Gyr. Near outer boundary F(s)~s^k, q~sqrt(ψ), hence local J~ψ^(2k+7/2) at cutoff, not a singular edge source. Cutoff birth rates differ by factor 4.26 despite similar total masses; physical unbound rates differ by factor 1.47. Bose enhancements 1.377–1.518 at these thresholds.

Code: continuation-04/code/lowered_bose_collision_escape.py. Independent direct CM weak-form check is being evaluated separately. This is an instantaneous source, NOT a measured evaporation rate, a halo lifetime, or a stationary open-halo solution. Actual transport, subsequent rescattering, collective evolution, density changes, and condensate formation are uncomputed. No cosmological relic, CLASS readiness or observational exclusion follows. Local checks only; remote CI not verified.
