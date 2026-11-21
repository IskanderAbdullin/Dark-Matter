# 19 — A hard-lowered Bose halo cannot be an exact CLOSED elastic Bose collision equilibrium

**Source:** Abdullin & Popov, JCAP 11 (2021) 055. This is a direct consequence of the established Uehling–Uhlenbeck / quantum Boltzmann H theorem, applied to our energy-lowered phase-space family. The theorem is not our discovery; the concrete phase-space application and counterexample are new to this project.

## Collisionless versus collisional stationary states

An isotropic f(E) is stationary under collisionless spherical Vlasov–Poisson because E is conserved along stationary orbits. This **does not** imply vanishing quantum Bose 2-to-2 collision operator C22. The latter has gain–loss factors f1 f2 (1+f3)(1+f4) and f3 f4 (1+f1)(1+f2) for admissible energy- and momentum-conserving collisions.

Define the Bose logarithmic collision invariant L(E)=log[(1+f(E))/f(E)]. Detailed balance requires L(E1)+L(E2)=L(E3)+L(E4) whenever E1+E2=E3+E4. Under a sufficiently connected collision kernel with positive scattering, the H theorem requires L(E)=a+bE (allowing momentum-dependent terms for a moving equilibrium). This corresponds exactly to a Bose–Einstein equilibrium occupation.

But the k-th order lowered distribution satisfies f_k(E)~C_k(E_t−E)^k as E approaches its finite escape energy E_t from below. Therefore

\[
\boxed{L_k(E)\sim-k\log(E_t-E)-\log C_k,\quad k\ge1.}
\]

This diverges at finite E_t and is not affine. Thus **no exact hard-lowered f_k is a closed elastic Bose detailed-balance equilibrium** when the collision kernel connects such energies. The distribution can still be perfectly stationary under the *collisionless* equation, and a physically **open** halo with escaping particles or sources/sinks could sustain a quasi-steady truncated distribution. This is not a physical impossibility theorem for Bose dark matter.

## Constructive 2-to-2 witness

At eta=4, take local dimensionless escape depths t1=.05, t2=.35 and t3=t4=.20 (sum .40 on both sides). At a spatial point with psi=1, dimensionless squared momenta are q1²=.95, q2²=.65, q3²=q4²=.80. Suitable 3D scattering angles conserve **both** momentum and energy. Evaluate Delta L=L1+L2−L3−L4:

| k | Delta L |
|---:|---:|
| 0: unlowered thermal BE | 0 |
| 1: C0 escape cutoff | 0.824817 |
| 2: C1 escape cutoff | 1.652053 |
| 3: C2 escape cutoff | 2.478586 |
| 4: C3 escape cutoff | 3.305254 |
| 5: C4 escape cutoff | 4.131894 |
| 6: C5 escape cutoff | 4.958523 |

Every lowered model fails this specific detailed-balance relation; the unlowered BE model obeys it exactly. Under the quantum H theorem, a nonzero imbalance on an open set of allowed collisions prevents a completely closed C22=0 entropy-production equilibrium. This test is **not a numerical full collision rate**, and does not quantify how quickly a weakly scattering halo evolves.

## Scientific consequence

The mathematical finite-mass positivity and classical linear stability literature do **not** establish kinetic equilibration of the same DF. To derive the energy cutoff physically requires open-system escape/accretion, a source/sink term and conservation laws, plus the full Bose C22/C12 exchange if condensates are present. A galactic rotation curve alone cannot determine Bose temperature or origin from a thermal relic.

References: quantum Uehling–Uhlenbeck gain–loss https://www.sciencedirect.com/science/article/abs/pii/S0378437199005506 ; quantum Boltzmann equilibrium and H theorem https://web.stanford.edu/~lexing/eQBE.pdf ; coupled GP/kinetic theory https://doi.org/10.1103/PhysRevD.108.083513 ; source JCAP https://doi.org/10.1088/1475-7516/2021/11/055.
