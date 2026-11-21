# Liouville, Casimir invariants and the nontrivial origin of galactic Bose–Einstein statistics

## The exact collisionless statement
For a sufficiently regular fine-grained distribution \(f(\mathbf x,\mathbf p,t)\) solving the time-dependent collisionless Vlasov equation
\[
\partial_t f+\frac{\mathbf p}{m}\cdot\nabla_x f-m\nabla\Phi\cdot\nabla_p f=0,
\]
the Hamiltonian phase flow preserves phase-space volume and \(Df/Dt=0\) along characteristics. With appropriate boundary conditions, for every suitably regular \(\mathcal C\),
\[
\boxed{I_\mathcal C[f](t)=\int d^3x\,d^3p\,\mathcal C(f)=\text{constant}}
\]
(including particle number and infinitely many other Casimir-like moments). The fine-grained histogram of occupancy values does not spontaneously turn into an arbitrary single-temperature Bose–Einstein occupancy just because the potential virializes.

## What phase mixing can and cannot do
Coarse-grained \(\bar f\) may change under violent relaxation while fine-grained \(f\) conserves its Casimirs; different incomplete mixing histories can lead to quasistationary functions of energy. Lynden–Bell statistical mixing maximizes a combinatorial entropy conditioned by occupied phase-level volumes, and is **not equivalent** to quantum Bose statistics of identical particles. Even though both use exponential-looking equilibrium functions, their entropy constraints and exclusion/occupation measures differ.

This is **not** a theorem that Bose–Einstein-shaped stationary \(f(E)\) is impossible in collisionless gravity. Such a function can be an exactly stationary Vlasov distribution if prescribed initially, or can be approached in coarse graining in special cases. The rigorous statement is that collisionless virialization by itself does **not derive a universal BE temperature and chemical potential**, and cannot generically discard all initial Casimir constraints.

## Why it matters for the source JCAP model
JCAP 2021 fits a late spatial density
\[
n(r)=\lambda_{\rm th}^{-3}\operatorname{Li}_{3/2}(z_0e^{-[\Phi(r)-\Phi(0)]/V_0^2})
\]
and identifies its physical mass parameter through \(n_0\lambda_{\rm th}^3=\nu\). If \(f(E)\) is merely a stationary shape chosen to fit \(V(r)\), the particle mass identification presumes actual **quantum cell normalization and population occupancy**, not merely an arbitrary classical two-parameter King-like \(f(E)\). A theoretically defensible particle-mass inference needs a model of phase-space occupation normalization and its evolution from cosmological initial data.

Key distinction:
- equilibrium of quantum collision term \(C_{22}[f]\to0\) permits BE under detailed balance and number conservation, but establishing the time to approach it requires rates;
- collisionless Vlasov dynamics can be stationary for \(f(E)\) without achieving thermal distribution or erasing fine-grained invariants;
- a condensate requires additional coherent degree of freedom and exchange; Lynden–Bell coarse graining should not be called Bose-condensation.

## Conditional phase-space feasibility check
If initial cold population is approximately smooth in a finite momentum width, its phase-space maximum \(f_{\rm init,max}\) is preserved by collisionless Vlasov. A proposed final fine-grained equilibrium with \(f_{\rm BE}(p\to0)\to\infty\) cannot be created from a strictly bounded initial \(f_{\rm init,max}\) by an *exact* collisionless map alone. But:
1. the divergence is an idealized pointwise infrared limit at exactly \(z=1\); finite resolution, weak interactions and finite-size cutoffs may regularize it;
2. coherent initially cold bosonic populations may have large or singular occupation, so cannot assume a finite initial maximum without specifying \(f(q)\);
3. collisional scattering, wave coherence, quantum kinetic noise and condensate exchange invalidate the simple Vlasov-only premise.
Thus it is an explicit **if-then necessary-condition test**, not an unconditional exclusion of the JCAP gas.

## New falsifiable route
1. Specify cold-production initial \(f(q)\) for Yin / Sakurai–Yin scenario and calculate \(f_{\max}\) for a physically meaningful finite-resolution and phase-space cell definition.
2. Construct a collisionless Vlasov evolution for a galaxy with time-dependent potential, test conservation of \(\int C(f)\), and compare coarse-grained \(\bar f(E,L)\) with ideal BE and Lynden–Bell based candidates.
3. If the desired BE limit requires phase-space population outside the reachable collisionless mixing hull, introduce **explicit** contact collision operator or coherent condensate exchange; calculate its rate and conservation.
4. Avoid identifying an apparent BE fit of \(V(r)\) as a particle *statistics measurement* without independently constraining velocity moments or formation process.

## Sources
Original source: Abdullin & Popov JCAP2021 055.
Chavanis, *Kinetic theory of collisionless relaxation for systems with long-range interactions*, https://doi.org/10.1016/j.physa.2022.128089.
2012 entropy in collisionless systems, https://doi.org/10.1088/0004-637x/748/2/144.
Proukakis et al. PRD108.083513, https://doi.org/10.1103/PhysRevD.108.083513.
Proukakis et al. PRD111.023505, https://doi.org/10.1103/PhysRevD.111.023505.
