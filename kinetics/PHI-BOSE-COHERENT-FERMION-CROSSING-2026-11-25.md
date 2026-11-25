# Coherent field coupled to stable charged fermions: exact mass crossing and local Dirac evolution
2026-11-25 UTC. Extends [coherent self-resonance](PHI-BOSE-COHERENT-RESONANCE-2026-11-24.md). Prescribed scalar trajectory, not a coupled fermion-scalar relic solution.

The off-diagonal real Yukawa interaction has instantaneous Dirac mass matrix M=[[M1,yφ],[yφ,M2]], eigenvalues (M1+M2±sqrt((M1−M2)²+4y²φ²))/2. The light signed eigenvalue crosses zero only for |φ|=sqrt(M1 M2)/y=44.72136 GeV for M1=.1, M2=.02 GeV, y=.001. A negative eigenvalue denotes a signed Dirac mass, not negative physical energy. Dark vector U(1) remains conserved; possible fermion production is particle-antiparticle symmetric.

At the previous normalized homogeneous T_SM=20 MeV benchmarks λ=4.698884e−13 and 3.82587e−8, φmax=1.18249 and .182763 GeV, m_light,min=.0199825 and .01999958 GeV. Bounds on |dm_light/dt|/m_light,min² are 2.19e−8 and 2.16e−8; |dθ/dt|/Δrest ≤4.61e−8 and 2.95e−7 with mixing θ=.5 atan(2yφ/(M1−M2)). These exclude a literal zero-mass crossing at this epoch and indicate adiabatic low-p motion, NOT a rigorous bound on all particle production or high harmonics.

For the pure quartic+quadratic adiabatic normalized branch, using exact action from previous report, the threshold corresponds to Acrit=102.902 and 29362.4, acrit=2.25530e−13 and 3.42573e−14. If gS∈[10.75,106.75], visible T at crossing lies in [.346,.743] GeV and [2.277,4.894] GeV respectively. This is a conditional range, not a solved SM thermal history. The pure scalar potential and lack of fermion/fragmentation backreaction can already fail at those high amplitudes.

At a crossing with φdot≠0, |dm_light/dt|=2y²φcrit |φdot|/(M1+M2). The isolated linear mass-crossing approximation has occupation exp(−πp²/|dm/dt|) and characteristic p*=sqrt(|dm/dt|). At amplitude 1.1φcrit, p*=.4483 and 7.572 MeV. Checks of the approximation: |m¨|/|m˙|^(3/2)=.052 and .878; consecutive crossing time Δt p*=28.42 and 1.68. For 1.001φcrit these become 2.045/34.54 and .975/.0577. Thus the larger-λ and near-touching cases are not reliably independent linear crossings. A formal isolated-crossing particle+antiparticle count 4 p*³/(8π³) cannot be multiplied by number of oscillations or interpreted as a relic.

## Exact prescribed-background Dirac check
For a fixed momentum/helicity, evolve the Hermitian 4×4 Hamiltonian H=[[M,pI],[pI,−M]] with φ=φmax cn(z|A²/(1+2A²)), z=m sqrt(1+2A²)t. Start from instantaneous negative-energy subspace at φ=0, propagate through the positive lobe to next φ=0, project onto positive-energy eigenvectors (ordered light, heavy). This is a precisely defined finite-time initial vacuum, not necessarily physical cosmological vacuum. Each occupation per spin and particle is bounded by 1; antiparticle count is equal. Two crossings interfere.

For λ=3.82587e−8 (momenta p/p*):
| φmax/φcrit | 0 | .25 | .5 | 1 | 2 |
|---|---:|---:|---:|---:|---:|
|1.001|0|.028944|.108605|.336047|.472601|
|1.1|0|.491532|.988123|.185122|.00069854|
|2|0|.412118|.409889|.057319|.00017633|

Heavy occupation at 1.1, p=p*, is .000141810. Largest norm error <5e−9, zero y control <2e−32, refined rtol 1e−11 changes occupations at 1.1,p* by <2e−13. At λ=4.698884e−13, 1.1φcrit,p=.25p*, one slower local point gives light .585880, heavy 3.27e−9, norm defect 1.29e−6; full small-λ grid was interrupted by local 20 s execution limit and is NOT claimed converged. The fast light occupation is not monotone in p or amplitude and the p=0 cancellation is consistent with two mass sign changes, not a contradiction of isolated single-crossing probability.

Reproduce python kinetics/phi_bose_coherent_fermion_crossing.py (default excludes expensive weak-λ Dirac integration). Local tests only, no new remote CI verified. Need true initial state, full cycle phase/Pauli evolution, UV and momentum convergence, expanding history, scalar self-fragmentation and backreaction, and eventual stable χ2 count including annihilation. CLASS-ready=false.
