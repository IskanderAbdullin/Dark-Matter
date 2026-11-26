# Strict on-shell no-go for a single evolving GLOBAL Bose Tφ(t),μφ(t) in the finite-P full Dirac Yukawa network (2026-11-24)

## Physical question and result

The full *vacuum-mass, hard parent-cutoff* four-channel Dirac Yukawa χ1,antiχ1,χ2,antiχ2+real φ collision kernel has a positive scalar first energy cell below the true physical pair-production / decay floor. At P=Ppair=2GeV the first scalar energy cell runs mφ=86eV to 99.9975186keV. **Every** actual on-shell χ1↔χ2+φ and χ+antiχ↔2φ event is outside this interval; its physical Bernstein-P2 number/energy moments are exactly collision-frozen, no matter what thermal/chemical occupancies are assumed.

Suppose one attempts to evolve ALL φ energy modes using a **single Bose-Einstein form**
  bφ(E,t)=[exp((E−μφ(t))/Tφ(t))−1]^−1,
with positive subcritical μφ<mφ, and allows both T and μ to vary in time. This *two-parameter thermal manifold* is NOT invariant under the physical restricted collision operator whenever another active energy cell has nonzero source. This is a rigorous scoped closure obstruction, NOT a physical 86-eV Bose dark matter/condensation no-go.

## Proof using physical cell moments and a full-rank 3×2 tangent

In the frozen finite-width cell [mφ,K], let Ba(E), a=0,1,2 be nonnegative Bernstein-P2 weights with ΣBa=1 and ΣEa Ba=E, and Dφ(E)dE=p²dp/(2π²). The three physically integrated thermal cell moments are

  Ma(T,μ)=∫[m,K] Dφ(E)Ba(E)bφ(E;T,μ)dE.

A genuinely frozen collision cell requires dMa/dt=0 for ALL a=0,1,2. The Jacobian of these three moments with respect to (μ,T), at positive subcritical Bose distribution, is

  A_{aμ}=1/T ∫Dφ Ba b(1+b),
  A_{aT}=1/T² ∫Dφ Ba(E−μ)b(1+b).

Both columns belong to the physical positive DOS susceptibility map, not invented capacities. The columns are LINEARLY INDEPENDENT because a finite populated energy interval has nonzero variance: the B0- and B2-weighted mean energies differ (B2/B0 monotonically grows in E). Thus rank A=2, so A(μdot,Tdot)^T=0 enforces μdot=Tdot=0. The active scalar source at a distinct energy is generically NONZERO for nonequilibrium chemical affinities, so it cannot simultaneously equal a thermal manifold tangent with both parameters stationary. No homogeneous BE(T(t),μ(t)) ansatz can exactly represent the full collision-time scalar dynamics of this hard-P network.

Independent physical full 5-population on-shell example at P2, common T1GeV, χ and antiχ chemical potentials all zero, **μφ=−0.2GeV** (admissible BE):
- First φ cell all three on-shell moment rates: **exactly zero**.
- Sum of active scalar energy-cell number sources: **+6.740027458965697e−12 GeV⁴** (nonzero).
- Total microscopic H=**+1.34800549179314e−12 GeV⁴** (positive).
- Singular values of the first-cell physical (μφ,Tφ) moment Jacobian: **2.90710083983365e−13 and 2.21404818158768e−18** (units depend on column units); **full rank TWO**, ratio ≈7.615×10⁻⁶, independently verified with positive 168-node physical Bose IR quadrature.
- Net Dirac U1 and shared on-shell total E event defects vanish at double precision.

Independent tests confirm the exact collision freeze and nonzero active source for hard P2/8/16, the physical two-parameter Jacobian rank, and that as μφ→mφ−δ with δ=1e−12,1e−14,1e−16GeV the **frozen** chemical response grows, while its actual microscopic source remains zero. Tests run in focused CI #45426807185, completed success.

## What this DOES and DOES NOT imply

This theorem does not assert that a physical nonthermal scalar cannot be dark matter or become condensate. It only invalidates an *exact*, globally uniform BE two-parameter closure for the limited vacuum Yukawa network at a STRICT finite parent cutoff. It strengthens the existing numerical warning that a hard cutoff creates disconnected IR basis DOFs: a thermal ansatz changes those modes even though the physical collision operator cannot.

In the true early Universe there is no literal hard cutoff in an infinite thermal FD tail. For the same GeV-temperature benchmark, near-rest scalar production is bounded by 10^(−23093)GeV⁴ for 86–86.1eV, so the physical continuum IR may instead be **astronomically slow, not exactly frozen**; this is a different timescale question. Additional elastic SM channels, full 2↔2, finite-temperature masses, expansion in q=a p and a true zero-mode would change the support. The no-go is explicitly not universal under those interactions.

**Important consequence for existing 5-field MB chemical relaxation eigenvalue pilots:** if those pilots constrain all 5 populations to equilibrium BE/FD shapes with a single μ,T per species, they characterize only a *projected thermodynamic response*, NOT the exact strong kinetic relaxation of all physical energy modes. The real 5-population energy-resolved quantum FD/BE solver is required before promoting a chemical rate into physical cosmic freezeout or BEC.

Reproduce:
- [Exact no-go / susceptibility numerical reference](phi_bose_frozen_ir_thermal_manifold_no_go.py)
- [Four scientific tests](../tests/test_phi_bose_frozen_ir_thermal_manifold_no_go.py)
- [Physical 5-field event solver](phi_bose_five_dirac_full_yukawa_implicit.py)
- [Machine-readable values](results/phi-bose-frozen-ir-global-BE-no-go-2026-11-24.json).

**NEXT_ACTION:** retain energy-resolved scalar fφ(q,t) in all future SM+FRW solvers; test how tiny but nonzero true thermal tails or additional scattering lift the disconnected modes and whether relaxation is slower than H. Do not enforce a single BE(T,μ) for all modes unless full kinetic detailed balance has been quantitatively demonstrated. No physical relic/BEC/CLASS conclusion from this pure closure theorem.
