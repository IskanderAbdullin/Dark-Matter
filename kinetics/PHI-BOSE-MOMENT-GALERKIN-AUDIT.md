# IR-safe Bose/Fermi inverse-fugacity moment Galerkin: physical p-integrated number, energy and Yukawa H (2026-11-22)

## Scientific question, competing alternatives and prespecified gates

The prior four-population χ1,χ2,real φ,one-lepton pilot has on-shell Yukawa detailed balance, quantum entropy and exact discrete charges, but **incorrect φ scalar number**: multiplying BE(E=mφ) at the rest-energy node by the full first radial-cell capacity yields a false condensate/abundance. At mφ=86e−9 GeV and illustrative T=1 GeV, prior N16/64/256 values exceed the true excited thermal gas by 87199×/1363×/22×. Simply swapping the rest-node BE value for its *cell average* changes the inverse fugacity in quantum F/R and can violate exact thermal chemical equilibrium. No physical Ω or BEC follows from the earlier bad counts.

Falsifiable target: a single **physical p-space integration + quantum inverse-fugacity interpolation** must yield the true thermal N and U for mφ=86e−9 at 3 grids, the massless ζ(3) limit, separate χ1/χ2/φ Yukawa N/E/H, and chemical fixed point without clipping or an invented condensate. A counterexample must show why the tempting first-cell average shortcut is not a valid collision replacement. This is an **instantaneous moment** reconstruction, not a validated nonlinear time solver, thermal plasma or cosmological relic.

## Mathematical construction and independent proof

Take a physical energy domain E∈[m,Emax], nonnegative piecewise-linear hat functions h_i(E) with

    Σ_i h_i(E) = 1,             Σ_i E_i h_i(E) = E.

At each species' energy nodes store the **inverse fugacity** η_i; at physical energy E reconstruct η(E)=Σ_i h_i(E) η_i. Then f(E)=1/(exp η(E)+1) for fermions and b(E)=1/(exp η(E)−1) for a real boson, demanding ηφ(E)>0 on the positive-measure support (massless critical ηφ(0)=0 allowed at the p=0 point of zero measure). The true physical radial measure is dΠ=g p²dp/(2π²), E=√(p²+m²). Define **moment coordinates**, not cell point occupancies,

    M_i(η) = ∫ dΠ h_i(E) n(E;η).

The hat identities give N=ΣM_i and U=ΣE_i M_i *on exactly the same quadrature* for any Fermi/Bose f. The occupation is not multiplied by a fictitious first-bin rest-energy volume. Thermal η(E)=(E−μ)/T is affine and thus represented exactly at any energy-grid spacing; at μ=0 even massless BE b(E) diverges only at p=0, a zero-measure point in the excited gas (no δ³(p) condensate added).

The thermodynamic entropy has derivative ds/dn = η. Therefore dS=Σ_i η_i dM_i. The inverse-fugacity-to-moment Jacobian is

    ∂M_i/∂η_j = −∫ dΠ h_i(E)h_j(E)n(E)[1±n(E)]

with plus for bosons and minus for fermions. It is **symmetric negative definite** for a finite strictly admissible field with independent hats. Numerical eigenvalue checks and selected independent finite differences confirm the implemented Jacobian. Its inversion for finite-time evolution, especially near a true Bose critical point, is NOT demonstrated here.

For each physically on-shell moving Yukawa event χ1(E1) ↔ χ2(E2)+real φ(Eφ), E1=E2+Eφ, use the same energy hats for quantum occupations and signed moment deposition. With any positive event weight W and J=F−R,

    F=f1(E1)[1−f2(E2)][1+bφ(Eφ)]
    R=[1−f1(E1)]f2(E2)bφ(Eφ)
    dotM1,i=−WJ h1,i(E1), dotM2,i=+WJ h2,i(E2),
    dotMφ,i=+WJ hφ,i(Eφ).

Then dot(N1+N2)=dot(N1+Nφ)=0 and dot(U1+U2+Uφ)=WJ(−E1+E2+Eφ)=0, while the full entropy derivative is exactly

    dotS = WJ[−η1(E1)+η2(E2)+ηφ(Eφ)]
         = W(F−R) ln(F/R) ≥ 0.

There is no post-step clipping, fixed external reservoir or BGK relaxation parameter. At common T and μ1=μ2+μφ, η(E) is affine and detailed balance holds on every event. The integral event weights in the proof pilot are **chosen synthetic positive quadrature weights**, not a physical Γ, actual Yukawa decay width or early-universe coupling solution.

## Independently computed numerical tests

Input: illustrative mφ=86 eV = 8.6×10⁻⁸ GeV, gφ=1, T=1 GeV, Emax=12.1 GeV; 36-point Gauss physical-p quadrature per hat interval; independent adaptive SciPy p²dp integration (epsabs 2e−13, epsrel 2e−11). No cosmological thermal history or observations.

| Uniform E intervals | Nφ from physical hat moments | Independent direct p²dp Nφ | Relative difference | Old erroneous rest-node Nφ |
|---:|---:|---:|---:|---:|
| 16 | 0.12174521280791657 | 0.12174521280791706 | −4.0×10⁻¹⁵ | 10616.026218 (×87199) |
| 64 | 0.12174521280791631 | 0.12174521280791706 | −6.1×10⁻¹⁵ | 165.995785 (×1363.47) |
| 256 | 0.12174521280791605 | 0.12174521280791706 | −8.3×10⁻¹⁵ | 2.713570 (×22.289) |

Same synthetic thermal energy Uφ=0.3283420088705699 and entropy Sφ≈0.43795566463 on all these grids. Independent direct p-integral confirms energy at tested precision. The massless BE thermal N up to 50T is 0.12179382823357307, versus exact ζ(3)/π² = 0.12179382823357308. The large point occupation at p=0 is **not** a finite condensate density.

For moving on-shell χ1=.1, χ2=.02, φ=86e−9 and T=.8, μ1=.005, μ2=.012, μφ=−.007, the chemical fixed point has |Rφ|≈1.47×10⁻²⁰ under 16 synthetic positive weighted events (roundoff). A χ1 inverse-fugacity perturbation +0.3 sin(3E) produces negative φ net rate Rφ=−1.16250902×10⁻⁵, **positive** microscopic entropy 1.91204548×10⁻⁶, dark-charge defect ≲2×10⁻²¹ and total energy defect 3.18×10⁻²². Entropy evaluated via independent event affinities and moment dual gradients agrees within ≈10⁻²¹. This is a structural counterexample to any assumption that φ is always produced, not a physical 86-eV production rate.

### Negative control: correcting only the first-cell occupation is unsafe

At N16 the physically integrated first half-cell Bose average is bbar≈3.49053664, while the real rest-point b(E=mφ)≈1.16×10⁷. Assigning this average to the old **rest-energy node** changes ηφ,0 from ≈8.6×10⁻⁸ to ln(1+1/bbar)≈0.251, even though E=mφ still physically requires ≈8.6×10⁻⁸ at T=1. For a synthetic common-T four-event Yukawa sample with arbitrary unit event weights this causes **spurious net φ production +0.739305655 and positive entropy 0.168685498**, whereas the consistent moment basis gives |Rφ|≈4.66×10⁻¹⁵ and negligible entropy. Fixing a scalar total number *alone* can therefore destroy exact quantum equilibrium: the moment reconstruction and source must be repaired together.

## Reproduction, literature and limits

Source: [phi_bose_moment_galerkin.py](phi_bose_moment_galerkin.py); [8 focused regression tests](../tests/test_phi_bose_moment_galerkin.py); [primary JSON](results/phi-bose-moment-galerkin-2026-11-22.json); [archive stage](../archive/stages/2026/11/2026-11-22-phi-bose-moment-galerkin.md).

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_moment_galerkin.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_moment_galerkin

Eight focused tests passed locally including independent physical number/energy integrals, massless ζ(3), Jacobian finite-difference/symmetric negative controls, two conservation/H cases and a wrong-cell-average negative control. The preceding false-equilibrium Scientific CI #45046227473 passed on a6c2e975; **new CI has not been claimed until checked after this commit**.

Methodological comparison: [Fast conservative and entropic numerical methods for the Boson Boltzmann equation](https://arxiv.org/abs/1009.2748), [A numerical scheme for the quantum Boltzmann equation with stiff collision terms](https://numdam.org/articles/10.1051/m2an/2011051/). Entropy-dual kinetic ideas preexist; these papers are context, not an independent reproduction of this particular Yukawa calculation. No external code was copied.

**Scientific boundary:** The new representation is NOT a full nonlinear four-species time scheme, complete SM/χ1χ2/φ finite-T microphysics, actual 86-eV FRW production, physical Γelastic/H, a proven Bose condensate, ΩDM or CLASS-ready. An actual BEC would need separately derived condensate δ³(p) dynamics and a controlled μφ→mφ limit. At b=0 the finite η inverse-fugacity representation is a limiting problem that is not solved here. The next high-value gate is an **entropy-stable implicit moment-space timestep**, solving M(ηnew)=Mold+dt collision(ηnew) with positivity/backtracking and exact Yukawa charges/energy; then compare at 3+ dt/momentum resolutions and rejoin the full χ1χ2 elastic+lepton+φ scheme only after verified physical normalization. CLASS-ready=false.
