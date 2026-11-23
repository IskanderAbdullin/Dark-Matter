# Nonoverlapping cellwise Fermi/Bose entropy-dual moments remove a specific false IR tangent violation (2026-11-23)

## Research question and predeclared falsification gate

In the previous [physical Bose tangent audit](PHI-BOSE-MOMENT-TANGENT-AUDIT.md), overlapping piecewise-linear energy hats produced a mathematically rigorous but **numerically induced** obstruction. On a quadratic N16 scalar energy grid with mφ=86e−9 GeV and χ1 colder than χ2,φ (Tχ1=.9,Tχ2=Tφ=1), the first scalar hat collision source was exactly zero, and **every** scalar hat had nonpositive source. Yet its neighboring source coupled via the off-diagonal physical susceptibility to ηφ(rest), giving dηrest/dt≈−2.9982×10⁻⁹. With initial ηrest=8.6×10⁻⁸, this predicted a local bound dt≈28.68; the old implicit solver passed dt20 and rejected dt30. It was **not** evidence of Bose condensation.

A falsifiable numerical alternative: use a discontinuous energy-cell P1 entropy-dual basis with **two local endpoint coordinates per physical cell**, no cross-cell hats. The interpolation remains affine within each cell (Σh=1 and ΣE_i h_i=E), with a duplicated but possibly discontinuous η at interfaces. Test whether (1) first-cell collision absence ⇒ both first-cell moment rates **and both inverse fugacity tangent rates** vanish exactly; (2) actual physical p²dp Bose density and thermal equilibrium survive; (3) the identical boosted on-shell Yukawa 1↔2 source preserves both dark charges, total energy, positive microscopic quantum entropy; (4) the full nonlinear implicit solver succeeds at old failing dt30 without η clipping; and (5) independent energy-grid and time checks honestly record residual nonconvergence.

## Formulas and physical measure

For a species with rest mass m, energy-cell edges E_j<E_(j+1), define **two nonoverlapping local basis functions** on each cell:

    h_(j,L)(E) = (E_(j+1)-E)/(E_(j+1)-E_j),
    h_(j,R)(E) = (E-E_j)/(E_(j+1)-E_j).

They vanish outside their own cell; the two endpoints are **not shared** with neighboring cells. At exact interfaces the selected right or left trace is a discretization convention. For each cell, reconstruct the fermion/boson inverse fugacity η(E)=h_L η_(j,L)+h_R η_(j,R), with f=1/(expη+1) and b=1/(expη−1), subject to ηφ,endpoints>0. The real scalar physical integrated moments are

    M_(j,a) = g/(2π²) ∫_(p_j)^(p_(j+1)) p²dp h_(j,a)(E(p)) b[η(E(p))],

with the **actual massive dispersion** E(p)=√(p²+m²). The same construction applies to both χ fermions, with the Pauli occupancy. The physical 2×2 per-cell susceptibility is

    J_(j,a)(j,b) = − g/(2π²) ∫cell p²dp h_(j,a)h_(j,b) n(1±n),

negative-definite for finite nondegenerate positive-measure cells and interior occupations. **Cross-cell blocks are identically zero**, unlike the previous overlapping continuous P1 Gram matrix. Therefore if C_(j,L)=C_(j,R)=0 and Jcell is nonsingular, dη_(j,L)/dt=dη_(j,R)/dt=0. This eliminates exactly the specific cross-cell artificial ηrest decrease, *without adding a condensate or clipping Bose occupations*.

Per on-shell physical Yukawa event χ1(P)↔χ2(k)+φ(q), E1=E2+Eφ; the full quantum flux is Jevent=f1(1−f2)(1+bφ)−(1−f1)f2 bφ. Deposit one positive-weight event extent W Jevent into the **same local hats used to evaluate** the three η values, with signs (−1,+1,+1). Then, exactly by partition of unity, energy reproduction and the quantum odds identity:

    dotN1+dotN2 = dotN1+dotNφ = 0;
    dotU1+dotU2+dotUφ = 0;
    dotS = Σ eta_i dotM_i = Σevents W(F−R)log(F/R) ≥0.

The fully implicit physical moment step uses M(ηnew)−M(ηold)=dt C(ηnew), a dense event Jacobian but a block-diagonal physical moment susceptibility. Under a valid nonlinear root, concavity yields Snew−Sold≥ηnew·(Mnew−Mold)=dt Hnew≥0. This is still **first order in time** and is only a conditional finite-event numerical test, not a claim of global existence.

## Independent local numerical checks

The self-contained [new module](phi_bose_cell_entropy.py) was executed in a clean sandbox (NumPy/SciPy); [seven focused tests](../tests/test_phi_bose_cell_entropy.py) passed locally in 0.92 s. It does not depend on the previous overlapping-basis code for its source or susceptibility; the old and new results are compared at the *same physical Yukawa tree-level event model*. All inputs are illustrative: mχ1=.1 GeV, mχ2=.02 GeV, **mφ=86e−9 GeV**, y=.001; χ1 thermal FD(.9), χ2 FD(1), φ BE(1); parent pmax=2 GeV, 12-point p Gauss ×12-point rest-frame decay cosine, physical parent rate g1 p²dp Γrest m1/E1 dc/(4π²), g1=4 (historic numeric convention). Fermion edges uniform E, scalar quadratic edges E=mφ+(Emax−mφ)(j/N)²; each cell has two independent inverse-fugacity endpoints. Each state uses positive physical p-space radial moment integrations, 24-order Gauss per cell.

The independent on-shell scalar source is **net destruction** Rφ=−1.16493012862968×10⁻¹¹ (illustrative natural rate units) with true microscopic H=+1.64846911073637×10⁻¹², matching the earlier moving Yukawa weak reference. The new source has exact paired dark charge and energy defects <~10⁻²². At N16, the **entire first scalar cell** has C_L=C_R=0, and its independently inverted susceptibility gives dη_L/dt=dη_R/dt=0; the previous *overlapping* N16 quadratic had Cfirst=0 but dηrest/dt=−2.9982×10⁻⁹.

The new **nonlinear implicit** solver, with no η clipping and exact physical integrated Bose moments, gives:

| φ energy cells N | dt (illustrative GeV⁻¹) | Newton iterations | ΔNφ | ΔS physical | minimum ηφ |
|---:|---:|---:|---:|---:|---:|
| 8 | 30 | 1 | −3.494744×10⁻¹⁰ | +4.945389×10⁻¹¹ | 1.90041×10⁻⁷ |
| **16** | **30** | **1** | **−3.494750×10⁻¹⁰** | **+4.945377×10⁻¹¹** | **8.6×10⁻⁸** |
| 16 | 1000 | 2 | −1.164681×10⁻⁸ | +1.648143×10⁻⁹ | 8.6×10⁻⁸ |
| 32 | 1000 | 2 | −1.164642×10⁻⁸ | +1.648078×10⁻⁹ | 8.6×10⁻⁸ |
| 64 | 1000 | 2 | −1.164616×10⁻⁸ | +1.648032×10⁻⁹ | 8.6×10⁻⁸ |

The mathematical **old N16, dt30 failure is removed for this prescribed finite event graph**, and dt1000 also succeeds here; all net φ production remains negative, total entropy grows, chemical equilibrium with independent nonzero μ1=μ2+μφ is a fixed point, and Fermi/Bose bounds hold. The 3-grid direct Bose number integration agrees with an independent radial Bose integral to relative <5×10⁻¹², with no fake rest-cell density or invented δ(p) mode. Seven tests include exact heat/charges/H on three energy meshes, analytic affine hat and physical susceptibility block structure, the old failing dt, and massless/massive physical radial integrals.

## **Critical negative:** no numerical or physical relic convergence yet

At the much longer **still synthetic** dt=10⁷ GeV⁻¹, the net scalar changes at N=8/16/32/64 are **−7.59298/−7.28080/−6.85008/−6.42321 ×10⁻⁵**. The change is not converged under N-refinement (N32→64 differs ~6.2%, N8→64 ~18.2%), even though every run conserves dark charges/energy, has positive entropy and remains Bose admissible. A larger parent pmax was previously known to change the physical inverse Yukawa weak source strongly: parent cutoff2 misses ~65.6% of |Rφ| relative to cutoff16 for an IR-curved spectrum. Here pmax2, angular 12² and the discontinuity at energy-cell interfaces are **not independently converged**; the apparent cure is structural **only**, not an accurate physical spectrum.

Also, discontinuous P1 has twice as many inverse-fugacity degrees of freedom per energy cell, and discontinuities can create new shape/threshold alias. There is no guarantee the full true quantum collision has an admissible finite-η solution for every dt, no rigorous global time error or zero-mode condensate equation, and no coupled χ1χ2–SM/e±γν/Z′ plasma or FRW H. A true BEC cannot be inferred from a near-zero Bose η or from any numerical time-step failure.

## Reproduction and scientific context

- [Executable code](phi_bose_cell_entropy.py)
- [Seven focused regressions](../tests/test_phi_bose_cell_entropy.py)
- [Raw N/time table](results/phi-bose-cell-entropy-2026-11-23.json)
- [Permanent stage record](../archive/stages/2026/11/2026-11-23-phi-bose-cell-entropy.md)

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_cell_entropy.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_cell_entropy

Related independent numerical method (not a validation of this specific Yukawa code): [entropy-stable discontinuous Galerkin Boltzmann moments, arXiv:1602.01312](https://arxiv.org/abs/1602.01312); [deterministic quantum Boltzmann–Nordheim and Bose condensation numerical analysis, JCP (2023)](https://www.sciencedirect.com/science/article/pii/S0021999123002929). No third-party code copied.

**Next action:** verify real Scientific CI of this new stage, then independently converge 3+ scalar N/nonuniform IR grids with actual parent pmax=4/8/12/16 and angular 12/24/48, compare DG-P1 and earlier continuous-hat physical η integral at equal mode costs and thermal/nonthermal states, and audit discontinuity effects near energy interfaces. Only thereafter couple the entropy-compatible φ mode equation to shared-lepton χ1χ2 Pauli elastic, e±γν/Z′ chemistry and FRW. **CLASS-ready=false; physical Γelastic, ΩDM, BEC, frozen fφ(q) are unproven.**
