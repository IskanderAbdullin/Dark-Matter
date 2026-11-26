# Implicit bounded quantum two-bin/network with finite SM heat — a REDUCED reference

**Purpose:** independently construct a **real time-evolving, Pauli-bound preserving, particle/energy/entropy-conservative reduced fermion scattering network with a finite radiation bath** as a numerical algorithm benchmark before implementing the physical relativistic 5D collision transition matrix. [Implementation](fermion_finite_bath_pairs.py) / [tests](../tests/test_fermion_finite_bath_pairs.py). **This is not a microscopic χe contact kernel. The pair rates below are chosen positive constants without scattering amplitude calibration.** No Bose production, cosmological Γ(T), 86-eV relic, FRW expansion or CLASS.

## Exact discrete finite bath reaction ledger

For ordered states E_i<E_j and corresponding fermion capacities V_i,V_j (the maximum particle numbers in each energy cell), (N_i=V_i f_i), with 0≤f_i≤1, choose a phenomenological Markov transition coefficient Γ≥0. The physically motivated detailed-balance *ratio* between upward and downward rates is exp(−ΔE/T). One convenient parametrization is

```text
J_up   = Γ exp(−ΔE/T_SM) f_i (1−f_j)
J_down = Γ f_j (1−f_i)
J      = J_up−J_down
ΔE=E_j−E_i>0.
dN_i/dt=−J,  dN_j/dt=+J, dE_dark/dt=ΔE J, dE_SM/dt=−ΔE J.
```

Both χ-number and combined dark/SM energy are conserved at the **event** level for *any positive* Γ. A true physical χe pairwise Γ must later be obtained from the invariant 5D collision operator, including electron recoil, Pauli electron blocking, phase-space weights, and the actual lepton bath. **Nothing here establishes such a microscopic Γ.**

Let (E_{SM}=A T_{SM}^4), (S_{SM}=(4/3)A T_{SM}^3), A>0 constant in the fixed-volume test (equivalent to ideal radiation EOS, no FRW work). The total entropy

```text
S_total = Σ_i V_i [−f_i ln f_i −(1−f_i) ln(1−f_i)] + (4/3) A T_SM^3.
dS_total/dt = J ln(J_up/J_down) >=0.
```

The identity uses `ln(J_up/J_down)=ln[f_i/(1−f_i)]−ln[f_j/(1−f_j)]−ΔE/T_SM` and follows from (x−y)ln(x/y)≥0; endpoints understood by continuity.

## Fully implicit bounded pair update

For any dt≥0 solve the **single scalar nonlinear** equation

```text
ΔN = dt * [Γ exp(−ΔE/T_new) f_i,new(1−f_j,new)
           −Γ f_j,new(1−f_i,new)]
f_i,new= f_i,old−ΔN/V_i
f_j,new= f_j,old+ΔN/V_j
E_SM,new=E_SM,old−ΔE ΔN
T_new=(E_SM,new/A)^(1/4).
```

The physically allowed bracket is `ΔN∈[−min(Vj fj,Vi(1−fi)), min(Vi fi,Vj(1−fj),ESM/ΔE)]`. On this domain the scalar residual has opposite boundary signs (unless a stationary endpoint) and the root exists. Brent finds the unique physical branch in the tested regimes. This **does not clip occupations** after updating: Pauli 0≤f≤1 and radiation ESM≥0 follow from the domain. Total number and energy are conserved directly in the definitions. The entropy increase is numerically checked for short and extremely stiff dt values; an implicit finite step inherits discrete entropy monotonicity because its post-step rate still has the detailed-balance sign and concavity of physical total entropy along the reaction coordinate.

A multi-level network is built by Lie-splitting **actual sequential implicit pair updates**, each of which preserves total number, total energy and nondecreasing entropy. The network is not globally time-exact: **time-order/splitting errors remain** even when all conservative residuals vanish.

## Independent tests and measured scope

Tests include 130 deterministic random positive masses, capacities, temperatures, rates and dt spanning 6 orders of magnitude; both directions of heat; exact FD(TSM,µ) fixed points for dt=10⁻4..10⁵; fully occupied/empty Pauli endpoints without clipping; 3-bin finite bath network split; standalone particle/energy residuals and entropy; fail closed on invalid inputs. For illustrative pair E=(.1,1), V=(3,5), f=(.9,.1), ESM=50, A=10, Γ=2, a very stiff dt=1e5 gives approximately f=(.491293,.345224), ESM=48.896492, TSM=1.487029 and entropy gain +1.960507, with |energy residual|≈0 in double precision. Those are **arbitrary test units**, not GeV-scaled dark-matter rates.

At 3-bin E=(.1,.6,1.3), V=(1,3,5), f=(.9,.1,.1), ESM=50, A=10, edges (0→1,Γ=2),(1→2,Γ=1), total dt=1, the 10-vs-100-substep occupations differ by >0.001 even though both conserve energy and number within 1e−12 and increase entropy. This regression illustrates that green conservation tests do **not** guarantee temporal convergence.

**Scientific limitation:** The conserved-bin capacities, phenomenological Γ, absence of actual e±/ν energies, elastic recoil coupling, expansion and 1↔2/Bose physics mean this is a deliberately reduced quantum thermalisation **numerical design reference**, not a solved physical χ+SM Boltzmann equation. Applying it as if the real `γ_el` were derived from the Yukawa or contact Lagrangian would be scientifically wrong.

## Next

Implement a transition kernel built from the **physical 5D scattering matrix element and final-state electron Pauli** and transfer its microscopic event rate into conservative q-bin capacities, reproducing continuum weak moments and explicitly tracking SM heat and entropy. Only after actual collision rates/moments, angular/q/time convergence, gauge/UV/SM completion, and coupling to 1↔2 BEST could one model late freezeout.
