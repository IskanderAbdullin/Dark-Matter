# Exact finite-parent-cutoff kinematic IR floor: refine the previous Bose-DG mesh interpretation (2026-11-23)

## New falsifiable question

The previous [DG-P1 within-cell tangent negative](PHI-BOSE-DG-TANGENT-COUNTEREXAMPLE-AUDIT.md) found that at fixed parent pmax=2 GeV, scalar quadratic energy meshes N24/32/40 have a negative first ηdot, yet N48/64 first-cell source and tangent vanish. The interpretation referred to “finite event quadrature misses the narrower first cell”. This is insufficient: **does any exactly on-shell event exist within those cells at the hard parent cutoff, even at infinite radial/angular resolution?** This can be settled exactly without expensive simulations.

## Exact massive boosted daughter minimum

With parent mass m1, fermion daughter m2, strictly positive scalar mass mφ, parent momentum p and two-body rest-frame daughter momentum p* and energy Eφ*, the Lorentz-boosted boson energy is

`Eφ(p,c)=[sqrt(p²+m1²) Eφ*−p p* c]/m1,` with `−1≤c≤1`.

The minimum at fixed p is Eφ(p,+1). **Numerically stable rationalization**:

`Eφ,min(p) = [p² mφ² + m1² (Eφ*)²]/[m1 (sqrt(p²+m1²) Eφ* + p p*)].`

The derivative is `dEφ,min/dp=[p Eφ*/sqrt(p²+m1²)−p*]/m1`; hence the minimum decreases with parent p until `pturn=m1 p*/mφ`, where its value is exactly mφ, and increases afterward. Thus for a finite parent cutoff `pmax < pturn`, the **global** daughter infrared support floor is `Eφ,min(pmax)`, independent of daughter-angle quadrature or collision signs.

At the actual benchmark m1=.1, m2=.02, mφ=86×10⁻⁹ GeV, rest p*≈.048 GeV and `pturn≈55,813.95 GeV`. The earlier parent pmax2/8/16 GeV regimes are all on the decreasing branch. Their strict minimum daughter energies are:

| Parent cutoff pmax | Exact allowed Eφ,min | Equivalent eV | Quadratic first-cell exclusion threshold N |
|---:|---:|---:|---:|
| 2 GeV | 0.00119925094 GeV | 1.199×10⁶ eV | **N≥43** |
| 8 GeV | 0.00029998829 GeV | 3.000×10⁵ eV | N≥166 |
| 16 GeV | 0.00014999855 GeV | 1.500×10⁵ eV | N≥329 |

The grid used by [original DG physical cell module](phi_bose_cell_entropy.py) has first scalar energy edge `Eedge=mφ+[sqrt(pmax²+m1²)+.2−mφ]/N²`. At pmax2 its N32 first upper edge is 0.0021509633 GeV **above** the exact Eφ floor and therefore potentially source-connected; N42 remains above at 0.00124866783, while **N43 is below** at 0.00119126953, as are N48 at 0.00095603147 and N64 at 0.00053780532 GeV. Consequently **no possible physical on-shell event** from any p≤2 GeV can enter first scalar cell N≥43 — even with perfect outgoing-angle and parent-momentum integration. The earlier N48/N64 Cfirst=0 is *exact finite-cutoff kinematic exclusion*, not merely finite quadrature under-sampling. The negative N32 tangent is a real numerical problem in a physically accessible first cell and remains unchanged.

For pmax8 and pmax16, by contrast, N64 first cell is accessible, so the empty-cell argument cannot be carried across cutoff changes. Mathematically the **full** untruncated parent distribution can access Eφ down to mφ only around p≈5.58×10⁴ GeV in this benchmark; the thermal T≈1 GeV parent occupation is fantastically suppressed there. We have **not** quantified the full probability or physical full collision source at finite temperature, nor included additional production/scattering channels or real cosmological evolution.

## Evidence and corrective decision

- [Independent executable](phi_bose_parent_cutoff_kinematic_ir.py) computes rationalized exact endpoint, turning momentum, physical cutoff-support floor and grid threshold. [Two focused tests](../tests/test_phi_bose_parent_cutoff_kinematic_ir.py) check against six direct boosted shells and 501 outgoing angles each, analytic pturn, and N42/43/48/64 and pmax8/16 boundaries; **2/2 local pytest passed** in ~0.06 s.
- [Numerical compact JSON](results/phi-bose-parent-cutoff-kinematic-ir-2026-11-23.json); [dated stage](../archive/stages/2026/11/2026-11-23-phi-bose-parent-cutoff-kinematic-ir.md).
- This is an exact **support** theorem conditional on a *hard parent cutoff and this on-shell 1→2 channel*, not a numerical convergence claim or physical BEC/no-BEC statement. A source-free first cell may still contain physical thermal φ particles; zero collision rate there is not a physical δ³(p) condensate.

**NEXT_ACTION:** normalize any future positive exact-energy event network and IR source comparison **only over cells with actual on-shell support**, and extend parent pmax and angle quadrature together. A scalar first-cell moment may be unexcited by this reaction for physical kinematic reasons, rather than a defective scheme. For nonzero accessible cells, the DG-P1 η tangent obstruction remains and must be repaired with a truly realizable exact-energy quantum kinetic closure. Physical φ masses/Γ, SM/FRW and CLASS remain separate open work.
