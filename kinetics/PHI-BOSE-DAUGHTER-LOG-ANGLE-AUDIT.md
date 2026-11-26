# Boosted light Bose daughter angular quadrature: exact log(Eφ) map and a previously hidden ≈10–17% Yukawa weak-rate bias (2026-11-23)

## Prior limitation and prospective falsifiable test

The moving, on-shell χ1↔χ2+real φ scalar Yukawa weak source was previously checked with isotropic mother-rest-frame Gauss cosine quadrature. Its [early numerical reference](MOVING-YUKAWA-BOSE-PAULI-WEAK-AUDIT.md) reported apparently stable mother *radial* n8/12/16/24/32 at fixed nangle16 and φ rate about −1.16×10⁻¹¹ in an illustrative colder-mother thermal state (Tχ1=.9,Tχ2=Tφ=1). **Radial stability does not imply angular stability:** boosted daughter φ can have Eφ≪T in the forward direction and its 1+bφ factor is approximately T/Eφ. A smooth thermal inverse fugacity ηφ=Eφ/T therefore has a strongly *angular* peaked quantum collision even though its η(E) is affine and trivially resolved on any energy mesh.

The predeclared gate: derive a nonnegative exact mapping from outgoing CM cosine c to log(Eφ), independently verify both the angular zeroth measure ∫dc=2 and the Rayleigh–Jeans surrogate ∫dc/Eφ in closed form for multiple mother momenta, compare 3+ cosine nodes on the same on-shell F/R Pauli/Bose weak source against high-order ordinary Gauss, and retain any nonthermal peak that defeats a low-order mapped rule. Do not treat one mapped result as physical 86-eV relic validation.

## Exact on-shell change of variable

For fixed mother momentum p, physical Lorentz boost of scalar Yukawa two-body decay gives

    Eφ(c) = a − b c,
    a=(E1/m1) Eφ*, b=(p/m1) p* >=0,
    E_min=a−b>0, E_max=a+b, −1≤c≤1.

The physical scalar mass mφ=86e−9 GeV is positive in this benchmark; E_min>0 for all finite boosts. For p>0 choose y=log[Eφ/E_min] with L=log(E_max/E_min), y∈[0,L]. Then

    Eφ=E_min exp(y),    c=1−(Eφ−E_min)/b,
    ∫_(−1)^1 dc G(c) = ∫_0^L (Eφ(y)/b) dy G(c(y)).

A standard Gauss node x∈[−1,1] uses y=L(x+1)/2 and the exact **positive weight** w_x (L/2) Eφ/b. At mother p=0 the angular dependence is absent and uniform cosine Gauss is exact.

Two independent exact controls:

    I0 = ∫dc = 2,
    I_RJ = ∫dc/Eφ(c) = log(E_max/E_min)/b.

The log mapping integrates I_RJ **exactly to floating precision at any Gauss order**, because the Jacobian cancels 1/Eφ. The transformation is an exact change of variables and does not alter physical Γrest, mother p²dp, Pauli/Bose F/R, scalar mass, inverse channel or gain/loss normalization.

## Measured independent analytic and quantum 1↔2 controls

Independent numerical implementation [source](phi_bose_daughter_log_angle.py) and [10 focused tests](../tests/test_phi_bose_daughter_log_angle.py) passed locally (0.66 s). For m1=.1,m2=.02,mφ=86e−9 GeV and *eight* cosine Gauss nodes, the relative uniform-c error in I_RJ is **−14.73% at p=1 GeV**, **−27.52% at p=2 GeV**, **−46.50% at p=8 GeV** and **−52.90% at p=16 GeV**. Mapped I_RJ errors remain ≤3.1×10⁻¹³ in these controls; mapped zeroth ∫dc=2 error at the hardest p=16/n8 is ~9.5×10⁻⁸ relative. No artificial angular rate normalization or postprocessing occurs.

### Full physical tree-level 1↔2 weak source, prescribed thermal Fermi/Bose state

Use y=.001, parent pmax=2/8/16 GeV, mother momentum Gauss nr48, real mφ=86e−9 GeV, g1=4 inherited numeric convention, identical fermion FD(Tχ1=.9,Tχ2=1), scalar BE(Tφ=1), and strict boosted two-body shell. The reported Rφ and H are *conditional* on these chosen distributions and cutoff (not cosmologically fitted abundances or rates).

| Parent pmax (GeV) | Uniform angular n12 Rφ | Log-Eφ n12 Rφ | Log-Eφ n64 Rφ | Uniform n12 error vs mapped n64 |
|---:|---:|---:|---:|---:|
| 2 | −1.164929972×10⁻¹¹ | −1.295453962×10⁻¹¹ | −1.295453955×10⁻¹¹ | **−10.08%** in |Rφ| |
| 8 | −2.654770155×10⁻¹¹ | −3.211995890×10⁻¹¹ | −3.212015333×10⁻¹¹ | **−17.35%** |
| 16 | −2.669725105×10⁻¹¹ | −3.232079803×10⁻¹¹ | −3.232104314×10⁻¹¹ | **−17.40%** |

Thus **the earlier source at the same physical pmax=2 and low angular order undercounted net inverse φ destruction by ~10%**, despite its radial stability and independent same-underresolved-quadrature comparisons. The earlier source is still valuable as an *algorithmic reference at a specified discretization*, not a numerically converged physical weak rate. This is a numerical correction, **not** the discovery or exclusion of a dark matter particle.

At pmax=2, mapped n8/n12/n16/n24/n64 thermal Rφ≈(−1.29545147,−1.29545396,−1.29545396,−1.29545396,−1.29545396)×10⁻¹¹, while ordinary n8/n12/n16/n24/n64 gives (−1.06792745,−1.16492997,−1.21777417,−1.26646138,−1.29516166)×10⁻¹¹. Mapped n12 and n64 agree at ~5×10⁻⁹ relative, on *this analytic thermal state*.

### Negative control: curved infrared Bose shape is not cured by one low-order mapping

For the physically admissible positive ηφ(E)=E+0.40 exp[−((E−.04)/.02)²], with identical fermion f and parent pmax16/nr48, mapped n8/n12/n16/n24/n48/n96 Qφ are

    −1.87121441, −1.78944195, −1.72620400,
    −1.71810098, −1.72986696, −1.72985519   ×10⁻¹¹.

At mapped n8 the conditional rate differs by **~8.17%** from mapped n96; mapped n48 vs n96 differs ~0.00068%, but this is not a global source convergence certificate. A narrower arbitrary scalar η shape, near ηφ=0 or different T may demand additional adaptive elements. Parent cutoff pmax2→16 for this sharp IR shape changes |Rφ| from ≈5.952×10⁻¹² to ≈1.730×10⁻¹¹ — **about 65.6% of the pmax16 rate is missed at pmax2** even with properly resolved outgoing angle. This confirms previous UV-tail negative evidence, while exposing its separate angular bias.

## Scope, reproducibility, next gate

- [Independent full-angle source](phi_bose_daughter_log_angle.py)
- [10 local regressions](../tests/test_phi_bose_daughter_log_angle.py)
- [Primary compact angle/cutoff numerical data](results/phi-bose-daughter-log-angle-2026-11-23.json)
- [Dated permanent stage](../archive/stages/2026/11/2026-11-23-phi-bose-daughter-log-angle.md)

    OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phi_bose_daughter_log_angle.py
    OPENBLAS_NUM_THREADS=1 python -m kinetics.phi_bose_daughter_log_angle

The derivation is original for this on-shell test. Related independent numerical context: [entropy-stable DG Boltzmann moment methods](https://arxiv.org/abs/1602.01312) and [Boltzmann–Nordheim quantum numerical analysis](https://www.sciencedirect.com/science/article/pii/S0021999123002929). These papers do not validate this specific 86-eV cosmology.

**Next critical calculation:** feed these mapped, physically weighted events to both the old overlapping and the new discontinuous P1 physical-p²dp nonlinear inverse-moment source with *identical* parent cutoff and 3+ p/c/IR energy/time grids. Resolve the interface-discontinuity error and the previously measured long-time Nφ grid drift. Then independently couple χ1+χ2 lepton vector Pauli, ν/e±γ/Zprime bath and FRW. No validated Γelastic, Bose zero-mode/BEC, ΩDM, frozen fφ(q) or CLASS-ready distribution.
