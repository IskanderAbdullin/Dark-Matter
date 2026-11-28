# Number and heat zeros separate in the local full-FD two-bath portal

Consider the same symmetric Dirac pair and charged leptons as the
[nonlinear two-bath reference](INDEPENDENT-NONLINEAR-FD-PORTAL-REFERENCE-2026-11-26.md).
Hold the two temperatures fixed, then vary the pair chemical
potential η=μ_D/T_D. No equation of state, expansion clock or
numerical trajectory is required. Define F as the net event rate
for annihilation (positive dark pair loss) and Q as dark-to-visible
energy transfer. The net dark fermion-number source is −2F.
An event of lab energy E has affinity A(E)=2η+βE,
where β=1/T_SM−1/T_D. Its flux J(E) has the sign of that affinity
by exact forward/reverse FD detailed balance. Then

    F = ∫ J(E),   Q = ∫ E J(E),   σ = ∫ A(E) J(E) = 2ηF + βQ ≥ 0.

These statements use full Pauli factors and a positive kernel,
not the Maxwell–Boltzmann cancellation or the mediator pole energy
as the deposited heat. Timelike F and Q are strictly increasing in
η: increasing dark fugacity increases the forward occupations and
decreases the reverse blocking factors at each event. If F=0 and
β≠0, use E*=−2η/β, so Q=∫(E−E*)J(E) has the sign of β strictly
for the retained continuum. Consequently the timelike heat zero
lies below the number zero when the dark bath is hotter, and above
it when the dark bath is cooler. Direct elastic heat adds a term of
the sign of β and cannot change F, so it moves any heat zero farther
from the number zero; the monotonic proof is for the timelike
moments only. The full moments cannot vanish simultaneously at
distinct temperatures because their collision entropy is strictly
positive there.

## Independent common-temperature limit

With r=T_D/T_SM and a common reference point r=1, η=0, integrate the
equilibrium forward/backward event weight to get M0, M1, M2, its
zeroth, first and second lab-energy moments. Then

    η_N′(1) = −M1/(2 T_SM M0)
    η_Q′(1) = −M2/(2 T_SM M1).

Their separation is −Var(E)/(2 T_SM ⟨E⟩). The leading elastic heat
conductance K changes only the second expression, to
−(M2/T_SM+K)/(2 M1). This limit is checked against a finite
difference of separate nonlinear roots at r=1±.0003, with no
change to the timelike or elastic kernels.

## Finite displacement tests

| T_SM (GeV) | T_D/T_SM | η at F=0 | η at timelike Q=0 | η at combined Q=0 |
|---:|---:|---:|---:|---:|
| .050 | .80 | 1.006244 | 1.099864 | 1.253308 |
| .050 | 1.20 | −.783773 | −.862620 | −.986247 |
| .300 | .80 | .985870 | 1.099657 | 1.258220 |
| .300 | 1.20 | −.808911 | −.939538 | −1.051229 |
| .6882487043 | .80 | 1.968986 | 1.986381 | 1.986936 |
| .6882487043 | 1.20 | −1.346761 | −1.357548 | −1.357682 |

All twenty unequal-temperature controls resolve separate heat and
number zeros. At .30 GeV and r=.80 the heat-only elastic channel
more than doubles the separation, from .113787 to .272350 in η.
At .68825 GeV and the same ratio the corresponding separation
moves only from .017395 to .017950. Thus a single chemical
threshold does not determine the heat sign, and how far apart
the thresholds lie depends on the portal state. The blocked
timelike exchange alone can even carry heat against the temperature
gradient when chemical disequilibrium supplies sufficient entropy
production; the elastic term by itself never does. These are local
charged-lepton vacuum controls. They establish neither an elastic
or chemical decoupling temperature nor any physical abundance.

Reproduce in the cloud from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_nonlinear_nullclines_reference.py
```
The [receipt](results/fd-portal-nonlinear-nullclines-2026-11-27.json.gz)
includes rate and entropy ledgers at each root, the two quadrature
orders, all five common-temperature derivative checks, and the
independent adaptive timelike quadrature at four selected controls.
