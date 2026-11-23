# Exact physical scalar shell capacity and midpoint no-go (2026-11-23 UTC)

**Scope:** independent analytic phase-space capacity component for a real massive scalar, not a complete 1↔2 strong Yukawa collision kernel. Units: natural ħ=c=1, momentum and energy in GeV; m_phi=86 eV = 8.6e-8 GeV is a historical benchmark only.

For one scalar internal state, physical occupation f(p) is normalized by n=∫d³p/(2π)³ f=∫_0^∞ p²dp/(2π²) f. For on-shell E=√(p²+m²), p²dp=E√(E²−m²)dE. Therefore the **exact positive capacity** of an energy cell [L,U] with U≥L≥m is

V_phi([L,U]) = [(U²−m²)^(3/2) − (L²−m²)^(3/2)]/(6π²).

Multiply by an explicitly justified internal degeneracy g_phi (g=1 for a single real scalar). For a cell crossing threshold replace L by max(L,m); if U≤m its physical capacity is zero. In a numerically narrow cell evaluate the difference using a stable factored form or compensated arithmetic. This capacity is independent of the collision matrix element and **must not** be replaced by an arbitrary event-graph node weight. Energy capacity for constant occupation is independently ∫_L^U E²√(E²−m²)dE/(2π²), not generally E_node V_phi; this distinction obstructs exact energy conservation from a single constant-energy bin representation.

## Reproducible independent numerical counterexample

Command (Python 3 stdlib):

```python
import math
m=86e-9
for U in (1e-6,1e-4,0.00860360):
    V=((U*U-m*m)**1.5)/(6*math.pi**2)
    E=(m+U)/2
    Vmid=(U-m)*E*math.sqrt(E*E-m*m)/(2*math.pi**2)
    print(U, V, Vmid, Vmid/V-1)
```

| Upper edge U (GeV) | Exact capacity (GeV³) | Midpoint capacity (GeV³) | Relative midpoint bias |
|---:|---:|---:|---:|
| 1e-6 | 1.6699867896979555e-20 | 1.3480304644408889e-20 | -19.2789743753% |
| 1e-4 | 1.688684520610624e-14 | 1.2676011889145414e-14 | -24.9355830859% |
| 0.00860360 | 1.0754485462607377e-8 | 8.065944720637698e-9 | -24.9992503251% |

As m/U→0 and L=m→0, midpoint rule evaluates E² at U/2, giving U³/(8π²), whereas exact integral is U³/(6π²): **−25%** relative error. This is an analytic negative control for coarse scalar energy cells, not an estimate of the error in an actual Yukawa reaction rate; Bose/Pauli factors and kinematic support still require their own integration. In particular, the previous first-cell cut-cell weak source cannot be corrected by simply multiplying by 4/3: its rate integrand is different.

## Boundary for the next stage

1. Derive parent and fermionic daughter capacities with correct degeneracies and physical phase-space normalization.
2. Match the **same invariant on-shell event measure** to all three species' weak N/Q moments; compare against independently converged continuous reference with angular occupations integrated inside each cut cell.
3. Check finite-volume reversible event detailed balance and entropy only after proving physical energy handling; arbitrary bin-centre E1=E2+Ephi must not be assumed.

**Negative result:** exact cell capacity by itself does not yield a positive conservative reversible nonlinear collision operator; no strong time evolution, BEC, Ω or CLASS claim follows.
