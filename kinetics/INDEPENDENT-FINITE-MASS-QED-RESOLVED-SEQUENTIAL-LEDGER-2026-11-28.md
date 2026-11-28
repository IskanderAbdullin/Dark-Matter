# Finite-time storage in a separately resolved charged pole

The [occupied roots and their spin residues](INDEPENDENT-FINITE-MASS-QED-CUT-LOCAL-SEQUENTIAL-SUBTRACTION-2026-11-27.md) identify a charged intermediate state on each of two distinct bare poles. This is a **separately resolved diagnostic** for each pole, with a phenomenological choice of step clocks. It provides a check on kinetics ledgers when an intermediate occupation is evolved, not an additive scattering channel or a prescription for regulating the coherent overlap.

## Exact two-step ledger

Let f be the occupation of a single internal fermion mode, and let R₁ = a(1−f)−bf and R₂ = cf−d(1−f). All four controlled coefficients are positive. Then

```text
df/dt = R1 − R2
H = a+b+c+d,  f* = (a+d)/H,  J* = (ac−bd)/H
f(t) = f* + (f0−f*) exp(−Ht)
integral R1 dt = J*t − (a+b)(f0−f*)(1−exp(−Ht))/H
integral R2 dt = J*t + (c+d)(f0−f*)(1−exp(−Ht))/H
```

Thus identifying the two currents at finite time discards the stored change of f. The nonnegative local entropy production is the sum of each directed rate difference times the log ratio of that edge's forward and inverse rates, including f and 1−f. Integrating gives the two externally weighted currents plus the change of binary fermion entropy. At stationarity the result reduces to J* log(ac/bd) ≥ 0. Thirty-six toy time controls compare against independent DOP853 integration and entropy quadrature.

For leg kinematics K carries positive energy ω, supplied by the free charged-bath lesser cut L; the inverse deposits that energy and uses the greater cut G. Scalar emission uses (1+n) and absorption uses n. At the parent intermediate pole, the charge-bath step precedes the scalar step; at the daughter pole the order is reversed. The fitted clock on each step normalizes its forward-plus-inverse coefficient, while retaining the ratio of the two microscopic statistical factors. Both factor orderings give the same form of external cycle affinity:

```text
logit(f_parent) − logit(f_daughter) + log(L/G) + log(1+1/n).
```

A common temperature, common conserved fermionic chemical potential, and neutral scalar with zero chemical potential give zero affinity since E_parent + ω = E_daughter + E_scalar. Each pair of the step rates, independently, then selects the same thermal occupation at the intermediate energy. Nonthermal clocks are not proportional to a computed width: their sums were chosen in the ratio .01, 1, or 100 and the comparison is over one resulting relaxation time.

## Archive replay

At all 48 independently occupied spacelike roots across the fine and third archives, the source takes the local charged and scalar occupations and the [exact free transverse cut pair](INDEPENDENT-FINITE-MASS-QED-TRANSVERSE-STATISTICS-2026-11-27.md). The diagnostic restores the occupied level's storage explicitly. The energy ledger includes the changed occupation times its intermediate energy, the lost parent, the produced daughter and scalar, and the energy taken from the bath. Maximum independent defects are 5.68e−17 in energy, 2.78e−17 in stored number, 1.67e−16 in integrated entropy, and 1.12e−16 in the common-thermal target. The maximum finite-time difference between integrated scalar output and its stationary prediction in the normalized tests is .103 of a mode occupation unit. It has no physical unit of scalar production rate.

Run `python -W error finite_mass_qed_resolved_sequential_ledger_reference.py` from `kinetics/`. See the [receipt](results/finite-mass-qed-resolved-sequential-ledger-2026-11-28.json) for both cumulative currents and the occupation stored at each test point.

The separately held bath factors do not reconstruct a finite-width intermediate propagator. In a self-consistent evolution, the cut already represents charged matter; adding its resolved sequential contribution to a collision term which already evolves those reactions without matching would double count. In particular, the two orderings must not be summed incoherently at the [coalescent bare poles](INDEPENDENT-FINITE-MASS-QED-CUT-POLE-COALESCENCE-2026-11-27.md). No full virtual term, interacting population, scalar NLO rate, or cosmological claim is inferred.
