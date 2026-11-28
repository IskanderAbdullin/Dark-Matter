# A zero energy transfer is not a zero dark-temperature drift

Extend the [separate FD number and heat nullclines](INDEPENDENT-NONLINEAR-FD-PORTAL-NULLCLINES-2026-11-27.md)
without adding any interaction or model species. Use the same ideal
symmetric Dirac population and thermal scalar closure but hold each
instantaneous state fixed for the collision evaluation. F is the
net annihilation event rate and Q is total energy transfer from
dark to visible (timelike plus direct charged-lepton elastic).
At fixed volume, n_D′=−2F and ρ_D′=−Q. From the exact FD energy
derivatives at fixed T_D and n_D,

    C_V|n T_D′ = −Q + 2χF,   χ = (∂ρ_D/∂n_D)|T_D > 0.

Here χ is the Fermi occupation-derivative mean single-particle
energy, not μ_D. When particles disappear at zero net energy
transfer, the remaining system can warm. When T_D′ vanishes with
F nonzero, energy still moves between the baths. No simultaneous
zero of both bath temperatures and pair evolution is possible away
from detailed balance in this restricted model.

Define G=Q−2χF. At common temperature and η=0, both F and Q vanish,
so derivatives of χ do not enter to first order. With the event
energy moments M0, M1, M2 and elastic conductance K from the previous
gate, the local temperature-null slope is

    η_T′(1) = −(M2/T_SM + K − 2χ M1/T_SM)/(2(M1 − 2χ M0)).

At finite displacements χ itself varies, and G is not a positive-
weighted rate. In particular no monotone-root theorem from the
pure timelike heat and event moments is available for G. Count
sign-change roots only in the stated scan interval; do not interpret
the fold as a physical phase transition, an evolving dark-sector
abundance, or multiple equilibrium states.

## Checks

Check temperature drift independently by perturbing the conserved
dark number and energy ledgers, then inverting the FD state to
recover T_D. For the reported folds combine the sign of the
local peak across the fold, opposite signs on both eta sides of
the peak, negative eta curvature, a positive ratio derivative
and a separate two-dimensional solution of G=∂G/∂η=0. Repeat
with shifted quadrature orders, thermal cutoffs and FD state
integration order. No expansion terms or scalar production are
included in this check.

## What the zeros do in a short time evolution

At T_SM=.050 GeV and T_D/T_SM=.80, the two temperature roots lie
at η=2.20966 and 3.74720; the combined energy zero at the same
temperatures is near 1.25331. At both temperature roots the dark
pair number is decreasing and Q>0. Under a fixed visible
thermostat, T_D′=0 at a root but η′=−2F/(∂n_D/∂η)<0.
Consequently T_D″=−G_η η′/C_V|n at that instant. The lower
and upper roots give a local temperature minimum and maximum,
respectively. Independent short integrations in conservative
number and energy variables, including the thermostat entropy
ledger, recover the signs and local curvatures. These are not
stationary points of the complete number-energy evolution.

## Two examples of a fold in the temperature-zero set

| Visible temperature (GeV) | Critical T_D/T_SM, with elastic | η at the fold, with elastic | Critical T_D/T_SM, timelike only |
|---:|---:|---:|---:|
| .015 | .85179139 | 2.5458495 | .71943782 |
| .100 | .80850028 | 2.8250974 | .64841552 |

Near each reported fold, raising the ratio through its critical value
creates two nearby sign-change roots in η. A ratio step of ±.0005
tests this directly. The inferred critical ratio changes by less
than 1.2e−7 between quadrature orders. At fixed ratio .80 no roots
were detected in [−5,5] at .015 and .100 GeV with elastic heat;
the local solutions lie on the two-root side already under timelike
exchange alone. The full scan is numerical rather than an enclosure
excluding all other roots. The thresholds depend on the adopted
ideal thermal states, the vacuum portal and included species.

Replay from the kinetics directory with no external equation-of-state
data or imports outside this pull request:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error fd_portal_dark_temperature_null_reference.py
python -W error fd_portal_thermal_turn_reference.py
```
Receipts: [root and fold checks](results/dark-temperature-nulls-2026-11-27.json.gz)
and [conservative time integration](results/fd-portal-thermal-turns-2026-11-27.json).
