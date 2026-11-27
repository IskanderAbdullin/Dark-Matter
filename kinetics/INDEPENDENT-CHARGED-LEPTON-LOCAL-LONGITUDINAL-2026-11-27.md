# Local longitudinal control with separate heat and momentum exchange

Use the [common-state FD portal coefficients](INDEPENDENT-CHARGED-LEPTON-COMBINED-FD-MOMENTUM-2026-11-27.md)
in a static thermal two-fluid moment test. The dark fermions and scalar
share a temperature and bulk velocity by closure. The visible bath
is a radiation fluid, with a *chosen* enthalpy ratio h_V/h_D of 1 or
10. There is no radiation anisotropic stress, gravity, expansion,
continual production or physical visible equation of state.

The perturbations are total dark fermion and antifermion density, dark and visible energy densities,
and the two longitudinal velocities. At common T and η_D=0, use
enthalpy h=ρ+P, heat capacity C_V|n, and χ=(∂ρ_D/∂n_D)|T. The
dark perturbations obey

    δT_D = (δρ_D − χ δn_D)/(C_V|n)
    δη_D = δn_D/(∂n_D/∂η_D)|T − χ δT_D/T²
    δP_D = (h_D/T) δT_D + n_D T δη_D.

Let F be net annihilation events, Q dark-to-visible energy flow, and
Δu=u_D−u_V. The number, energy and momentum sources on the dark side
are −2F, −Q, and −K_p Δu. The visible energy and momentum gains
have exactly the opposite signs; there is no visible-number variable.
The pair/heat coefficients are the equilibrium energy moments and
elastic energy conductance, not a single damping parameter.

## Independent structural checks

The streaming generator is skew in the thermodynamic quadratic
metric, and exchange decreases that norm. For the dark thermodynamic
part the metric is δn_D²/n_η + (C_V|n) δT_D²/T²; add the analogous
visible energy term and h_i u_i²/T for each velocity. The loss matrix
is independently assembled from the positive number/heat moment
matrix in (2δη_D, (δT_D−δT_V)/T²), plus K_p |Δu|²/T. The discrete
solver checks this balance, stationary common boosts at zero wave
number, total energy and momentum, analytic relative drift damping,
and nonincreasing entropy norm under the exact matrix exponential.
At low wave number the acoustic speed approaches that of one common
thermal fluid; in the opposite limit the separate chemically frozen
dark and radiation sound speeds reappear.

## A heat-only coefficient misestimates relative motion

Measure time in units Γ_Q⁻¹, with Γ_Q=K_Q/h_D, and use k/Γ_Q=1,
h_V/h_D=10, and initially opposite momentum density but zero
density perturbations. In the two illustrative replacement controls,
keep the exact energy and number sources while either setting
K_p=K_Q or dropping only the elastic contribution to K_p.

| Common T (GeV) | Relative velocity remaining after Γ_Q⁻¹, full | Set K_p=K_Q | Drop elastic momentum only |
|---:|---:|---:|---:|
| .050 | .5265 | .2635 | .7703 |
| .300 | .5416 | .2630 | .7710 |
| .6882487 | .7961 | .2628 | .7964 |

These conditional finite-time amplitudes can depend on the initial
perturbation, wave phase and thermal closure. They demonstrate why
the energy coefficient cannot determine longitudinal momentum
damping even in this simplified model. They are not a matter transfer
function or a bound on the physical dark-matter cutoff.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  charged_lepton_portal_local_longitudinal_reference.py \
  results/charged-lepton-portal-combined-momentum-2026-11-27.json
```
[Full response receipt](results/charged-lepton-portal-local-longitudinal-2026-11-27.json.gz).
