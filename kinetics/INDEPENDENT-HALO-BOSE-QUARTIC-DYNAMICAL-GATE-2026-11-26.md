# Cross section needed for rapid local halo kinetics

The [local shear and speed-shape collision brackets](INDEPENDENT-HALO-BOSE-LINEARIZED-KINETIC-CLOCKS-2026-11-26.md)
scale linearly with the vacuum cross section if the
held halo, microscopic angular kernel, and dispersions
are unchanged. Translate the ratios into a comparison
for the constant real-scalar Born quartic
$V=\lambda\phi^4/4!$, for which
$\sigma/m=\lambda^2/(128\pi m^3)$ in the
nonrelativistic limit with the identical-final factor.
This does not use large gross stimulated event counts
as a substitute for kinetic response.

For a local fluid description on the circular
dynamical timescale, rates in every pertinent
nonconserved direction would need to exceed
$\Omega_c=V_c/r$, usually by a substantial factor.
A chosen test direction with a smaller Rayleigh
quotient is an obstruction to that uniform local
timescale ordering for this fixed collision operator.
Matching the tested rates alone is not sufficient.
This comparison is distinct from slow convergence
over an assumed galaxy age.

| Radius, peak occupation 7,000 | $\sigma/m$ for tested shear to equal $\Omega_c$ | $\sigma/m$ for tested isotropic shape to equal $\Omega_c$ |
|---:|---:|---:|
| $.1R$ | $17.2\ {\rm cm^2/g}$ | $41.9\ {\rm cm^2/g}$ |
| $R$ | $27.3\ {\rm cm^2/g}$ | $54.9\ {\rm cm^2/g}$ |
| $10R$ | $242\ {\rm cm^2/g}$ | $367\ {\rm cm^2/g}$ |

The [separate upstream quartic ledger](https://github.com/IskanderAbdullin/Dark-Matter/blob/90c141722376dab64e4cd10422d545474817787c/kinetics/PHI-BOSE-QUARTIC-CURVATURE-HALO-SCATTERING-TENSION-2026-11-24.md)
uses an illustrative reference $5\ {\rm cm^2/g}$
and a conditional static plasma curvature value
$\lambda=1.08516\times10^{-6}$. With the same
Born convention these correspond to about $5$
and $1006\ {\rm cm^2/g}$ for an $86$ eV scalar.
At $5\ {\rm cm^2/g}$, the two test rates at
$.1R$ are only .29 and .12 of the circular frequency;
the much larger static-curvature choice passes that
local comparison but does not validate the assumed
halo profile, early-time plasma, or external
observational constraints. These are independent
benchmarks, not a new cluster likelihood or a
universal limit on self-interactions.

~~~sh
python -W error halo_bose_quartic_dynamical_gate_reference.py \
 results/halo-bose-kinetic-clock-comparison-2026-11-26.json
~~~

[Local orbit-clock and quartic unit conversion ledger](results/halo-bose-quartic-dynamical-gate-2026-11-26.json).
