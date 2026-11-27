# Positive-energy escape cannot follow a fixed-shape lowered Bose sequence

The [isolated particle-energy flux](INDEPENDENT-HALO-BOSE-ESCAPE-ENERGY-LEDGER-2026-11-26.md)
removes positive energy when a particle escapes to
infinity. This supplies an instantaneous consistency
test for an isolated quasi-static sequence of
one-component lowered Bose equilibria.
No feedback evolution is computed.

At fixed particle mass and fixed dimensionless
occupation and escape depth $(F,W)$, the central
mass density scales as $\rho_0\propto V_0^3$
from the phase-space measure. The normalized
Poisson scale therefore scales as $R\propto
V_0^{-1/2}$, and the total bound mass has
$M\propto V_0^{3/2}$.
Each finite profile has vanishing edge pressure;
the isolated virial identity is $2K+U=0$.
Hence $\mathcal E=K+U=-K\propto-M V_0^2$
and
$$
 |\mathcal E|\propto M^{7/3}
 \quad\text{along a fixed-shape sequence.}
$$
Such a sequence becomes **less** bound in total
energy as it loses particles. Positive-energy
escape instead removes energy and initially
increases $|\mathcal E|$ in isolation. These
signs are incompatible. No adjustment of just
the dimensional scale keeps both the mass flux
and energy flux consistent at fixed $F,W$.
A physical remnant would need shape changes,
departures from this equilibrium family, or
additional physics such as a non-isolated energy
budget. The tidal-removal example does not obey
the isolated positive-energy premise.

The driver independently integrates both kinetic
energy and the full isolated potential-energy
integral rather than assuming the virial identity
to get both sides. It also looks at two algebraic
first steps inside a wider family: vary W while
holding F, or change F at fixed W. The latter
inversion is ill-conditioned at large F, where
density and total mass are close to their critical
limits; the result is an exploratory diagnostic,
not a recommended evolution. Varying W is
better conditioned, but also does not assert
that the gas remains on this family.

For $W=5,F=7000$, the isolated clock gives
$d\ln|\mathcal E|/dt=+2.16\times10^{-4}\,\mathrm{Gyr}^{-1}$,
where fixed shape would give $-8.98\times10^{-4}\,\mathrm{Gyr}^{-1}$.
If only $W$ can change at fixed $F$, the algebraic
local tangent has $dW/dt\simeq-0.00261/\mathrm{Gyr}$
and $d\ln V_0/dt\simeq+0.000438/\mathrm{Gyr}$.
All three sampled widths give a decreasing $W$
and increasing $V_0$ in this assumed sequence.
The largest independent virial integral defect
is below $2.7\times10^{-6}$; the fixed-$W$
inversions have condition numbers above $10^5$,
so they are not used to infer a physical change
in central occupation.

~~~sh
python -W error halo_bose_lowered_energy_family_reference.py \
 results/halo-bose-escape-energy-ledger-2026-11-26.json
~~~

[Virial checks and conditional tangent data](results/halo-bose-lowered-energy-family-2026-11-26.json).
