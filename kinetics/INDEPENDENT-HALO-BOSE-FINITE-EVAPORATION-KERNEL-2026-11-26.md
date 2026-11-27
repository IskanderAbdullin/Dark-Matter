# Finite lowered Bose profiles leak under direct elastic collisions

The [finite lowered Bose halos](INDEPENDENT-HALO-BOSE-FINITE-LOWERED-MIXING-2026-11-26.md)
are steady solutions of the **collisionless** spherical
equations. Cutting the energy distribution at the
escape boundary does not make it an elastic thermal
fixed point. Here compute the initial collisional
loss of bound particles into an empty unbound
bath, with the same illustrative $R,V_0$ and
vacuum cross section as the [orbital benchmark](INDEPENDENT-HALO-BOSE-ORBITAL-ELASTIC-CLOCK-2026-11-26.md).
Keep the potential and distribution fixed.

At dimensionless radius with local remaining escape
energy $W-y$, the incoming speeds $v,u$ are both
bound. The outgoing kinetic energies are
$(v^2+u^2)/4\pm Vg c/2$.
The upper daughter exceeds the escape energy only
when $|c|>C$, with
$C=[W-y-(v^2+u^2)/4]/(Vg/2)$.
Since both incoming particles were bound, both
outgoing particles cannot simultaneously escape.
The vacuum ejection probability for an isotropic
contact collision is $1-\min(1,\max(0,C))$.

For final-state statistics, give the escaped leg
zero outgoing occupation and keep $1+f$ of the
surviving bound partner. Its Bose-plus-lowering
factor has an elementary polar integral, using
$\int dx/(e^x-1)=\ln(1-e^{-x})$.
Independent polar quadrature checks this identity.
The global rate is the local pair rate
$n^2\sigma\langle g P_\mathrm{esc}\rangle/2$
integrated over halo volume and divided by
total bound count. The incoming identical-pair
division by two is needed because an ejection
removes only one bound particle. Separate radial
and momentum/scramble controls check the result.

| Escape depth, peak occupation 7,000 | Vacuum initial fractional loss / Gyr | With final bound-side Bose factor / Gyr |
|---:|---:|---:|
| 2 | .003193 | .006158 |
| 5 | .001209 | .001530 |
| 8 | .0001269 | .0001404 |

The rates here use the illustrative $1.25\,{\rm cm^2/g}$
vacuum cross section and scale linearly only while
the held profile is unchanged. Bose enhancement of
direct escape is at most about two in these arms,
much smaller than tagged low-speed central encounter
enhancement. The deepest model retains its particles
more effectively in this immediate sink diagnostic.
The largest combined grid and scramble difference
was under .24%; it is not an uncertainty for actual
self-consistent halo evolution.

This is a direct large-energy-transfer collision
kernel, not a solution for diffusion in escape
energy or a time-dependent evaporation lifetime.
Returning particles, evolving temperature,
tidal stripping, slow transport, baryons and
a feedback response of the potential are not
included. If there is no initially empty
unbound region, the inverse flux would need
to be calculated too.

~~~sh
python -W error halo_bose_finite_evaporation_kernel_reference.py
~~~

[Full-volume vacuum and stimulated initial losses](results/halo-bose-finite-evaporation-kernel-2026-11-26.json).
