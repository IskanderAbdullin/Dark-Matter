# CLASS long-wave normalization recovered, 2026-11-25

[Independent control 45624656360](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45624656360) uses the same actual 480-cell energy-closed kinetic spectrum and equal-density CDM reference. We varied gauge, cold quadrature (uniform ascending vs transformed descending), and replaced the narrow cold tracer by true CDM while retaining the hot tracer.

| configuration | P/P_CDM at k=.01 | at .5 | at 1 Mpc^-1 | sigma8 |
|---|---:|---:|---:|---:|
| Newtonian two tracer | 17.8525 | 2.9081 | 1.3862 | 2.1533 |
| Newtonian hybrid | 1.0139 | .6775 | .3117 | .8104 |
| synchronous two tracer | .99986 | .67641 | .31081 | .80777 |
| synchronous hybrid | 1.00001 | .67700 | .31116 | .80803 |
| synchronous uniform cold | 1.00055 | .67943 | .31260 | .80912 |

Background omega_m=.14243911 and omega_ncdm=.12006910 match between gauges for two tracers. The pathological result is specifically associated with the Newtonian gauge/ultracold tracer combination, not a physical enhancement or simply quadrature error. The exact internal origin is not yet established. Synchronous is the accepted numerical branch because the large-scale normalization, independent cold-CDM replacement and alternate quadrature all agree. The hybrid is valid only on scales where cold free streaming is negligible.

Switch production transfer runs to synchronous gauge and assert max absolute deviation from unity below k=.03 Mpc^-1 <5% on every run. Do not use earlier Newtonian half-power values or conclusions. New full transfer and mesh convergence runs are required. The perturbations are still conditional on adiabatic initial conditions after freeze-out, collisionless species, and stable massless daughter radiation; no observational likelihood is implied.
