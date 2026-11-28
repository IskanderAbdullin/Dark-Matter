# Tidal removal and escape to infinity carry different particle energies

The [initial fixed-halo ejection](INDEPENDENT-HALO-BOSE-FINITE-EVAPORATION-KERNEL-2026-11-26.md)
uses a sink at the specified cutoff $E_t$.
In a spherical isolated vacuum continuation with
potential zero at infinity, that cutoff is still
negative. The threshold for positive energy is
above it by $GM/r_t$ per mass. In the dimensionless
log-radius Poisson solution this offset
is $dy/d\ln r$ at the tidal radius.
Use the same initial halo DF for the two thresholds.
The isolated-infinity comparison does not grant
a particle escaping to a tidal reservoir the
same outcome as an unbound particle.

For an outgoing high-energy particle with
dimensionless speed-energy parameter $a$ and
polar range $|c|\in [C,1]$, the tested excess
over a chosen threshold $q$ is $a+b|c|-q$,
where $C=(q-a)/b$ restricted to the allowed range.
The bare polar integral of this excess is
$b(1-C)^2/2$ when $0<C<1$.
The quadrature independently checks that identity
and weights the same linear excess with the
surviving bound-side Bose factor. For the tidal
case convert the moment to the zero-at-infinity
convention by subtracting the exterior offset.
Neither final escaped distribution is injected
back into the fixed initial bath.

The recorded energy is the one-particle orbital
energy in the held potential. In the infinitesimal
isolated-removal limit, this is the correct
particle contribution to the total energy change:
the first variation of the Newtonian
self-gravity energy $U=\frac12\int \rho\Phi$
is $\int \delta\rho\,\Phi$, not half that.
But the finite-time potential change or energy
exchange with external tides is not solved.
A tidal energy selection is not itself a proof
of outward crossing of the geometric radius
with arbitrary angular momentum.

| Escape depth | Energy-selected tidal fraction / Gyr | True isolated escape / Gyr | Mean tidal-selected energy relative to infinity, in $mV_0^2$ |
|---:|---:|---:|---:|
| 2 | .006157 | .0008144 | −.224 |
| 5 | .001530 | .0003850 | −.188 |
| 8 | .0001404 | .00006704 | +.150 |

All three examples use peak occupation 7,000.
Only about 13%, 25% and 48% of the tidal-selected
events have nonnegative energy at infinity.
Consequently the mass and the sign of the particle
energy flux need not agree between boundary
prescriptions. The finest independent scrambles
differ by at most .52% on the measured moments;
the coarse-to-fine difference on the deepest
energy moment reaches about 5.1%.
The lower-order tidal rates agree with the previous
independent angular implementation within .1%.

~~~sh
python -W error halo_bose_escape_energy_ledger_reference.py
~~~

[Two-threshold held particle and energy rates](results/halo-bose-escape-energy-ledger-2026-11-26.json).
