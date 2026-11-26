# Momentum-resolved check of late elastic exchange

A source-free fixed-comoving-shell evolution starts from the separate
thermal shapes inferred by the earlier late closure at b=17.714657
relative to the 12 MeV snapshot. It evolves projected number and
kinetic-energy moments of both species without keeping them thermal
after initialization. Only the below-pole nonlinear angular elastic
operator is included. The initialization is inherited from the thermal
test and is not a separate validation of the earlier preparation.
The Hubble comparison includes only freely redshifted photons as before.

A common-temperature detailed-balance reference factors the exponential
collision residual. It need not itself match the actual two-species
energy or shape. Newton starts from the current nonthermal physical
state, and the residual is normalized by the actual old or trial counts.
Separate shell counts, pressure work and the equal-and-opposite energy
ledger are checked. Initial scalar count loss at the high cut is below
1e-12 and fermion loss is about 1.47e-5.

At b=410 the elastic shell runs give scalar moment-equivalent
temperatures 10.86 to 11.14 keV and fermion temperatures 0.459 to
0.469 keV. The 32-interval separate thermal calculation from the same
starting temperatures gives 10.69 and 0.479 keV. Under the scalar
moment diagnostic, the count occupies about 0.00427 to 0.00462 of
continuum saturation, versus 0.00483 in the matched thermal closure.
On the 12-shell, 16-step test the initial comoving kinetic energy
normalized scalar heat exchange is -0.0431; the matched thermal
16-step value is -0.0468. Rate/cumulative heat is more sensitive
to angular and shell refinement than the bulk normal-state comparison.

The integrated scalar count in the first bin (q/T_R below .005) is
about 2.01e-7 with or without exchange on this mesh. It is a finite-volume
count, not a p=0 occupation or a test of singular condensation.
Below-pole exchange remains weak enough to give distinct effective
temperatures at this stage, for both closures. These tests do not solve
the full late cosmology or other possible number/energy-transfer processes.

## Reproduce

Prepare the external `EOS.dat` with the [shared benchmark setup](INDEPENDENT-BENCHMARK-SETUP.md).

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1
B=results/12-quantum-48-2026-11-26.json
S=results/12-quantum-48-shape-2026-11-26.json
A=results/late-thermal-anchor-2026-11-26.json
python -W error late_two_temperature_expansion_reference.py EOS.dat "$B" "$S" --steps 36 > anchor.json
python -W error late_soft_two_moment_test_reference.py EOS.dat "$B" "$A" --steps 16 --mesh 12
python -W error late_soft_two_moment_test_reference.py EOS.dat "$B" "$A" --steps 12 --mesh 16 --orders 4 3 2 10
python -W error late_two_temperature_expansion_reference.py EOS.dat "$B" "$S" --steps 8 16 32 --start 17.71465708816752 --restart-temps 0.0002888504771421622 0.00018112220246041004
```
