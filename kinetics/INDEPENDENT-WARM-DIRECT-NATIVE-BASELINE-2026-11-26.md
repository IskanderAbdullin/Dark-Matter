# Separate background effects from the native response of a small scalar component

Use the [positive direct vacuum shell input](INDEPENDENT-WARM-DIRECT-SOURCE-CLOCK-AND-PSD-2026-11-26.md)
with its conditional 2.0992e−6 scalar share of dark
mass when the held stable count is used. This note sets
up an independent matched CLASS test and bounds the
homogeneous background before assigning any matter-power
difference to the physical component.

At fixed present density the nonthermal scalar replaces,
rather than augments, the corresponding CDM density.
Bound its kinetic energy in each positive q shell by
monotone endpoint energies. Use the compensated identity
sqrt(1+x²)−1=x²/(sqrt(1+x²)+1) instead of subtracting two
nearly equal energy integrals. The tiny sampled overflow
is enclosed by a conservative q/T0≤1300 for the finite
source quadrature: total lab energy is no larger than
max(60T,V+30T)+75T throughout 12–50 MeV. This is a bound
on the represented numerical source, not an absolute
momentum cutoff on a thermal bath.

If all those particles are already deposited at 50 MeV
in the reference model, the early relative Hubble
difference lies in [1.6059e−8,1.8333e−8] there. Near
equality it is [6.255e−11,8.152e−11]. The positive
relative equality scale shift is [2.502e−10,3.261e−10]
under the stated rest/kinetic equality definition.
The conformal time shift to reference equality is
between −8.616e−8 and −6.831e−8 Mpc. The separate
relativistic-density equivalent, expressed in standard
neutrino units, is [1.8457e−7,2.1071e−7]. These are
conditional background comparisons and do not presume
that particles really exist before their production.

## Matched response preparation

Use synchronous gauge, fixed cosmological parameters,
omega_cdm+omega_ncdm=.120, scalar rest mass 86 eV,
and a common scalar comoving input. The other dark
matter is a CDM proxy. Normalize input PSD tables
positively, preserving the template count and squared
momentum for both tested smoothing widths. The CLASS
cut q/T0=29 contains more than .99999997 of each
preserved squared moment. A Fermi distribution with
the direct shape's same squared momentum tests whether
that moment alone determines the response at equal
abundance.

To separate a physical effect of order a few parts
per million from numerical solver changes, retain both
a zero-species baseline and a trace-species baseline
with exactly the same non-cold phase sampling, as well
as the physical fraction and deliberately amplified
fractions. Use tighter tolerance and hierarchy/momentum
sampling as independent controls. The preparation does
not itself certify small power differences, numerical
convergence, a full nonthermal production model or an
observational constraint.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_direct_vacuum_native_class_reference.py EOS.dat \
 "$D/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz" \
 "$D/warm-direct-vacuum-scalar-psd-2026-11-26.json" \
 "$D/warm-direct-vacuum-fixed-metric-response-2026-11-26.json" native > native-spec.json
python -W error warm_direct_vacuum_background_response_bound_reference.py EOS.dat \
 "$D/warm-direct-vacuum-scalar-psd-2026-11-26.json" "$D/warm-direct-vacuum-fixed-metric-response-2026-11-26.json"
~~~

[Background difference bounds](results/warm-direct-vacuum-background-response-bounds-2026-11-26.json);
[response input/specification](results/warm-direct-vacuum-native-class-run-spec-2026-11-26.json).
