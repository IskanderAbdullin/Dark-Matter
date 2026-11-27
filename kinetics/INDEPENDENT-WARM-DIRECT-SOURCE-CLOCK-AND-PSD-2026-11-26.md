# Separate temperature residency from the direct scalar momentum distribution

The [vacuum lepton and measured pion bands](INDEPENDENT-WARM-LATE-COLD-EXACT-TRAVEL-BANDS-2026-11-26.md)
share a 12–50 MeV adiabatic radiation history. Resolve
where the cold yield is produced by splitting the time
integral at 25, 35, 40, 45, 47.5, 49 and 50 MeV. Each
term is positive. The actual travel cut is set at each
birth; for the equal-temperature label the same comoving
cut is used throughout.

| Selection | Fraction from 45–50 MeV | Fraction from 49–50 MeV | Yield if duration above 45 MeV is reduced to one fifth |
| --- | ---: | ---: | ---: |
| p/T≤.1 at equality | .85567 | .29156 | .31547 of standard |
| Path ≤.3 Mpc | .99739 | .64847 | .20209 of standard |
| Path ≤1 Mpc | .83533 | .27509 | .33173 of standard |

Two independently sampled temperature/rate grids agree.
At the same peak temperature, replacing H with the
prescribed H_rad sqrt(1+(T/30 MeV)²) reduces these
three yields to .5353, .5217 and .5367 of standard;
using 20 MeV instead gives .3893, .3775 and .3906.
These parameters specify clock controls, not simulated
reheating or an independently verified extra fluid. For
these controls the geometric selection and entropy curve
are held fixed. A different expansion during subsequent
travel would require a new distance calculation. The
point of the comparison is that the peak alone does
not fix the deposited yield.

## Positive momentum input for the later response

Preserve the entire redshifted scalar spectrum from
these two direct sources as positive comoving shell
counts. At each kinetic sample the unpolarized scalar
energy is uniform over the parent-dependent lab
interval. Intersect that interval exactly with the
momentum cuts, rather than binning a sampled decay
angle. Use independently refined time, invariant
energy, boost and initial thermal statistics quadratures.
No thermal shape is imposed on the produced scalars.

The total scalar yields are 6.14738e−6 (leptons) and
5.52000e−7 (measured pion interval) in the finer
quadrature. The maximum difference in either normalized
CDF between the two controls is 7.82e−6. At q/T0=.1
the independently integrated cold value is recovered
within .018% (lepton) and .00019% (pion), despite the
very small counts there. The q/T0>18 overflow is
stored in a separate positive bin rather than removed.

The result is a conditional *input* for collisionless
tests. It assumes both channels are vacuum sources,
contains no pre-existing phase, omitted hadronic
bands or rescattering, and does not establish a relic
density or independently generated gravitational metric.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
D=results
python -W error warm_late_cold_source_temperature_localization_reference.py EOS.dat KLOExSec.data \
  "$D/warm-late-cold-visible-exact-recoil-2026-11-26.json" "$D/warm-late-cold-exact-travel-bands-2026-11-26.json"
python -W error warm_direct_vacuum_scalar_psd_reference.py EOS.dat KLOExSec.data \
  "$D/warm-late-cold-visible-exact-recoil-2026-11-26.json" "$D/warm-late-cold-visible-kloe-recoil-2026-11-26.json"
~~~

[Temperature residency and clock choices](results/warm-late-cold-source-temperature-localization-2026-11-26.json);
[positive comoving input spectra](results/warm-direct-vacuum-scalar-psd-2026-11-26.json).
