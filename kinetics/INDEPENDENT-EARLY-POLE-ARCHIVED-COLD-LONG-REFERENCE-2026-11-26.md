# Longer cold-window exposure does not establish freeze-out

Repeat the [matched cold-window comparison](INDEPENDENT-EARLY-POLE-ARCHIVED-SHARED-COLD-WINDOW-2026-11-26.md)
at two longer exposures, with the wider photon/fermion upper support
4/5, actual archive splits near the hard aperture, and the same
initial states for the decay-only and combined runs. The exterior
fermions and photons are not evolved here, and the scalar
spectators retain their contribution to the Hartree energy.
Both rate coefficients continue to be the finite-box vacuum
references, not real in-medium production rates.

| Fixed-scale exposure | Fine gain, decay-only / combined | Third gain, decay-only / combined |
| ---: | ---: | ---: |
| .0003 | 2.2771e−5 / 1.9987e−5 | 2.0278e−5 / 1.7369e−5 |
| .0008 | 5.6964e−5 / 4.0826e−5 | 4.8015e−5 / 3.2280e−5 |

At .0003 the suppression is 12.23% (fine) and 14.35% (third).
Doubling from 8 to 16 scalar cells on the fine archive changes
the value to 12.18%; the absolute gains move by about .5%.
At .0008 the suppression rises to 28.3% and 32.8% on the
two eight-bin references. Nevertheless the combined scalar
count increases in every reported subinterval, including the
last one; increased suppression is not a production stop or
freeze-out. The strict total energy, number, and entropy checks
remain active.

**Separate redshift risk.** Hold all original comoving occupations
and the conformal scalar pole fixed, and move only the bare and
asymptotic hard dispersion with scale. The representative hard
parent aperture shifts by roughly 0.35 of its starting half-span
by ln(a/a₀)=.0008 and by 1.76–1.77 half-spans by .004. At the first
of those times the frozen cold-source change is about +2% on the
fine archive but +24% on the third. These are not predicted
time-dependent rates: they are an indicator that the saved local
fermion troughs and the moving aperture matter. A longer fixed-scale
box cannot substitute for joint redshift, heat transfer and actual
medium evolution.

Reproduce the exposure comparisons with
early_pole_archived_shared_cold_window_reference.py --grid fine/third
--wide 1 --refine-aperture --exposure .0003 or .0008.
The static aperture comparison is
early_pole_frozen_archive_redshift_aperture_reference.py.
Raw: results/early-pole-archived-cold-long-reference-2026-11-26.json.gz.
