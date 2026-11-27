# The early spectra are still being produced

Retain all three comoving momentum distributions from the
[conservative live-parent calculation](INDEPENDENT-EARLY-POLE-MASSIVE-DAUGHTER-ENERGY-2026-11-26.md)
and record both occupation and scalar number per shell. The visible
bath here freely redshifts, and H includes that bath and all tracked
particle energies. Shell brackets do not assume an internal
sub-cell occupation shape.

At x=0.05 with no parent mixing, the computed scalar yields on the
third grid (97 scalar, 156 daughter, 176 parent cells) are
0.000346861 (massless fermion daughter), 0.000364093 (mass ratio
0.06), and 0.000380008 (mass ratio 0.30). The massive benchmarks
are 5.0% and 9.6% above the exact third-grid massless null. This
refines the comparisons on the smaller grids.

All three third-grid cases put more than 99.998% of the scalar
count below q=10⁻⁵, yet have a momentum RMS near 0.0046–0.0049
in the same units q=p/T. A very sparse hot tail carries the second
moment, and the self-consistent Hartree contribution keeps the
homogeneous scalar equation of state close to one third. The
small momentum of most particles does not establish dust-like
stress or a physical matter transfer function at this time.

Most importantly, the scalar number is still changing: at x=0.05
the net d ln Y/d ln a is 1.71–2.19 on the third grid. For mass ratio
0.06 the next-coarser value is 4.46; the sign and order-one rate
matter here, not convergence of the instantaneous derivative.
On that second grid an illustrative extension to x=0.10 increases
the count by about 5.9 times. Do not export x=0.05 as a frozen
relic to CLASS. The extension remains the same conditional
effective-kinetic model, not a completed freeze-out calculation.

~~~sh
python -W error early_pole_spectrum_snapshot_reference.py --ratio .06
python -W error early_pole_spectrum_snapshot_reference.py --ratio .06 --scalar 95 --fermion 125 --parent 108
python -W error early_pole_spectrum_snapshot_reference.py --ratio .06 --xend .1
~~~

Three matched third-grid snapshots:
[massless](results/early-pole-spectrum-third-null-2026-11-26.json.gz),
[ratio 0.06](results/early-pole-spectrum-third-m06-2026-11-26.json.gz),
[ratio 0.30](results/early-pole-spectrum-third-m30-2026-11-26.json.gz).
Independent checks: [refined ratio 0.06](results/early-pole-spectrum-fine-m06-2026-11-26.json.gz),
[fast mixing](results/early-pole-spectrum-fine-fast-2026-11-26.json.gz),
[extension to x=0.10](results/early-pole-spectrum-fine-extend-2026-11-26.json.gz),
and [matched massless extension](results/early-pole-spectrum-fine-null-extend-2026-11-26.json.gz).
