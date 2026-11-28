# Closed-path spacelike fluctuations in the soft equilibrium reference

This is a finite prescribed-current calculation. It is not a
self-consistent neutral Yukawa vertex, a physical linewidth, a new
cosmological trajectory, or a source ready for CLASS. We ask whether
separating two charged lines gives the same phase covariance as treating
them independently.

We close the two straight real-time paths with a spatial endpoint link.
The Fourier transform on a segment from a to b is the displacement times
exp(i K·midpoint) times sinc(K·(b−a)/2). Contracting with K gives the
difference of endpoint phases divided by i. These cancel on the triangle;
the cancellation is less than 3.6e−14 relative in an independent complex
current check. The same test with the endpoint connector omitted fails;
the identical-path test returns exactly zero. Pure-gauge longitudinal
pieces proportional to K therefore decouple from the closed test current.
The specific path geometry is still a choice and is not derived from the
production operator.

We project its temporal component and spatial components perpendicular
to the exchanged momentum on the classical symmetrized equilibrium
spacelike hard-thermal-loop spectra. The fixed-cutoff rectangular check
at k/mD ≤ 20 and r mD=.6010 has the independently integrated static
longitudinal limit .158560; the finite-time fine result at t mD=160 is
within 0.4% of it. Taking the infinite-cutoff limit first would compare
instead to .161032. The squared current enters before the frequency
integral, with both spacelike sides included. Normalization is chosen to
recover the static reference of [Laine et al.](https://arxiv.org/abs/hep-ph/0611300);
see also [Caron-Huot](https://arxiv.org/abs/0710.5726) for real-time HTLs.

For the triangular test at t mD=996, parent speed .9982 and final
separation .6010/mD, the normalized spacelike contributions are:

| k cutoff / mD | Longitudinal | Transverse |
| ---: | ---: | ---: |
| 8 | .06760 | .02446 |
| 20 | .07761 | .02929 |

At cutoff 20, the last two quadratures change these terms by about
2.5e−6 and 1.1e−8, respectively. The cutoff variation is much larger than
the quadrature variation. The local Doppler approximation with a
continuously growing separation gives similar but not identical terms
and does not implement the full endpoint condition. A model at rest
cannot replace the moving transverse contribution.

The spacelike reference omits timelike collective poles and cannot by
itself establish their importance. It also omits hard-transfer matching,
fermion recoil, changes of the path by collisions, in-medium production
vertex matching and nonthermal gauge statistics. No coupling threshold
follows from this comparison.

Run from kinetics/ with OPENBLAS_NUM_THREADS=1 and
early_pole_finite_wilson_htl_reference.py --level 0, 1 or 2.
The closed-current, static shape, and local Doppler scripts provide
independent checks. Raw combined data are in
results/early-pole-finite-wilson-spacelike-2026-11-26.json.gz.
