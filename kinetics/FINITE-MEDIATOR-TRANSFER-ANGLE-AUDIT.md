# Independent transfer-angle coordinate for full 5D finite-mediator Pauli heat (2026-11-22)

## Pre-registered numerical gate

The existing boosted-CM five-dimensional integral uses direct Gauss-Legendre integration over scattering cosine c, retaining the nonconstant vector spin trace, the outgoing boosted chi/electron energies, all four Fermi factors and signed energy transfer. A finite spacelike vector mediator multiplies the trace by [mV²/(mV²−t)]². For mV small relative to CM momentum p*, direct cosine quadrature undersamples the forward region. Hypothesis: a second *coordinate system* for the same integral gives compatible Q at adequate angular resolution and better forward resolution. Failure of agreement is a numerical error, not a physical heat signal.

Define x=|t|=2p*²(1−c), u=x/(mV²+x), u_max=4p*²/(mV²+4p*²). Then c=1−mV²u/[p*²(1−u)] and dc in the reversed angular orientation is mV² du/[p*²(1−u)²]. The propagator is exactly (1−u)², cancelling the explicit Jacobian denominator in its product; the full spin trace, boosted final energies and Pauli factors must still be evaluated at c(u), and **cannot** be replaced by their vacuum angular averages. At p*=0 the limiting collision phase-space measure vanishes; implementation uses guarded division.

Implemented as optional `angular_scheme="transfer"` in `kinetics/independent_massive_chie_5d.py`, leaving `cosine` default and contact benchmark unchanged. Tests `tests/test_finite_mediator_transfer_quadrature.py` compare both coordinate rules at mV/T=1,3, massless/massive leptons and thermal fixed-point cancellation. Both share radial and azimuth grids and the spin-trace implementation, so this is **not** a fully independent matrix-element calculation or a continuum accuracy certificate.

## Status and negative controls

Previous HEAD 23386aa scientific verification [run 44975219960](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44975219960) completed successfully. New implementation/test commits 2032662 and 6ffb9bf were pushed; their CI and numerical results were not yet observed at this checkpoint. No numerical Q agreement or convergence is claimed before actual tests. A single fixed-mass or equal-temperature test cannot bound the error of signed Q, since positive and negative event contributions can cancel. A valid error budget needs independent angle, radial, azimuth, cutoff, mediator-mass and bath-temperature scans, with absolute as well as relative errors.

Next executable check: `python -m pytest -q tests/test_finite_mediator_5d.py tests/test_finite_mediator_transfer_quadrature.py`, inspect real public CI, then add a compact reproducible angular/cutoff scan and independent numerical t-integration. No physical gamma_el, 86-eV relic, BEC, Omega or CLASS conclusion follows.
