# Live number and energy moments in a hard 1↔2 aperture

The frozen [subcell moment check](INDEPENDENT-EARLY-POLE-GAUGE-VOLUME-MOMENTS-2026-11-26.md)
is extended to a reversible, evolving decay window. Each hard fermion
cell keeps count and a centered bare-energy moment, using a local
Fermi function from the two moments. The full gauge energy is
evaluated from these reconstructed distributions with a fixed
spectator photon moment. Differentiate through the local fit:
the gauge moment derivative is the weighted covariance row times
the inverse two-by-two occupation Hessian. It is not a shift at a
single cell center. The scalar mode uses an optional live Hartree
energy derivative. With it, the interaction energy is kept only once.

Each active parent/daughter cell pair uses the exact on-shell triangle
domain for the piecewise variational dispersion. A positive quadrature
in the parent momentum deposits count and centered energy in both hard
cells. The forward and inverse occupations are evaluated from those
same local fits. The independently accumulated entropy production
agrees with entropy evolved from the reconstructed states. The
closed-energy model has zero collisional work to rounding; using bare
events despite the stored gauge energy gives a small nonzero work
that accounts for the measured continuous-time energy drift.

The local Pauli fitting helper independently tests discrete quadrature
moment feasibility up to 0.1% of either extremal energy at fixed count;
the numerical moment error is below 9e−16 in that control. This does
not guarantee that arbitrary larger steps remain inside the moment
domain. The live benchmarks instead test the occupied states throughout
a finite short exposure.

**Initialized from the archived x=.05 spectra.** In each grid, only the
hard radial window 1.30–1.56 and the scalar bin containing q/T=4.1e−6
are evolved. Other scalar bin counts are fixed spectators but their
contribution enters the live Hartree gap and energy. Outside-window
charged contributions and the hard photon moment remain fixed. The
exposure is divided by the initial occupation of the chosen scalar
bin, separately in each archive. It is a numerical matched clock,
not a physical cosmological time.

| Frozen initialization | Bare-energy event gain | Variational event gain | Final energy drift, variational |
| --- | ---: | ---: | ---: |
| Fine archive, 33 subcells per species | 9.9242e−5 | 1.0307e−4 | 0 within rounding |
| Third archive, 48 subcells per species | 7.6089e−5 | 8.5612e−5 | −3.5e−18 |

The dressed/bare difference is about 3.9% and 12.5% in these two
volumes; the bare-event energy drifts are −5.41e−11 and −4.15e−11.
Refining the same archived cell representation from factor 2 to
factor 3 changes the dressed gain by .29% and .86%, respectively.
This separates local reconstruction refinement from the mismatch
between the original global grids and scalar pole/occupancy.
The synthetic broader volume is more sensitive to cell closure;
its coarse variations should not be called a physical correction.
Neither comparison fixes the microscopic vertex, expansion, photon
pair rate, or actual late scalar spectrum.

Reproduce with early_pole_cell_pauli_entropy_reconstruction_reference.py,
early_pole_live_pauli_gauge_moment_decay_reference.py and
early_pole_archived_pauli_moment_window_reference.py.
Results: results/early-pole-live-pauli-gauge-moment-decay-2026-11-26.json.gz.
