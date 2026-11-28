# A conditional scalar-only equilibrium capacity, not condensation dynamics

Use the scalar number and conserving Hartree energy from the
[early spectra](INDEPENDENT-EARLY-POLE-SPECTRAL-SNAPSHOT-2026-11-26.md)
as fixed inputs to an auxiliary equilibrium problem. Ask whether a
spatially homogeneous, ergodic scalar-only sector can match those two
moments if its elastic scattering is fast and conserves scalar
particle number. None of those dynamical assumptions was used in
the production solver.

At each positive scalar quasiparticle mass m, integrate a Bose
distribution at its limiting chemical potential μ=m. Use
E−m=p²/(E+m) when evaluating low momenta, to avoid cancellation.
The Hartree gap includes the excited inverse-energy moment and,
on the branch with a formal zero mode, n₀/m. Energy includes the
quartic counterterm and m n₀ for the zero mode. First match number
on the saturated regular branch to obtain its reference energy.
Regular solutions at nine fugacity offsets on the same number
slice have larger energies in these tests. Then solve separately
for a zero-mode population and positive excitation temperature at
the actual energy. Orders 96 and 192 check momentum quadrature;
the weak-coupling, small-vacuum-mass limit also reproduces the
analytic ideal-gas saturation constants in an independent check.

At x=0.05, the three retained third-grid spectra have actual scalar
energy about 0.047% of the saturated regular reference at the same
number. The formal zero-mode fraction that matches both moments
is 99.984–99.986% in this auxiliary model. This is not a measure
of how much is already condensed, how fast it would condense, or
the stability of a coherent field. A real scalar's particle number
does not have an exact symmetry protecting it. The driven system
may remain nonergodic and the Hartree ansatz need not model a
coherent condensate. The closure should not be used as a CLASS
input or as evidence for late Bose stars.

As a separate spectral observation, the fraction with q≥0.1 in
these early runs is about 8–11 parts per million, yet it supplies
more than 99.98% of the second momentum moment. A converged
early total count cannot replace a convergence test of its hot tail.

W. Yin's [decay production study](https://arxiv.org/abs/2301.08735)
motivates the early burst; the frozen equilibrium diagnostic above
does not follow from that decay calculation. Sakurai and Yin's
[scattering study](https://arxiv.org/abs/2410.18968), Sec. IV.4,
discusses possible later redistribution and a zero mode in a distinct
reaction. Neither is a computation of the relaxation rate or the
coherent phase for the present model.

~~~sh
python -W error early_pole_scalar_condensate_capacity_reference.py --ideal-control
python -W error early_pole_scalar_condensate_capacity_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Raw capacity results](results/early-pole-scalar-condensate-capacity-2026-11-26.json).
