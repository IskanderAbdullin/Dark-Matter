# Stationary reversible radial jumps versus their diffusion limit

Build an explicitly reversible momentum-jump generator with reflected or
omitted out-of-domain links and a stationary radial measure p² exp(−E/T).
Detailed balance is enforced by symmetric link conductances in the
weighted representation. With diffusion generator F, local detuning d,
readout regulator epsilon and real normalized probe u, solve
y=(epsilon−F+i(d−omega))⁻¹u. The spectral readout Re(u†y) equals
y†(epsilon−F)y and is nonnegative. This positivity follows from
stationarity and a dissipative generator, not from treating an initialized
two-path coherence as a spectrum.

**Imposed nearest-neighbor witness.** With radial step size h, we
scale the per-link rate as h⁻² for a matched diffusion refinement.
For the broad probe on a .15 T radial domain, the spectral central
inverse-height proxy (unit area, not HWHM) is 7.0255, 7.0325, 7.7316,
and 3.1979 at the four imposed rate normalizations .1, 1, 10, and 100.
The 27/39 cell difference is below .003, and 8/16 angular nodes
at fixed 39-cell radial resolution agree to roundoff at angular
redraw rate 2. The unnormalized line-area check gives π.
With the narrower .008 T endpoint weight, the trend differs:
the large-rate proxy depends materially on the reflected domain
(9.59–19.73 across three domains), while the low-rate results agree
to a few parts in 1,000. A width extracted from a selected momentum
window is not automatically a transport coefficient.

**Nonlocal soft-kick witness.** Form symmetric radial jump taps from
the spherical soft HTL transfer kernel of the [kick comparator](INDEPENDENT-EARLY-POLE-SOFT-KICK-CHARACTERISTIC-2026-11-26.md).
Each transfer is split between adjacent discrete jumps with positive
weights conserving its second moment; bulk normalizations agree
within 2e−13. The comparison local generator places the same second
moment on nearest neighbors. Both use the same stationary weights and
endpoint probe. The angular dipole proxy is only about 1e−5 in units
of the inverse angular time for the strongest benchmark. It is omitted
in the matched radial calculation; the omission in the corresponding
toy angular redraw changes a unit-norm correlation by no more than
4.5e−5 over four such times.

For g=.005, cutoff 20 mD, and a .008 T endpoint weight, the magnitude
of the stationary correlation at one inverse angular half-span is
.42060 without collisions, .40410 with nonlocal jumps, and .39602 with
matched diffusion. At two half-spans they are .02844, .02485 and .02252.
The 83/119-cell radial refinement changes the one-time collisional values
by about 6e−5 and 1.8e−4 respectively; comparison with 49 cells exposes
larger discretization effects. With a common positive readout regulator
.35, the central spectral values are .7780, .7621 and .7543 on the
119-cell grid. The regulator is a numerical readout, not a vertex width.
The difference shows that preserving the bulk momentum variance alone
does not preserve a selected stationary resonance response.

These are prescribed classical stochastic generators in an equilibrium
reference. They do not include the full two-line quantum operator,
self-consistent nonequilibrium polarization, hard transfer matching,
actual flavor heat exchange or the scalar production rate. Do not add
the same classical soft field again as an independent dephasing term.

Reproduce from kinetics/ with OPENBLAS_NUM_THREADS=1 and
early_pole_reversible_momentum_stationary_line_reference.py or
early_pole_reversible_soft_jump_comparison_reference.py (--level 0–2).
Raw controls are in results/.
