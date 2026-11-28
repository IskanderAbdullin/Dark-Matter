# A common dressed equilibrium and its collision Jacobian

This is a reversibility test of the existing
[joint packet model](INDEPENDENT-EARLY-POLE-LIVE-JOINT-HARD-PACKETS-FRW-2026-11-26.md),
not a claim that the archived populations are at equilibrium.
At fixed scale factor the two pair channels require a common
chemical potential for photons and both charge-symmetric fermions.
Decay then requires zero scalar chemical potential. The common
chemical potential need not be zero because this truncated model
preserves total photon-plus-charged count.

Fit the photon and fermion populations using the actual variational
energy coefficients of their cells, with the scalar Bose occupancy
computed from its own Hartree energy. The manufactured state is
self-consistent to rounding. With scale snapshots eta=0 and .028,
two different temperatures and two nonzero common chemical potentials,
decay alone, both pair kernels together, and all three channels
together give moment derivatives no larger than 1.3e−15. The
chemical constraint is imposed only for the fixture.

To test the derivative through moving energy remapping, write each
reaction's net flux as its positive logarithmic-mean mobility times
the entropy-force difference. The general gradient perspective is
discussed by [Maas and Mielke](https://arxiv.org/abs/2004.02831).
At equilibrium, s·grad E=s·grad C=0 for each event and
grad S=beta(grad E−mu grad C). Differentiating the moving balance yields
J=L(H_S−beta H_E), where L sums the positive equilibrium rank-one
event mobilities and C is the linear auxiliary count. A derivative
holding the remapped directions fixed would omit the energy Hessian.

For a direction in scaled moment coordinates, compare a direct
finite-difference derivative of the full code to a response with
fixed equilibrium geometry in which only the entropy force is
perturbed by (H_S−beta H_E) times that direction. At the finer
packet resolution and smallest test step, both comparisons agree
within 3e−9 relative. Omitting energy curvature leaves discrepancies
3.56e−4 and 2.28e−3 in the same two cases. Left null tests at this
step leave less than 8e−15 for energy and 1e−16 for count.

An independent full finite-difference matrix on the coarser packet
fixture projects out the two conserved constraints. Its constrained
entropy metric is positive; the similarity-transformed collision
matrix is symmetric within 4.4e−11 relative and all represented
relaxation rates are positive in the test. The smallest rates
are about 3e−5 in the imposed clock and vary with support. This is
not an estimate of a physical equilibration time, a proof of a
fully connected infinite network, or a derivative of the expanding
solution.

Reproduce from kinetics/ with early_pole_dressed_joint_equilibrium_reference.py
(--level 0 or 1). Raw: results/early-pole-dressed-joint-equilibrium-2026-11-26.json.gz.
