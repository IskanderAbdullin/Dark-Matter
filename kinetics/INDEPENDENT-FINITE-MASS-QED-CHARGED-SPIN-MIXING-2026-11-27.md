# Separate the charged spin-projector change from energy and residue

The [dispersive bridge](INDEPENDENT-FINITE-MASS-QED-DISPERSIVE-CUT-BRIDGE-2026-11-27.md)
checks only the on-shell energy projection of the charged self-energy.
Here retain the full rotationally symmetric Dirac Hamiltonian matrix
of the two bare charged cuts, at the central Born event on each frozen
archive. The charged particle and antiparticle are populated with the
same tabulated function. There is no photon bath in this comparison.

After multiplication by gamma^0 and azimuthal integration, expand
the self-energy matrix on I, alpha·p_hat and beta. With signed internal
energy t and polar propagator K²=A+B cos(theta), the components are

    -2 t J1,       2 k Jmu,       4 m J1,
    J1 = integral du/(A+B u),
    Jmu = integral u du/(A+B u) = (2-A J1)/B.

A small B/A series avoids cancellation in Jmu. Direct gamma-matrix
angular integration checks all three components separately before
summing either energy branch. Radial refinement resolves the lines
through the archive steps and the integrable coincidence log.

Let P_±=(I±H0/E)/2 be the energy-normalized vacuum Hamiltonian
projectors. For an on-shell Hamiltonian perturbation W their
positive-energy spin-projector derivative is

    dP_+ = (P_- W P_+ + P_+ W P_-)/(2E).

Thus the energy projection and spinor change are different
operations on W. Independently diagonalize H0+hW at positive
and negative h, and replay the scalar-current trace
Tr(P_{+,daughter} beta P_{+,parent} beta). This trace uses fixed
vacuum-event momenta and energy-normalized projectors. It has no
pole-derivative residue or shifted two-body aperture.

The two negative-frequency branch contributions to the relative
spin trace are about +6.30/+6.28 (fine/third), the positive
contributions about −3.94/−3.95, in units g²/(4 pi²). The
combined fixed-geometry change is +1.49e−6/+1.48e−6 at g=.005.
The held-matrix finite diagonalization reproduces the linear
change at that coupling to better than 7e−12 in the relative trace.
This is an isolated contribution, not the net thermal correction
to scalar production. The sign of the full vertex, wavefunction
residue correction, occupation reweighting and plasma treatment
cannot be inferred from the spin trace alone.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_charged_spin_mixing_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-charged-spin-mixing-2026-11-27.json).
