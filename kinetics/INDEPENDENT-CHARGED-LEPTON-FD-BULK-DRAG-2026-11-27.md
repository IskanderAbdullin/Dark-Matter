# Elastic momentum drag cannot be substituted by heat exchange

Use the same finite-vector e/μ/τ elastic process as the
[nonlinear two-bath FD reference](INDEPENDENT-NONLINEAR-FD-PORTAL-REFERENCE-2026-11-26.md),
now at common temperature with an infinitesimal relative bulk velocity.
The stable-pair fugacity is fixed in each control; elastic exchange
preserves dark fermion number. For one event transferring four-momentum
Δp, the linear departure from detailed balance for a small relative
velocity is Δp·Δu/T in spatial components. Isotropy gives
⟨Δp_i Δp_j⟩=δ_ij |Δp|²/3. The scalar thermal-force and vector
bulk-force cross coefficients vanish by parity.

The heat coefficient K_Q multiplies the fractional temperature
difference. The momentum coefficient K_p multiplies relative velocity
in units c=1. They use the same positive equilibrium scattering
weight, but their moment weights are (ΔE)² and |Δp|²/3, respectively.
In a mostly-minus metric |Δp|²=(ΔE)²−t. A separate lab momentum
reconstruction uses the longitudinal CM boost and azimuth-resolved
transverse kick, without forming t in that check. It reproduces the
drag moment, and the thermal component matches the existing
energy-only response to within 3e−13. Altering quadrature order and
thermal cutoff shifts the combined responses by at most 2.0e−6
relative on the nine examples. As a further independent check,
boost the exact incoming and outgoing dark Fermi distributions
along six symmetric directions with the visible bath fixed.
Halving the speed from .001 to .0005 reduces the drag-slope
discrepancy by a factor of four in both tested states; the smaller
probe agrees to 5.4e−7 relative.

| Common T (GeV) | Dark η | K_p / K_Q | τ_Q / τ_p under shared dark bulk |
|---:|---:|---:|---:|
| .005 | 0 | 3.4134 | 5.9471 |
| .050 | 0 | 2.6257 | 3.0727 |
| .300 | 0 | 2.4824 | 2.9138 |
| .6882487 | −1 | 2.3602 | 3.9893 |
| .6882487 | 0 | 2.4046 | 2.8227 |
| .6882487 | 2 | 2.6277 | 1.5253 |

The last column makes an additional explicit closure: fermions and
the scalar share a dark bulk velocity and temperature, the visible
bath is held fixed, and no annihilation or expansion is included.
It compares K_p/(ρ+P) with K_Q/(T C_V|n), not two coefficients
with identical thermodynamic susceptibilities. The table is not a
claim about actual relative velocity, scalar free streaming,
cosmological kinetic decoupling or a linear transfer function.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error charged_lepton_elastic_fd_bulk_drag_reference.py
```
[Full result](results/charged-lepton-elastic-fd-bulk-drag-2026-11-27.json).
