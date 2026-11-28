# Instantaneous affine shear of the same frozen scalar states

Complement the [uniform-pressure test](INDEPENDENT-EARLY-POLE-DAUGHTER-SCALAR-ADIABATIC-2026-11-27.md)
with an independent deformation. Starting from isotropy, scale p_z by
exp(η) and p_x,p_y by exp(−η/2). The momentum-space determinant is one:
occupations and density are transported without resetting, and the
bare mass and scalar coupling stay fixed. At η=0, the angular average
of the quadrupole P₂(cos θ) vanishes, so the derivatives of I, the
mass gap, and the mean pressure vanish to first order.
The Hartree counterterm enters both principal pressures equally.

At η=0 use angle averages ⟨μ²⟩=1/3 and ⟨μ²P₂(μ)⟩=2/15.
If P_perp is the pressure along **one** transverse principal direction,
the slope of P_parallel − P_perp is

    S = Σ n p²/E · [1 − p²/(5E²)], so (12/5) P_kin ≤ S ≤ 3 P_kin.

The stationary energy derivative for any η equals the instantaneous
principal-pressure difference: differentiate the kinetic term at fixed
m, then use stationarity to eliminate the implicit gap derivative.
Hence ρ″(0)=S. Independent 8- and 16-point angular quadratures
re-solve the gap at each deformation, and fourth-order differences
check both the pressure slope and work curvature. The exact p=0
fixture has zero S despite its nonzero Hartree pressure.

| Frozen state | S / P | Fraction of native S in the retained cold aperture |
|---|---:|---:|
| Zero bare mass | .03817 | .00211–.00212 |
| At last local daughter aperture | .03979–.03983 | .00206–.00208 |
| At last represented hard support | .33616–.33623 | .000665–.000668 |
| Bare threshold P/ρ=.001 | 1.21376 | .000265–.000267 |

S is an infinitesimal slope, so S/P can exceed one without implying
a finite pressure anisotropy larger than P. The local-minus-native
difference is positive in these four controls at the same bare mass,
but only about 9.4e−12 to 1.4e−11 of S at zero bare mass and
4.3e−12 to 6.5e−12 at last represented hard support. The stable
subtracted result separates the retained occupation term from the
cold and spectator gap terms. At zero bare mass the gap offsets
roughly a third of the occupation term. These conditional numbers
describe the inherited spectator population and an instantaneous
prescribed affine map, not a propagating shear modulus, collision
relaxation rate, perturbation closure, or transferred cosmological
anisotropic stress.

## Replay

Use the kinetics input paths and paired held histories from the
[uniform-pressure test](INDEPENDENT-EARLY-POLE-DAUGHTER-SCALAR-ADIABATIC-2026-11-27.md).
Split the three cases, in the same order, from
`results/early-pole-daughter-scalar-adiabatic-2026-11-27.json.gz`
instead of the bare-stress receipt (same JSON shape: one list of
three cases). In the three replay commands replace the executable
with `early_pole_daughter_scalar_affine_shear.py` and pass the
corresponding scalar-adiabatic receipts as the fourth argument.
Preserve the angular mesh argument for the second case.
Output: [combined receipt](results/early-pole-daughter-scalar-affine-shear-2026-11-27.json.gz).
