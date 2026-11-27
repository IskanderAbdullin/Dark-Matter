# The locked sound branch does not inherit the relative-velocity drag at leading order

The [static longitudinal control](INDEPENDENT-CHARGED-LEPTON-LOCAL-LONGITUDINAL-2026-11-27.md)
shows appreciable changes in a finite-wavelength relative-velocity pulse
when replacing its momentum coefficient by the thermal one. A different
question is the damping of the *common*, long-wavelength acoustic pair.
At full chemical equilibrium η_D=0 and with the two fluids at the same
temperature, the thermal pressure variation along the locally
equilibrated slow state obeys δP_i/h_i = δT/T for each fluid.
Consequently that variation accelerates the two fluids equally: it
does not directly feed the relative velocity.

Whiten the generator with the positive entropy metric. Its evolution is
−C−ikB with C and B symmetric and C positive semidefinite. Split slow
conserved coordinates (common energy and momentum) from fast
coordinates (two scalar relaxations and one relative velocity). If v is
a normalized eigenvector of the streamed slow block, its attenuation is

    Re λ = −k² (B_FS v)ᵀ C_F⁻¹ (B_FS v) + O(k⁴).

The fast projection B_FS v lies entirely in the scalar block for
this local equilibrium. Changing K_p at fixed thermal and chemical
coefficients therefore leaves this k² term unchanged. The leading
dependence on relative drag enters at fourth order: the resolvent
path must go from scalar to relative velocity and back. With g the
scalar fast projection, the difference between two choices of
relative decay rate Γ_p is

    Δ Re λ = k⁴ || B_vector,scalar C_scalar⁻¹ g ||² Δ(1/Γ_p) + O(k⁶).

The rate and wave number are expressed in the same Γ_Q units
as the preceding longitudinal tests. The result depends on the
full-equilibrium and tightly thermalized closures; it does not
assert a drag cancellation in an expanding, chemically displaced,
collisionless or more complete cosmological system.

## Checks

Across all ten thermodynamic controls the k² coefficient is reproduced
by direct small-k spectra and an independent positive entropy-loss
quotient, even when direct subtraction of nearly imaginary eigenvalues
would be poorly conditioned. Changing the relative drag leaves the
locked speed and scalar damping coefficient unchanged to numerical
precision. With h_V/h_D=10, the k² coefficients are
.00105868 at .015 GeV and .0000232382 at .050 GeV.
For the full-minus-heat-proxy real eigenvalue, the projected k⁴
coefficients are +.000105796 and +.00000337104, respectively.
Check both responses directly at three shrinking wave numbers and
two enthalpy ratios.

Replay with the same combined momentum receipt:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  charged_lepton_locked_acoustic_damping_reference.py \
  results/charged-lepton-portal-combined-momentum-2026-11-27.json
```
[Full acoustic receipt](results/charged-lepton-locked-acoustic-damping-2026-11-27.json.gz).
