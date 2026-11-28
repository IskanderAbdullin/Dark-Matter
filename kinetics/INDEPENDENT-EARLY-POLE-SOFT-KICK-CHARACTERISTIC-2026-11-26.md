# Soft kicks, momentum diffusion, and the same neutral phase

We independently integrated the classical longitudinal and transverse
hard-thermal-loop spectral weights at omega=k_parallel for an eikonal
test charge. With transverse cutoff c in Debye units, the single-line
second moments divided by g² T m_D²/(4π) should be
ln(1+c²) transversely and .5 ln(1+2c²) longitudinally. These are the
cutoff-integrated leading expressions of
[Ghiglieri, Moore and Teaney, Eqs. 22–29](https://arxiv.org/html/1509.07773).
We use only their leading Abelianized soft input, with charge squared
one and equilibrium two-Dirac m_D²; no QCD NLO term enters. The direct
two-dimensional sum agrees with both expressions through the five
cutoffs 2–40 to better than 4e−6 relative.

**The neutral phase and the soft radial kick are not separate inputs.**
For a common test-particle kick, differentiate the resonance detuning
with respect to the hard radial momentum: its derivative is the parent
radial velocity minus the daughter's projection. That same velocity
difference separates the two prescribed paths. For symmetric
independent soft kicks with a linearized detuning, a kick at time s
builds phase alpha q_parallel(t−s). Averaging over kick time gives
a log-amplitude t times the kick integral of
1−sinc(alpha q_parallel t). This is precisely the previously evaluated
local Doppler phase integral when the cut, common velocity and path
projection are matched. The implementations agree within .1% across
the spherical controls. Adding that phase and a radial kick damping
as if independent would double count this local soft effect. This
does not identify the initial, final, or collective-mode terms in a
finite contour with a collision operator.

A Gaussian momentum-diffusion approximation retains only the small
argument term and instead gives alpha² qhat_parallel t³/6. Its
phase exponent is always at least as large in this symmetric-kick
comparison. For the same spherical cutoff k ≤ 20 m_D and g=.005:

| Time in inverse angular half-spans | Soft kick exponent | Diffusion exponent |
| ---: | ---: | ---: |
| 1 | .05214 | .07888 |
| 4 | 1.2156 | 5.0486 |

At one inverse half-span, the radial RMS change is .688 times the
radial distance across half the angular window, even though roughly
76% of its second moment comes from transfers that individually
span more than that distance. The transverse RMS deflection is
.00473 radian for the same probe. At g=.003 the kick/diffusion
exponents at one half-span are .00827/.01022. The cutoff on transverse
momentum alone allows arbitrarily large longitudinal transfer in the
formal soft sum rule; its exponent is larger, and is not physical
support beyond the soft region. Changing the integration geometry
must be described before comparing numbers.

The drag proxy from Eq. 30 evaluated with the spherical soft moments
is about 1.39e4 H_ref at g=.005 under the 86 GeV radiation clock, but
only 8.2e−6 of the angular half-span in frequency units. This is a
single test-charge comparison, not a common-temperature or chemical
equilibration rate for the evolved nonthermal species. Fast compared
with cosmic expansion does not imply diffusion is valid over the
microscopic resonance width.

The shot-noise test assumes independent symmetric kicks, no energy
drag, no hard matching and no in-medium transition vertex. Its
finite-time characteristic is not a stationary spectral function.
Reproduce from kinetics/ with early_pole_eikonal_soft_kick_sum_rules_reference.py
and early_pole_kick_characteristic_reference.py (--nodes 8 or 12).
The integrated exponents agree between the angular/radial grids
within 3e−9 absolute; the fractions across a sharp threshold are
rougher, up to .016 absolute.
