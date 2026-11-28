# Fixed-mode angular aperture, broadened radial support, and thermal balance

This is an intentionally normalized regulator experiment on a frozen pole model,
not an inferred in-medium width. At x=.05 and scalar q/T=4.1e-6 the on-shell
parent momentum interval is [1.4290163, 1.4426010] T with the illustrative shared
gauge squared mass shift. A parent at the midpoint can vary its energy defect
by only 8.19997e-6 T by changing its direction.

For d = E_P - E_b - E_s, the angular transformation is exact:
p² dp dcos/(E_P E_b) becomes p dp dd/(q E_P). The aperture at any p is the
regulator-CDF difference at the two exact energies. The spin amplitude and
other common prefactors are held constant. Integrating over parent p is
essential: a reduction at a selected shell is not by itself a reduction
in the total weighted seed.

| Regulator width / full angular span | At midpoint, Gaussian / Lorentz | Seed versus narrow value, Gaussian / Lorentz |
| ---: | ---: | ---: |
| 1 | .38292 / .29517 | .999962 / .989303 |
| 5 | .07966 / .06345 | .999025 / .947372 |

Widths denote Gaussian standard deviation or Lorentz half-width, respectively.
The seed integrals above use p/T up to 4. Extending to 12 yields .989495 and
.948331 for the two Lorentz cases; a cut at p/T=1 excludes the pole, but the
Lorentz tails still generate .00287 and .01426 of the uncropped narrow seed.
Gaussian tails there are negligible at these widths. Thus even normalized
widths do not define a correction without specifying the line shape and
the represented momentum domain. 96 and 192 probability-space quadrature
nodes agree within 4e-17 in the three reported raw integrated weights.
The narrow Gaussian tends to the independent on-shell one-dimensional
integral.

**Detailed balance matters before comparing signs.** The uncorrected
fermionic reverse/forward ratio at common T is exp[(E_s+d-Δμ)/T]. Smearing
only the delta therefore leaves nonzero net stimulated exchange at the
thermal point Δμ=E_s. In the broad Lorentz case, truncating at p/T=1 gives
a spurious negative stimulus of about −2.20e-8 in the reported units.
If the energy defect is exchanged with a thermal reservoir, a toy pair
of weights e^(d/2T) for forward and e^(−d/2T) for reverse restores
exp[(E_s−Δμ)/T] pointwise. This gives exact cancellation at the gate and
opposite signs at ±10% of E_s around it. It does not prescribe the
microscopic exchange kernel. In particular the common-T factor does not
apply to the earlier separate-temperature fermion control.

Reproduce from kinetics/ with early_pole_finite_resonance_aperture_reference.py
(--nodes 96 or 192); both raw results are in results/. This does not validate a
cosmological coupling threshold, extend the trajectory in x or define CLASS input.
