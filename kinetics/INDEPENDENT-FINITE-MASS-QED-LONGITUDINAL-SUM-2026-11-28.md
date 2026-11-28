# Longitudinal Coulomb cut and the Debye static sum

Extend the [frozen longitudinal timelike pole](INDEPENDENT-FINITE-MASS-QED-FROZEN-LONGITUDINAL-POLE-2026-11-27.md) and the [positive free spacelike density source](INDEPENDENT-FINITE-MASS-QED-SPACELIKE-DENSITY-2026-11-28.md) into one collisionless retarded Coulomb check. Use the same two massive isotropic charged species, held endpoint histograms and mirrored antiparticles. Work in a normalization where the inverse Coulomb propagator is k²−Π₀₀, the static Π₀₀ is −m_D², and the positive spectral density at positive frequency is −2 Im D_R.

For each nonzero velocity boundary, with jump weight W = −Δf p E, the dimensionless angular kernel is

```text
u/(2v) log[(u+v)/(u−v+i0)] − 1,  u = ω/k.
```

The real part is verified against independently evaluated angular principal values, to 8.16e−16 in Debye units on both input histograms. Minus its imaginary part agrees with the positive radial density formula, including the density threshold term, to 4.45e−16 relative. The exact free finite-recoil imaginary response stays within 8.45e−6 of the soft value over the sampled k≤.01 and phase speeds.

For each sampled k below the measured timelike pole endpoint, split the cut quadrature at every histogram velocity jump. The real logarithms at the exact boundary do not carry a discrete spectral weight. Within each panel a sine-squared map clusters nodes at both ends, with independent 36- and 128-point rules. The real pole contribution uses its independently differentiated residue. In these units the static identity is

```text
k² integral_0^infinity dω ρ00(ω,k)/(π ω) = m_D²/(k²+m_D²).
```

It follows from the subtracted retarded dispersion relation with high-frequency limit 1/k² and static limit 1/(k²+m_D²), as also used by [DeRocco, Galanis and Lasenby, Appendix D](https://arxiv.org/abs/2201.05167). Their application is different; here this is a numerical check on two frozen portal spectra. The same cut-plus-pole coordinates also reproduce a direct imaginary-frequency angular evaluation at three distinct frequencies, rather than only the static endpoint.

| k | fine cut | third cut | fine cut below held scalar energy | third cut below held scalar energy |
| --- | ---: | ---: | ---: | ---: |
| .0002 | .0038271 | .0038208 | .0038271 | .0038208 |
| .0008 | .0560677 | .0559837 | .0560677 | .0559837 |
| .00128 | .1250707 | .1249092 | .0922339 | .0916437 |
| .00256 | .2826554 | .2825050 | .0851739 | .0846727 |
| .004 | .2948627 | .2949926 | .0493849 | .0491344 |

These are contributions to the same dimensionless **inverse-frequency spectral sum**, not collision probabilities or ordinary area under ρ. Maximum absolute defects are 3.04e−14 for the static sum, 3.01e−14 for the three imaginary-axis checks, and 3.81e−12 under changing the quadrature order.

Replay via `python -W error finite_mass_qed_longitudinal_static_spectral_reference.py` from `kinetics/`. Raw values and retarded residues are in the [JSON receipt](results/finite-mass-qed-longitudinal-static-spectral-2026-11-28.json). This does not assign a statistical state to the pole, use k beyond the documented pole end, infer any physical width, or add the longitudinal cut to the scalar collision term.
