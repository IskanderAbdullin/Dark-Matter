# Expansion-era Bose production: hot-tail bottleneck — 2026-11-25

Building on the [direct daughter-energy solver](DAUGHTER-ENERGY-FV-2026-11-25.md), solve the collision equation in a radiation-dominated background with constant effective degrees of freedom, external FD mother, initially empty daughters, massless production kinematics. This is not yet the closed cosmological model requested in the roadmap.

Use q=p/T, x=M/T, T proportional to 1/a, H=H(M)/x². Then dq/dt=0 between collisions, d/dx=(Gamma/H(M)) x C_x, where C_x is the previous dimensionless kernel with mass ratio x. The exact allowed cell area for qr>=x²/4 is computed by integrating r_hi-c/q on the threshold-crossed part and r_hi-r_lo on the fully allowed part. Symmetric balanced closure gives exact discrete thermal balance and daughter-number conservation. BDF uses ln x.

## Numerical result

For m=86 eV, g*s=106.75 during production, g*s today=3.909, T_CMB=2.7255 K, omega=0.12 corresponds to n_phi/T³=0.2381183. We tune the width parameter, not the final PSD normalization. At x0=0.1, kappa=Gamma/H(M)=0.1347:

| cells | omega_phi h² | mean q | rms q | fraction q>.1 |
|---:|---:|---:|---:|---:|
| 320 | .119786 | 1.01261 | 1.87161 | .430191 |
| 480, stricter time tolerance | .120101 | 1.01115 | 1.87016 | .429523 |

The daughter-energy kernel reproduces the direct quadrature area to normalized max error <1e-11 at x=.1,1,5. Discrete thermal balance and number/Pauli/positivity checks pass. The frozen density changes <1e-6 relative from x=20 to 60. The qualitative burst at x=.2 initially has only 0.13% above q=.1; later spontaneous production/inverse processes change the final distribution to 42.95% hot by that operational threshold. This threshold is not an observational warm-DM bound.

With the specified entropy dilution, a stable, collisionless massless chi2 would add Delta N_eff≈.07035; the final daughter energy (both species treated relativistically at production) is ~3.15% of the assumed bath radiation energy. These are diagnostic quantities, not energy-conservation proofs or CMB constraints. The external mother bath supplies energy and maintains its own FD distribution; feedback, thermal masses, branching fractions and changes in g* are omitted. The stable free-streaming chi2 assumption is optional and must be stated in any subsequent cosmological test.

## Reheating sensitivity (not a matched fit)

At x0=.2, .05, .02, .01, the coarse scans locate near-target abundance examples with hot fractions ~64%,31%,18%,12%, respectively (abundances not identical; see data). Earlier initial time can improve the coldness but drives the IR occupation to ~1e15 and introduces severe physical assumptions about reheating, daughter masses, and interaction rates. This is not a no-go theorem. No unique x0 or microphysical width has been inferred. At fixed x0=.1, weak couplings .001/.01 yielded hot fractions >.998 and underproduction, while .2–1 gave excess number with considerable hot tails.

[Executed code](radiation_fixed_bath.py), [selected local numerical outputs](results/radiation-fixed-bath-2026-11-25.json). Eight additional near-target and reheating parameter scans and multiple coarse exploration runs executed; values retained selectively, code can regenerate full history and cell occupations with --keep-psd. No CLASS calculation yet.

## Checkpoint

Export the frozen spectrum with moment-preserving positive reconstruction in two momentum coordinates for one collisionless 86-eV species. Check the numerical number/energy moments against the raw cells, then attempt linear CLASS with the actual predicted abundance (no arbitrary normalization), including chi2 radiation only if the stable collisionless assumption is chosen. Compare to equal-density CDM and verify the transfer with independently resolved PSDs. A full thermalization/closed energy-budget and observational likelihood remains beyond this benchmark.
