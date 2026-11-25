# Finite mother chemical relaxation, fixed radiation bath — 2026-11-25

A key limitation of the [fixed-bath FRW model](RADIATION-FIXED-BATH-2026-11-25.md) is its inexhaustible thermal parent. As a next controlled approximation, enforce a finite mother comoving number N1=n1/T³ and instantaneous *kinetic* equilibrium at externally specified T with dimensionless chemical potential mu determined by N1 via the Fermi integral. The same pair-area collision operator uses f1=FD(E_i+E_j-mu). Let lambda=Gamma_chemical/H be a prescribed constant; dN1/dln a= -dNphi/dln a + lambda(Neq-N1). For lambda=0, N1+Nphi and N1+Nchi2 are conserved; for nonzero lambda the reservoir can replenish or remove mothers. This is a phenomenological number-relaxation closure, not a derivation of chemical rates or the full mother distribution. The bath remains external; no total-energy conservation claim.

At x0=.1, kappa=.1347, xend=40, initial parent in equilibrium, final results:

| lambda | 400-cell omega_phi h² | mean q | q>.1 fraction |
|---:|---:|---:|---:|
| 0 | .183782 | 2.16125 | .739785 |
| 10 | .122703 | 1.18962 | .496969 |
| 100 | .120079 | 1.03131 | .437454 |
| fixed mother, 480 cells | .120101 | 1.01115 | .429523 |

lambda=1 on 240 cells gives omega=.16300 and q>.1 fraction=.68921. Finite parent number suppresses early cold production, and the late decays of the remaining parent population can increase the hot fraction/abundance. At lambda=0, the final number equals the initial N1=.364681857 to ~1e-13; it is not a matched-abundance cosmology. The 240→400 cell changes in reported mean and number for lambda=10,100 are <.6%, and hot fractions differ by <.15 percentage points. This supports the sensitivity to chemical replenishment in the chosen model, not a universal no-go result.

The analytic Jacobian is checked against a directional finite difference. Mother number is permitted to undershoot zero numerically by <3e-13 at extremely late times, with negative values treated as empty in the Fermi inversion; no claim of exact positivity. [Executed code](finite_mother_bath.py), [selected numerical outputs and commands](results/finite-mother-bath-2026-11-25.json). No CLASS run for these spectra.

## Checkpoint

Next close the total stress-energy budget: use a comoving temperature reference proportional to 1/a rather than fixing the actual bath T, evolve the bath temperature from total FRW continuity, and calculate H from the evolving radiation/mother/daughters. This retains kinetic equilibrium and the phenomenological chemical source but measures the impact of energy feedback. Ultimately replace those closures with physical collision kernels and finite masses.
