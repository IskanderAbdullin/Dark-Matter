# Separate slow scalar support and fast bulk kinetic modes

These are near-equilibrium finite-support diagnostics and a manufactured nonlinear implicit trial. They extend the joint nonthermal short calculation, but do not determine a full cosmological collision spectral gap, kinetic decoupling temperature or abundance.

## Linearized prompt-return collision map

At common continuous BE/FD equilibrium, the derivative of the prompt-return collision rate with respect to the local shape-exponent perturbation is a positive covariance of the event moment coordinates centered at each parent lab energy. The sign of the generator on actual population perturbations is the opposite: the susceptibility maps shape exponents to minus the corresponding moments. Whiten the matrix with the exact per-cell FD/BE number-energy susceptibilities. The two particle numbers and total energy give three physical null directions; additional tiny numerical eigenvalues are not declared physical zero modes. Validate the positive Gram, physical nulls, and a centered finite difference of the separate nonlinear rate.

For this sample (50 MeV reheating, 45 MeV visible snapshot), the bulk temperature-contrast response is about 8.35–8.93e7 H across the tested meshes by Rayleigh quotient. A direct untruncated thermal heat derivative and separate heat capacities reproduce it to roughly .32% or better. Most of this bulk-direction weight is in fast modes. That finding does not apply to cold scalar perturbations: in elastic-free prompt return, the p/T=.005 probe is unresolved from zero in these projected boxes and p/T=.05 is highly sensitive to cutoff/mesh. No minimum physical rate is inferred from those values.

## Below-pole elastic contribution

A separate joint two-species elastic Gram uses invariant mass, total lab energy, and both lab boost angles, retaining exact energy and occupancy factors. Both constituent cells are represented at both ends. The elastic integral excludes the real pole (s/M^2≤.82 or .72) and is added to prompt return without double counting. Its projector preserves the same three invariants. In the sampled equilibrium case a p/T=.005 scalar-shell perturbation has an added elastic Rayleigh response ~105–107 H, and the p/T=.05 coarse probes ~60 H. This is a specific finite-shell Rayleigh response, not an eigenvalue for the actual nonthermal distribution.

For the frozen linear .005-shell probe over .02 Hubble times, resonance alone leaves the initial overlap essentially one, while adding the below-pole elastic process leaves .119–.123. An independent exact spectral composition of each suboperator (Strang split) agrees within 1.6e-5 in entropy-metric norm on the sampled grids, and improves when the time subdivisions are refined. This is a local controlled linear result and should not be extrapolated to source-depleted or nonthermal backgrounds.

## Manufactured nonlinear implicit step

The backward Euler residual is written in the local BE/FD shape exponents. It advances both species' moment variables with one shared flux. Use the analytic occupancy and parent-fixed-point Jacobian and a line search that retains normal BE support. A centered numerical derivative verifies the assembled non-equilibrium Jacobian to 7e-11 relative. The manufactured strictly populated test reaches .25e-6 Hubble times in one admissible implicit step, with count and energy residual below 1.4e-15 and increased entropy. Refined BE sequences and a separate explicit Heun evolution over a common shorter window differ by .057% and .030% in the recorded L1 moment measure for the corresponding step counts. This checks weighted conservation and stiffness handling; it is not proof that an arbitrary empty tail or a full-expansion physical initial state is admissible at that step size.

## Reproduce

Use prior direct-lepton generated spectra and the same [central visible EOS](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat).
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error kinetics/empty_lepton_pair_spectrum_reference.py EOS_TABLE > RAW_SPECTRA.json
python -W error kinetics/joint_resonant_kinetic_linear_modes_reference.py EOS_TABLE RAW_SPECTRA.json
python -W error kinetics/joint_below_pole_elastic_linear_modes_reference.py EOS_TABLE RAW_SPECTRA.json
python -W error kinetics/joint_resonant_implicit_step_reference.py EOS_TABLE RAW_SPECTRA.json
```

## Manufactured FRW redshift and collision split

On a fixed comoving grid, free redshift transports the existing quadrature occupations in physical momentum and recomputes the shifted kinetic moments without changing comoving number. The compensated energy check independently integrates the kinetic pressure trace `rho−3P` through each half step. A separate implicit collision substep updates both spectra with opposite energy transfer. This isolates expansion, collision entropy, and the small entropy added by reconstructing the limited moment basis.

The manufactured test uses only 2.3e-6 of the actual resonant amplitude to check the split across a finite expansion a=1.0725, starting from a bounded perturbation near equilibrium. Zero collisions close the independent redshift ledger; with collisions, the three- and six-split scalar heat values are -8.81e-4 and -9.01e-4 of initial energy, with positive collision entropy. Number, total energy and equal-opposite collision heat residuals are below 1.3e-14. The reconstruction-only entropy change is separately below 2.1e-10 of initial entropy. This small manufactured coupling and one closed momentum support do **not** establish the convergence or admissibility of the actual physical cosmological update.

Run `python -W error kinetics/joint_resonant_frw_split_reference.py EOS_TABLE RAW_SPECTRA.json` after the same spectrum generation.
