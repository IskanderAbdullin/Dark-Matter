# 30. Momentum-constrained positive entropy-consistent orbit/local collision flux
2026-11-25 UTC. Harmonic external trap; no claim of physical converged collision operator.

Paper 29's independently sampled overlap does not conserve physical orbit-cell p. Two modifications:
1. Expanded event stoichiometric ν=(-1,-1,+1,+1); X,Y positive forward/reverse products. K=sum w Λ(X,Y) ννᵀ with logarithmic mean Λ. g=K∇S. Per shell project K→K-Kp(pᵀKp)^(-1)pᵀK, retaining PSD, N/E nulls and entropy sign. Evaluate Kp without enumerating orbit quartets by Λ=integral X^uY^(1-u)du and factorized conditional moments. Gauss ng8 passes smooth moving BE and entropy relation; at low occupation logmean quadrature identity error ~1e-9 absolute. No boundary positivity proof for this correction, hence prefer the following approach.
2. For each expanded event, d=ν·p, tilt Jλ=X exp(-λd/2)-Y exp(+λd/2), choose λ separately per shell so Σ w d Jλ=0. The residual is monotonically nonincreasing, derivative -Σ w d²(X exp(-λd/2)+Y exp(+λd/2))/2. For strictly positive f and any nonzero d there is a finite unique root; degenerate d=0 needs no tilt. With zeros a finite root need not exist; code fails explicitly on bracket exhaustion. At finite λ both forward/reverse factors remain positive and loss for a leg vanishes at f_i=0, so quasi-positivity is retained. Entropy production Σ w Jλ log(X/Y)=Σ w Jλ [log(X/Y)-λd]≥0 by the momentum constraint, including limiting interpretation at zeros. Factorize using a_i=f_i exp(λp_i/2), b_i=(1+f_i)exp(-λp_i/2); no extra orbit-quartet expansion. N,E exact because each expanded event conserves them; p is conserved in aggregate per shell, not eventwise. This is a modified numerical operator, not evidence of correct continuum relaxation.

Moving BE control f=BE(.5+e-.1p) has zero constrained source, λ=-.1 in every active shell. Not an equilibrium of full streaming. Uncorrected operator has spurious drag. Across grids from paper 29:
|NE,K,nL,nphase,radsub|uncorrected max abs df moving BE|tilted max abs df|phase perturbation raw p derivative|tilt L1 correction / raw L1|entropy rate raw→tilted|
|---|---:|---:|---:|---:|---:|
|6,4,8,32,1|.0330393|1.12e-15|-.00213930|.28795|.00196122→.00143173|
|6,8,16,64,2|.0167089|2.91e-15|-.000581509|.16654|.00103596→.000873949|
|6,16,32,128,4|.0135509|3.68e-15|-.000150547|.08803|.000737807→.000694996|
|12,8,16,64,2|.0190451|4.23e-15|-.000381753|.11489|.000762841→.000670935|

Onsager projection corrected max ≤3.54e-15 for moving BE and produces similar smooth rates. Exact tilt has no logmean quadrature error; raw λ0 matches independent paper 29 source ≤2.2e-19 count rate. Shell p residual <1.4e-17, N/E <5e-17. Random low occupations down to 4e-12 and ~20% exactly empty cells: min df on zeros positive in all tested grids (≥4.4e-11), finite roots |λ|≤2.66. Reported entropy rate at zeros uses 1e-300 floor solely as a diagnostic; exact one-sided entropy slope can be infinite.

Actual Strang-split nonlinear evolution, same two half streaming remaps for all controls, t=.1,dt=.02:
|grid|weighted L1 corrected vs no collision|weighted L1 corrected vs raw collision|
|---|---:|---:|
|6,8,16,64,2|.0006555053|.0001215426|
|6,16,32,128,4|.0005383590|.0000427724|
|12,8,16,64,2|.0005005912|.0000614417|
Coarse t=.2 also positive and collision entropy nondecreasing. No rejected steps; N/E drift ~roundoff, physical p increment per collision step <1.3e-18. The correction is still material (~8–29% initial source) and not demonstrated to vanish uniformly, especially near criticality or empty states. Do not extrapolate to halo survival. Reproduce: python lowered_bose_orbit_tilt.py --NE 6 --K 16 --nL 32 --nphase 128 --radsub 4 --t .1; analogous lowered_bose_orbit_momentum.py. Remote CI not run.

Next: port the orbit/local coupling to the actual lowered Bose potential, initially bounded-gap and finite negative-binding cutoff; demonstrate density and exterior source against independently integrated initial signed C22; then open/unbound branches and evolving Poisson. The correction is optional and must be tracked against raw rates.
