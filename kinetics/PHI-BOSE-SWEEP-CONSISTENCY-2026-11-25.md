# Linear sweep consistency and post-production redistribution
2026-11-25. Companion to [narrow sweep](PHI-BOSE-NARROW-SWEEP-2026-11-25.md); same conditional homogeneous branch and pure quartic. Not a nonlinear prediction.

## Formal failure of the fixed background
Given vacuum initial modes, physical hard-mode flux k³H/(2π²) times N=expm1(2I), and four cold particles consumed per two hard particles, the instantaneous formal fractional cold depletion per Hubble time is Q=√3 λ expm1(2I)/(2π² A²). Q is dimensionless. It is not the actual fraction lost in a Hubble time: I changes rapidly (I∝a^-10 in constant-g radiation). Integrate D(Ti)=∫ Q dln a from Ti to T=100 eV along the unmodified normalized branch. D≥1 signals that neglect of backreaction cannot be consistent with the formal linear prediction; it does NOT determine final depletion or when an earlier initial coherent state fragmented. For the two cases:
| λ | criterion | T (keV) | I | q (eV) | Q | D |
|---|---|---:|---:|---:|---:|---:|
|4.698884e−13|Q=1|404.643|11.0885|6.38735e−8|1|.00479464|
|4.698884e−13|D=1|413.349|13.9129|6.24437e−8|265.334|1|
|3.82587e−8|Q=1|5.43491|4.16921|6.43702e−6|1|.0148046|
|3.82587e−8|D=1|5.69247|6.62424|6.14578e−6|118.084|1|

The very large earlier linear exponents are not predictions of such large particle densities. A cold homogeneous branch cannot be extrapolated through them unchanged. The small-λ example's crossing modes are below the earlier chosen rms reference; the large-λ example's are above it. Within this model q at the Q=1 criterion equals the reference 1.14255e−6 eV for λ≈5.73080e−10, at I=1 for λ≈3.55379e−10. These are not LSS limits or a final fraction.

## Distinguish passage through the band from subsequent scattering
At small A the residence time Δt≈A⁴/(8H). For an initially hard particle p=√3m, the unoccupied-final vacuum cross section on cold particles yields τ_band≈Γ_el Δt. At the D=1 criteria τ_band≈4.05e−9 and .001839, respectively. This vacuum estimate is small even in the large-λ case. It does not rule out Bose-enhanced scattering or coherent wave effects.

After the passage, p∝a^-1, E=m sqrt(1+u²), u=p/m, σv=λ² u/[64πm² sqrt(1+u²)(1+sqrt(1+u²))]. With cold density ∝a^-3, constant-g RD and no prior scatter, the future integrated optical depth has a closed-form check τ=(Γ_el/H)_initial × 2 ln(3/2). Proof: substitute u=√3 a_i/a in ∫Γdt and integrate u/[sqrt(1+u²)(1+sqrt(1+u²))]. Numerical quadrature with massive electron entropy and finite endpoint T=100 eV yields τ≈.72056 and 4.58492e7 at the D=1 criteria. Thus the large-λ free-propagating hard-tail picture is inconsistent with even the vacuum cold-target estimate, while the small-λ branch is not guaranteed collisionless either. This is a diagnostic on a fictitious unchanged cold reservoir, not a many-body relaxation computation. The large optical depth is not the actual number of collisions once the cold state is depleted.

Reproduce python kinetics/phi_bose_sweep_consistency.py. Six crossing/normalization checks from companion and independent local quadrature for both criteria, no remote CI claim. Finite electron-decoupling toy, no initial seed distribution, no nonlinear backreaction, no Bose final-state factors, no fermion production, no transfer calculation. CLASS-ready=false.
