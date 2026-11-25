# 3D classical scalar fragmentation: mixed and quartic backgrounds (2026-11-25)

## Scope

Independent nonlinear fixed-Minkowski classical real scalar calculation, NOT a physical cosmological relic calculation. Potential V=m²φ²/2+λφ⁴/24; define A=√(λ/12) φ_max/m, u=φ/φ_max and scale z=m√(1+2A²)t and spatial coordinates by the same frequency. Equation u_zz−∇²u+a u+b u³=0, a=1/(1+2A²), b=2A²/(1+2A²). Homogeneous solution cn(z,p), p=A²/(1+2A²). Previous exact Lamé band in dimensionless k² is (3(1+A²)/(1+2A²), lower+3p²/(2√(1−p+p²)+2−p)). Radiation FRW expansion omitted locally; valid only if H Δt≪1, not across thermal epochs. No fermions, quantum initial state or metric perturbations.

Periodic box side 32. At t=0 mean u=1, u_z=0; add random Gaussian Fourier coefficients in this band with real-space rms 10⁻⁴, seed 42. This is an explicitly imposed classical seed, not vacuum fluctuations. Identical discrete initial modes across grids. Symmetric kick–exact linear spectral rotation–kick split, project cubic using strict |n_i|<N/4 to exclude cubic aliasing. Energy E=⟨u_z²+|∇u|²+a u²⟩/2+b⟨u⁴⟩/4. E_mean same expression for spatial mean (no gradients). D=1−E_mean/E includes interaction energy and is not a condensate or particle number fraction. U is k>2 share of nonzero-mode kinetic+gradient spectral energy, not f(p) or an observational hot fraction.

## Independent linear control

For A=2.72085 exact band k²=(1.59490037,1.79647282); discrete shells 1.61923197, 1.65778511, 1.69633826, 1.73489140, 1.77344454. Direct integration of the exact linear mode fundamental matrix over 4K gives Floquet μ_z=.02010029,.02846429,.03059165,.02806826,.01930636 respectively. Starting with actual seeded shell variance and zero velocity, matrix powers predict sigma after n periods. N=48, 128 steps/period ratio sigma_sim/sigma_lin at periods 4,8,12,16,20,24,28 = 1.00180,1.00486,1.00734,1.00931,1.01075,1.01008,.99834. This is independent of the nonlinear solver and verifies physical early amplification; later comparison invalid.

## Actual 64-period nonlinear integrations

| A | N | steps/period | first sampled D>0.1 (z) | final D | final U | max sampled abs(E/E0−1) |
|---:|---:|---:|---:|---:|---:|---:|
|2.72085|48|64|270.5682|.37928|.27233|.00356189|
|2.72085|48|128|270.5682|.38482|.27095|.00086437|
|2.72085|64|128|270.5682|.37542|.30282|.00085252|
|120|48|128|215.0710|.57359|.46487|.00091676|
|120|64|128|215.0710|.55103|.50744|.00092061|

Final z=468.00990 (A=2.72085) and 474.63935 (A=120). N48/128 vs N64/128 D at periods 24,32,40,48,56,64:
- mixed: .00073/.00073, .02024/.02022, .27533/.27413, .31137/.30343, .29511/.29895, .38482/.37542.
- quartic: .01300/.01300, .33799/.33716, .28722/.28581, .39652/.39131, .45462/.44814, .57359/.55103.

Nonlinear depletion of homogeneous component and onset are robust within these resolution/step tests. The final UV proxy differs by 3.19 and 4.26 percentage points between grids, so no converged particle spectrum or thermalization result. Energy diagnostic stroboscopic, not maximum substep error. Further finite-volume, seed amplitude, larger cutoff and late-time diagnostics are necessary. Large classical fields do not establish physical BEC, nor eventual Ωh². No physical matching to fermions or expanding quantum kinetic system.

## Reproduction

From kinetics: `OPENBLAS_NUM_THREADS=1 python phi_bose_classical_lattice.py --N 48 --spp 128 --A 2.72085 --periods 64 --save mixed48.npz`, repeat with N=64 and A=120. Checkpoints resume their stored configuration. A=2.72085 and A=120 are rounded previous coherent λ examples at T=20 MeV. Source commit 5eaab2d. Local runs complete; remote CI not checked.
