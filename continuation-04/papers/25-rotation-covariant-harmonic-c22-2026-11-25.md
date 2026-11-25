# 25. Rotation-covariant energy-resolved angular C22 and evolving local quadrupole
2026-11-25 UTC. Extension of paper 24, same normalization. No physical spatial halo or condensate claim.

## Exact angular geometry on an energy lattice
On each energy-conserving quartet let S=E1+E2=E3+E4, Δ_i=E_i−E_partner. P magnitude lies in [a,b]=intersection of the two triangle ranges; prior isotropic kernel is uniform in P on this interval. Along P the velocity projection is A_i=(P²+Δ_i)/(2P), transverse squared B_i²=E_i−A_i². Within a pair q_i·q_j=(P²−S)/2. For cross-pair legs the relative azimuth is uniform, so <dot>=A_i A_j and <dot²>=A_i² A_j²+B_i² B_j²/2. For a common uniform global rotation, exact harmonic weak covariances for q^l P_l with l=0,1,2 are respectively 1,dot/3,(3dot²−E_iE_j)/10. Integrate these over uniform P analytically using m2=(a²+ab+b²)/3, m4=(a⁴+a³b+a²b²+ab³+b⁴)/5, m−2=1/(ab), m−4=m2/(a³b³); a=0 implies all Δ_i=0 so inverse terms vanish.
Perturbation δf_i=E_i^(l/2) x_i P_l. Weak matrix K_ji accumulates W sign_j (∂J/∂f_i) covariance, mass diagonal D_i E_i^l/(2l+1). Use paper 23 W with identical energy-pair transitions included at half weight; these vanish for isotropic f but contribute angular scattering. This preserves SO(3) covariance by construction and discrete momentum via Σ_j sign_j q_j=0. Thermal detailed balance: K diag(f(1+f)) symmetric nonpositive (MB f instead of susceptibility). Nonthermal f is not stationary; instantaneous eigenvalues cannot prove stability or slowest continuum mode.

At Emax12, ψ3.9:
|N|MB shear|BE shear (μ=−.3)|linear shear|quadratic shear|
|---:|---:|---:|---:|---:|
|60|.56757868|.75787607|1.11679027|.96328107|
|120|.56590406|.75174705|1.10146310|.94790837|
|180|.56549739|.75011725|1.09654972|.94278102|
|240|.56532862|.74941940|1.09396356|.94037477|
Independent continuous values from paper 24: MB .565685425; linear 1.089020505; quadratic .935518853 (n48). N240 errors −.063%,+.454%,+.519%; no claim Emax12 sufficient for arbitrary angular tails. In contrast cubic velocity lattice coordinate directions have 12–14% error. l0 weak N and E residuals <9e−15, l1 momentum sum residual <4.5e−15. Thermal matrix symmetry relative error <2.5e−16; exact discrete thermal boost residual <3e−17 for l1 at N240. No positive real eigenvalues >1e−9 at sampled N60–240. Nonthermal l2 largest real eigenvalues nevertheless change from −.98828 to −.66373 (linear) and −.87435 to −.56499 (quadratic) N60→240. Thus finite-grid slow modes are not converged despite good pressure moment; no spectral gap claim. Units raw collision normalization; multiply by prefactors in paper 24 for physical Gyr−1.

## Evolving isotropic distribution plus first-order quadrupole
Evolve f0 nonlinearly with symmetric scalar C22 and x by exact l2 linearized K(f0(t)). Initial regular bounded perturbation δf=.01 E f0(E) P2, distinct from shear. Matrix-free angular derivative independently matches assembled K action (N60 initial max difference 8.67e−18). Base SSPRK2 positivity and entropy acceptance; check angular reconstructed minimum as f0+min(E x,−E x/2). This is first-order in anisotropy, does not include O(δf²) backreaction or higher nonlinear harmonics. Time τ dimensionless (τ=pref × Gyr).
|model|N|Emax|dt|initial Q rate|Q(1)/Q(0)|Q(1)/frozen single-rate Q|base outside count(1)|
|---|---:|---:|---:|---:|---:|---:|---:|
|linear|60|12|.02|1.1230633|.32577798|1.00152529|.00494112|
|linear|120|12|.02|1.10939546|.33041641|1.00199588|.00427913|
|linear|120|12|.01|1.10939546|.33039349|1.00192637|.00427923|
|linear|160|16|.01|1.10939546|.33039347|1.00192631|.00427922|
|quadratic|60|12|.02|.96944558|.38025031|1.00252321|.00175776|
|quadratic|120|12|.02|.95512386|.38595305|1.00308906|.00147364|
|quadratic|120|12|.01|.95512386|.38593604|1.00304486|.00147365|
|quadratic|160|16|.01|.95512386|.38593604|1.00304485|.00147365|
N/E conserved near roundoff; entropy increases; tested linearized reconstructed occupation nonnegative. At τ=1 final maximum relative angular amplitude ~.0095. Frozen exponential happens to approximate integrated Q to .2–.3% for this one smooth perturbation and finite-gap local problem, but energy-dependent and nonthermal modes not certified. Initial N60 vs120 changes Q ratio ~1.4–1.5%, base outside counts ~13–16%; do not assert full PSD convergence.

Reproduce: python lowered_bose_harmonic_c22.py --N 240 --Emax 12 --model linear; repeat MB,BE,quadratic. python lowered_bose_harmonic_evolve.py --N 160 --Emax 16 --model linear --dtmax .01. Local NumPy/SciPy only, CI not run. Next: conservative spherical streaming at fixed total binding, exact well-balanced finite-volume Liouville flux, exterior/outflow and angular momentum error before C22 coupling.
