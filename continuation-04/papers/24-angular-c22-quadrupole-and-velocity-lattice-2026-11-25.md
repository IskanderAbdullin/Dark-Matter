# 24. Continuous l=2 collision response and conservative angular velocity-lattice gate
2026-11-25 UTC; base d10947c188fc64fcdbb61a6de5846640339e74a0. Local finite-gap halo benchmarks from paper 23. No spatial evolution or BEC claim.

## Continuous pressure response
Use q in paper 23 units, E=q², h=q² P2(qhat·z). Apply a volume-preserving shear q_z→exp(ε)q_z, q_x,y→exp(-ε/2)q_x,y in F(ψ−E). Then δf/δε=−2 E F'(ψ−E)P2; dQ/dε=−∫q⁴f dq by integration by parts (Q=∫d³q/(4π) h f). Compute the exact initial derivative of the weak collision moment, not a relaxation eigenvalue of an evolving nonthermal gas.
For ordered pairs the weak moment is 1/4 ∫q1²dq1 q2²dq2 dc/2 dv/2 R_rel J Δh, with J=f1 f2(1+f3+f4)−f3 f4(1+f1+f2). Set P=2√A cosθ,R_rel=2√A sinθ, E1,2=A±B u,E3,4=A±B v, B=PR_rel/2. Pair exchange and sign symmetry reduce to 0≤u≤v≤1; measure becomes 4 A^(5/2) cos²θ sin³θ dA dθ du dv. Average over arbitrary global rotation and relative azimuth exactly: <P2(a·z)P2(b·z)>=P2(a·b)/5, <(n·n')²>=(uv)²+(1-u²)(1-v²)/2. Hence the contribution of each differentiated leg i is −3 R_rel² F'_i (∂J/∂fi)[<(n'·qi)²>−<(n·qi)²>]/10. Integrate with split logarithmic Gauss maps. For a classical MB f=exp(−q²) independent Gaussian relative-velocity integration gives rate 2^(5/2)/10=0.565685424949238; n24 Amax25 returns 0.5656854249499415 (~1.2e−12 relative). This control is not the Bose thermal fixed point.

Rates in physical Gyr⁻¹ use paper 23 prefactors 1.530710824252 (linear), 1.654075234687 (quadratic); central I_c=1.087257505036,1.006167553208 independently integrated. Selected raw dimensionless values:
|ψ|linear n32|linear n48|quadratic n32|quadratic n48|
|---:|---:|---:|---:|---:|
|1|.00225252|.00225346|.00055779|.00055776|
|3|.158944|.158750|.103634|.103632|
|3.9|1.08648|1.08902|.93500|.935519|
|3.99|1.79523|1.81824|1.65476|1.65486|
|3.9999|2.18253|2.17215|2.05567|2.04782|
For n32 at ψ=0.1,0.5,1,2,3,3.9,3.99,3.9999,3.999999, physical linear rates are 2.60e−6,.00036653,.00344796,.04092253,.24329791,1.66307905,2.74798395,3.34082154,3.42858564; quadratic 7e−8,.00005061,.00092263,.02046440,.17141951,1.54656152,2.73709822,3.40023673,3.47291060. Near criticality convergence is not uniform (n16–32 percent-scale differences), no critical limit claim. Even ψ=3.9 damping time ~0.60–0.65 Gyr, slower than typical deep orbit period ~0.073 Gyr (paper 23), while outer damping is orders slower. A local instantaneous shear response cannot justify isotropic closure or be identified with slowest mode.

## Angularly resolved reversible collision network
For integer velocity nodes inside sphere |v|≤K, group all distinct unordered pairs by exact vector P and sum |v_i|²+|v_j|². All pairs in a group have equal relative speed. Uniform outgoing unordered-pair averaging yields a nonnegative symmetric transition network. Equal pairs have zero relative speed and no transitions. With cell w=h³/(4π), let P_g=<f_i f_j>, S_g=<f_i+f_j>. Each pair deposits wR[(1+f_i+f_j)P_g−f_i f_j(1+S_g)] into both legs. Complexity O(number of pairs), not O(number of quartets). Number, three momentum components and energy conserved identically; Bose BE and MB fixed points, entropy sign by detailed balance. Finite velocity sphere changes the continuum angular measure. Explicit positivity and entropy-checked SSPRK2 gives actual nonlinear local anisotropic evolution, not spatial transport.

Independent orientation control: central difference of ±.01 volume-preserving shears. Compare continuum rates (MB .565685, linear 1.08902, quadratic .935519 at ψ=3.9). qmax4 MB and qmax3 Bose, velocity nodes K3/5/7/9 = 123/515/1419/3071:
|K|model|coordinate-axis|face-diagonal|body-diagonal|
|---|---|---:|---:|---:|
|3|MB|.321453|.503542|.597907|
|5|MB|.446057|.535821|.565748|
|7|MB|.482524|.544913|.565709|
|9|MB|.499113|.549060|.565710|
|3|linear|.409058|1.295745|1.601243|
|5|linear|.842926|1.073936|1.156348|
|7|linear|.914347|1.052424|1.099944|
|9|linear|.949352|1.054511|1.088743|
|3|quadratic|.321482|1.085899|1.405964|
|5|quadratic|.700448|.920800|.996140|
|7|quadratic|.770981|.901363|.945607|
|9|quadratic|.807306|.904181|.936621|
At K9 coordinate rate remains 11.8%,12.8%,13.7% low for MB/linear/quadratic; body direction looks accurate by coincidence of cubic angular symmetry, not sufficient. Averaging over l=2 cubic irreps (2 diagonal,3 off-diagonal) reduces but does not remove error. Initial above-cutoff production at K9 is .00512184/.00133248 dimensionless vs prior paper 23 continuum-scale .00956739/1.53071=.00625029, .00268895/1.65408=.00162565 (finite energy-grid reference, not exact). K7 gives .00441672/.00107159: no monotone convergence guarantee. BE fixed point max |df/dt|≤2.1e−16 on sampled grids; finite initial step conserves N,E,vector P to rounding; K9 t=.1 entropy increases .00280756/.00164328 and Q remains finite; no physical time convergence claim. Gross conservation alone demonstrably does not validate an angular solver. Do not couple this lattice to Poisson for quantitative halo inference until angular quadrature improved.

Reproduce in continuation-04/code: python lowered_bose_quadrupole.py --model MB --psi 25 --n 24; python lowered_bose_quadrupole.py --model linear --psi 3.9 --n 48; python lowered_bose_velocity_lattice.py --K 9 --qmax 3 --model linear --tmax .1. Repeat quadratic, MB --qmax 4. Local Python/NumPy/SciPy calculations; remote CI not run. Next: derive angular harmonic weak operator with exact rotational covariances, independent energy-dependent modes and momentum null ledger, before spatial streaming coupling.
