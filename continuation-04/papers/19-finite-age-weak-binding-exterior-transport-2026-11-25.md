# 19 — Finite-age weak-binding exterior transport
2026-11-25. Independent halo direction, following papers 16–18. Frozen initial spherical potential, prescribed initial C22 birth rate held constant, no subsequent scattering or feedback. This is a first-excursion **tracer** calculation, not total halo evolution or a prediction of an observed halo. Benchmarks eta=4, critical linear/quadratic, R=.6 kpc, V0=44.1 km/s, sigma/m=1.25 cm²/g.

## Differential source and exact exterior geometry
Let e0=Qt/rt, physical binding ep=−E/V0² and e=e0−ep. For 0<ep<e0/2, P=psi+e. Differentiate the cumulative integral of paper 16 at fixed incoming speeds and relative angle. Its upper outgoing integration endpoint is independent of P, lower endpoint psi+P−S; boundary motion gives zero because the outgoing interval has zero length on the angular/kinematic boundary. Thus the positive differential birth integrand is the existing incoming weight times [1+F(psi+P−S)]/B, integrated over the same allowed domain. S(ep) is the resulting rate per initial halo particle per Gyr per dimensionless binding. Cumulative finite differences and integral of S give independent checks.

Physical potential depth U=psi+e0 inside and Qt/r outside. Nonnegative spherical density implies U(r)≤Qt/r and U(r)≥Q(r)/r. At birth L²/2 ≤ r²(U−ep) ≤ Qt r−ep r² ≤ Qt rt−ep rt² for ep≤e0/2. The effective potential has a single minimum because Q(r)r is nondecreasing; the orbit can reach rt for every direction (apart from degenerate measure-zero cases). This is not true automatically for all less weakly bound energies. Explicit isotropic angular integration agrees (fraction 1).

For a=Qt/(2ep), x=rt/a≤1, Kepler period T_K=2pi Qt/(2ep)^(3/2) in R/V0 units. At outward crossing eccentric anomaly E0=acos((1−x)/ecc); exterior residence time T_ext=T_K[1−(E0−ecc sin E0)/pi]. Maximize the mean anomaly over allowed ecc∈[1−x,1]: its maximum occurs at ecc²=1−x (differentiate), giving Mmax=asin(sqrt(x))−sqrt(x(1−x)). Hence
T_K(1−Mmax/pi) ≤ T_ext ≤ T_K.
For the full interior trajectory y=r², y''=4(U−ep)−2Q(r)/r ≥ 2(e0−2ep). Convexity and 0≤y≤rt² bound any continuous interior residence before first outward crossing by Dmax=2rt/sqrt(e0−2ep), also in R/V0 units; proof: midpoint ≤rt²−a_min dt²/8 and midpoint≥0. No instantaneous crossing assumption required.

At 13.8 Gyr, binding threshold below which even the minimal first exterior excursion exceeds the age: ep=.14022191 (linear), .10012458 (quadratic). Initial gross rate into these weakly bound energies .0002277143 and .0001174630 per Gyr, respectively 44.09% and 33.54% of true-unbound creation rate. S(0)=.00136193/.00101911 for n24/nr96; finite differences agree within .000124 relative. Conservative lower bound for constant interior births known to have crossed and not returned at age T is integral S(ep) max(0,T−Dmax(ep)) over T_ext,min≥T, ep<e0/2. At T=13.8: .00279606 and .000298869 of initial halo; not a bound for the actual time-varying C22 source.

## Finite-age residence and divergence of a stationary source
If entry at rt itself is constant, exterior first-excursion population is integral S(ep) min(T,T_ext) dep; upper and lower bounds follow directly from the period bounds. For S(0)>0, mean residence diverges: integral S T_ext dep ∝ integral ep^(−3/2) dep. As T→∞ the prescribed first-excursion population grows as 3 S(0) C^(2/3) T^(1/3), C=2pi Qt (R/V0)/(2^(3/2)), since threshold ep_c=(C/T)^(2/3). A finite-mass stationary exterior cannot be inferred from a constant source with no losses. In the hypothetical stationary outbound branch the large-r number per radius ∝∫_0^(Qt/r) S(0)/sqrt(2(Qt/r−ep)) dep ∝ r^(−1/2), hence density ∝ r^(−5/2); the same divergence. Actual finite-time front, tides, sinks, and rescattering matter.

## Direct trajectory kernel
At each birth radius r and ep sample isotropic velocity direction (mu=|cos(theta)| and both radial signs). Let y=r², L²=2y(U−ep)(1−mu²). Locate inner turning point y_p by bisection of 2y(U−ep)−L². Integrate dy/[2 sqrt(2y(U−ep)−L²)] using y=y_p+(y_hi−y_p)z², interpolate potential from the self-consistent initial Poisson solution. First outward arrival delay d=t(p,rt)∓t(p,r). Exact Kepler dwell as above. For constant births in [0,T], tracer population still awaiting first outward crossing, outside on first excursion, and first-returned is weighted by min(T,d), clip(T−d,0,T_ext), and max(T−d−T_ext,0), respectively; sum exactly T. Multiplication by −ep gives the same tracer orbital-energy ledger. These are not full pair-collision or gravitational-energy ledgers: the other collision partners are not evolved here.
For any radius Rext≥rt, use Kepler anomaly at Rext if apocenter reaches it. The time interval outside Rext is [d+t_R,d+T_ext−t_R]. Its length overlapped with [0,T] gives the cumulative first-excursion tracer mass. Finite front bound follows from dr/dt≤sqrt(2Qt/r): r≤[rt^(3/2)+1.5 sqrt(2Qt) T/(R/V0)]^(2/3).

### Refined local numerical results
n=24, nr=80, ne=32, angular=16, orbit=40, potential points=2400. Fractions normalized by initial halo count (all same particle mass). Only selected ep∈(0,e0/2), not all above-cutoff or true-unbound births.

| model | T/Gyr | born | awaiting exit | outside first excursion | first returned |
|---|---:|---:|---:|---:|---:|
|linear|1|.0006789653|.0002566847|.0004222806|0|
|linear|5|.0033948264|.0002566847|.0030410909|.0000970509|
|linear|10|.0067896528|.0002566847|.0050111607|.0015218075|
|linear|13.8|.0093697209|.0002566847|.0060253046|.0030877316|
|quadratic|1|.0001187741|.0001187741|0|0|
|quadratic|5|.0005938703|.0002437017|.0003501686|0|
|quadratic|10|.0011877407|.0002437017|.0009440390|0|
|quadratic|13.8|.0016390822|.0002437017|.0013953804|0|

At 13.8, cumulative fractions at Rext/rt=1,1.5,2,3,5,10,20: linear [.006025305,.005623850,.004988811,.003530789,.001862326,.000156241,0]; quadratic [.001395380,.001165704,.000826016,.000148717,0,0,0]. rt=24.1793/78.1069 kpc; absolute parabolic front upper 306.632/327.879 kpc. This is the first-excursion tracer distribution, not a full outer density profile. In spherical symmetry exterior matter alone adds no force inside rt, but shifts potential relative to infinity; magnitude from this subset alone ≤ outside fraction times e0 in V0² units, while core mass redistribution is uncomputed.

Number and tracer-energy ledger errors ≤7e−16; independent cumulative source closure 8.24e−5/1.235e−4 at refined grid. Coarse n16/nr48/ne20/na10/ng24/npot1500 outside at 13.8=.006025746/.001395768 (relative difference .0073%/.0278%). Tiny far-edge bins have higher relative errors. Independent Cartesian solve_ivp trajectory control at r=.05rt, ep/e0=.1/.45, mu=.2/.8, both signs, both profiles (16 trajectories) agrees for arrival delay to <5.36e−10 relative. Each trajectory passes analytic delay and dwell bounds; cumulative radial mass monotone and zero beyond analytic front. A higher n32/nr96/ne18 residence run completed for linear (source closure 3.49e−5); the following quadratic run timed out, so no refined quadratic result at that resolution claimed. All runs local, remote CI not verified.

Code: continuation-04/code/lowered_bose_exterior_source.py, lowered_bose_exterior_residence.py, lowered_bose_exterior_delay.py, lowered_bose_exterior_transport.py. Reproduce: python continuation-04/code/lowered_bose_exterior_transport.py --model linear --n 24 --nr 80 --ne 32 --na 16 --ng 40 --npot 2400 (repeat quadratic). The independent Cartesian control uses acceleration -Q(r) x/r³ from the original Halo.sol, not the interpolated potential. No CLASS readiness, actual mass loss, final equilibrium, or condensation claimed.

Next: conservative signed initial C22 energy-bin source (all four legs, including deep/core sinks and true-unbound states); then couple finite-age exterior boundary transport to the evolving interior DF and Poisson field, including return and subsequent scattering. Do not replace returned tracers by an artificial absorbing sink.
