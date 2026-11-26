# 20 — Signed C22 source and whole-halo first Born response
2026-11-25. Follows paper 19. Critical lowered Bose halos, eta=4, same physical parameters. First Born/Duhamel approximation with fixed equilibrium potential and initial collision source. No evolving collision operator, no self-consistent response potential, no claim of equilibrium or real halo evolution.

## Conservative signed source
Use paper 17's reversible center-of-mass measure on A∈(0,psi), theta∈(0,pi/4) with reflected angular measure; incoming u≤u_c and u≤v≤1. Split v at u_c and apply logarithmic maps to A,u,v near the Bose pole. Every weighted event deposits (-1,-1,+1,+1) into binding bins s=psi−q², with paired gain-minus-loss factor P_u(1+S_v)−P_v(1+S_u). Alongside count C_i deposit binding moment M_i using actual s of each leg. The sum of counts and s moments is zero event-by-event, including particles above imposed cutoff and physical escape; local kinetic-energy conservation follows. Hard bin boundaries are -Qt/rt and 0, and near the deepest state. Numerical count/energy ledgers do not imply fine-bin convergence.

n=12/20/28, nr=48 mapped:
|model|n|cutoff /Gyr|true escape /Gyr|deep gap<.2 /Gyr|
|---|---:|---:|---:|---:|
|linear|12|.002670043|.000515550|.000161251|
|linear|20|.002662328|.000516142|.000162734|
|linear|28|.002659313|.000514286|.000163799|
|quadratic|12|.000613945|.000348215|.000254188|
|quadratic|20|.000621620|.000353758|.000254186|
|quadratic|28|.000624643|.000349527|.000254527|

Paper 16 targeted threshold references .002663467/.000516493 and .000625507/.000350195. Paper 17 deep gap .2 reference .001383669*.11788635 and .001534494*.16690623, respectively. Hard-bin discontinuities cause nonmonotonic percent/subpercent errors: no high-precision claim. Global count and binding errors <3.2e−17 per Gyr; largest local ledger relative to integrated absolute event weight <1.73e−13. The first batched n28 quadratic run timed out, subsequently completed individually.

## Actual signed tracer streaming
For each radius and bin represent C,M by two signed weights at bin endpoints a,b, clipped to physical s≤psi: W_a=(bC−M)/(b−a), W_b=C−W_a. This preserves local count and energy exactly. Sample isotropic mu∈[-1,1] by Gauss quadrature, v/V0=sqrt(2(psi−s)), then integrate planar Cartesian orbits with velocity Verlet in initial spherical potential (includes true-unbound trajectories). The entire signed source is propagated, not just ejecta. With unperturbed Liouville propagator U0(t), delta f_Born(T)=integral_0^T U0(t) C[f0] dt. Constant source approximation, no additional sink, all weights signed. At every t sum all tracer count, and monitor signed orbital-energy error. Accumulate smooth Gaussian radial mass probes and sharp enclosed masses via trapezoid in time. In particular delta mass at t=0 is zero because all local counts cancel.

Two local resolutions: coarse n12/nr32/na6/dt=.002 Gyr, fine n16/nr48/na10/dt=.001 Gyr. Fine 34900/29140 signed tracer trajectories. At T=1 Gyr, fractional initial total mass:
|model|Gaussian l=.5R|l=R|l=2R|l=5R|l=10R|delta M(<rt)|
|---|---:|---:|---:|---:|---:|---:|
|linear|.000182147|.000628926|.001331889|.000669087|−.001239813|−.001386615|
|quadratic|.000213953|.000623517|.001064375|.000335210|−.000788918|−.000001372|

Coarse smooth values linear [.000181924,.000628195,.001331392,.000670425,−.001238148], quadratic [.000214499,.000624435,.001065772,.000334976,−.000792433]. Smooth signs robust, subpercent differences. Sharp enclosed values at r/R=1,2,5,10: linear [.000373888,.001255582,.001829137,−.001113963], quadratic [.000426058,.001171623,.001017844,−.000945899]; percent-level resolution sensitivity. Tiny quadratic outside fraction has appreciable relative error; no precision claim. Linear mass beyond 2rt at 1 Gyr .000133191. Unlike first-excursion tracers in paper 19, here source includes negative weights and all energies. Global mass ledger ≤7.5e−19; initial signed orbital energy ≤2.3e−18, fine final weighted energy drift 2.50e−9/−2.49e−8, mean absolute per-tracer dimensionless energy drift .00012215/.00010480, about 4× improved relative to coarse. This is numerical integration error, not physical energy nonconservation.

**Important missing term:** gravitational response is of the same formal first order in sigma as the signed perturbation, not higher order. A small fraction of total mass does not establish that the local response potential is negligible. Next solve the linearized Vlasov–Poisson Duhamel/Volterra response (spherical), including delta Phi acting on f0; then address C'[f0]delta f and evolving field. No nonlinear stability claim, actual lifetime, or condensate formation.

Code: continuation-04/code/lowered_bose_signed_source.py, lowered_bose_signed_response.py. Reproduce: python continuation-04/code/lowered_bose_signed_source.py --model linear --n 28 --nr 48; python continuation-04/code/lowered_bose_signed_response.py --model linear --n 16 --nr 48 --na 10 --dt .001 --T 1; repeat quadratic. Local calculations; CI unverified.
