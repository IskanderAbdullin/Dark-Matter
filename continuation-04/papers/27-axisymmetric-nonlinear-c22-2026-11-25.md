# 27. Fully nonlinear axisymmetric UU: exact kernel and conservative angular projection

2026-11-25. Local dimensionless collision operator, not self-consistent halo evolution.

## Exact kernel
Write E=q² and p=q_z. The transverse radius is a=sqrt(E-p²); d³q/(4π)=dE dp dφ/(8π). Integrating two uniform transverse circles with radii a,b gives transverse vector density 1/[π² sqrt((U-P_perp²)(P_perp²-L))], L=(a-b)², U=(a+b)². Convolving incoming and outgoing circles gives I=1/π³ ∫_L^U dx / sqrt((x-l)(x-L)(U-x)(u-x)), where l=min L_i, L=max L_i, U=min U_i, u=max U_i. It is zero if U<=L and is 2 K(m)/(π³ sqrt((U-l)(u-L))) otherwise, with 1-m=(L-l)(u-U)/((U-l)(u-L)). Coincident endpoints are integrable logarithmic singularities, not finite point values. The ordered weak kernel is π I/32 in dE1 dE2 dE3 dp1 dp2 dp3 with E4,p4 determined by conservation; integrated over axial momenta I has marginal 2 min(q_i)/π. Thus it reproduces the earlier isotropic ordered weak kernel min(q_i)/16.

Independent Sobol integrations (2^18 points, seed 0) of the marginal for E=[1,2,1.4,1.6], [.03,3,.9,2.13], [.5,.5,.5,.5], [.1,1.9,.3,1.7] gave relative errors -9.29e-5, +.001043, -6.83e-5, -1.198e-4 respectively. Angular covariances were checked against analytic paper-25 values; near-cancellation relative errors are noisy.

## Direct lattice (negative result)
Exact parabola-cut cell volumes, symmetric nonnegative quartet weights and shared E,p integer sums enforce N,E,p conservation, Bose detailed balance and a discrete H theorem. On Emax=12 the cell volume sum is 13.85640646055=12^1.5/3. Center-only quadrature generated spurious initial Q=.0115; cutcells remove initial Q for isotropic f but not spurious Qdot. Cutcell NE16,K8,m4 gave lowered-linear shear rate 1.146993, Qdot=-.0002284 and exterior derivative .00719849; NE24,K12,m4 gave 1.048100694, -.00031223, .004711795. NE16,K8,m6 seed1 gave 1.01333068, -.000226867, .00676405. Therefore singular-kernel sampling remains uncontrolled. An additional experimental per-energy quartet IPF correction (not the preferred method) reduced active isotropy error .0762 to 5.54e-10 at NE12,K6 but left 27/203 unsupported groups and shear rate 1.33702 rather than continuum 1.08902.

## Alternative: continuous geometry then exact conservative joint rounding
At each E quartet sample P uniformly in its allowed interval, direction of P uniformly on sphere, and independent azimuths of incoming/outgoing pairs. This is the exact conditional angular measure from paper 25. For p/dp=x,y within each pair, first choose shared integer sum T by stochastic rounding of x+y. Conditional on T, round each leg up/down so its expectation is its original p/dp and the two integers sum to T. If f,g are fractional parts and f+g<=1, the one-up branch has first-up probability f/(f+g); if f+g>=1 it has probability (1-g)/(2-f-g). All branches are enumerated with their probabilities. Positive, common quartet weights are reflection-symmetrized. Use angular hat masses D=ΔE/4 integral of hat(p) dp on [-sqrt(E),sqrt(E)]. Exact energy marginal from paper 23 and optional IPF balancing give an isotropic limit at the energy-grid level. The discrete center momentum is conserved; reconstructed physical moment differs in boundary hats. Finite rounding changes angular relaxation; IPF is a numerical correction, not an exact finite-node physical kernel.

Results, lowered-linear psi=3.9,Emax=12:
NE,K,m | shear rate | isotropy L1 | unsupported groups | exterior derivative
16,8,8 |1.26987351|4.56e-10|6/444|.00716851
16,16,9 |1.1658892744|7.4088e-10|0/444|.00716851098
24,16,9 |1.1607420203|8.2306e-10|91/1378|.00497698915
Continuum rate 1.089020505. At NE16,K16 the maximum marginal error was 1.20e-5, mean 8.06e-8; mean weight change .14669. NE24,K16 has unsupported groups mostly inactive in this initial condition; they must not be ignored for evolving distributions. BE residual <=9.3e-17 and N,E,p residual <=7e-17 in the two refined runs. NE16,K8 evolved at tau=.1 with nonnegative f, entropy gain .0027964. No convergence claim. Need finer angular resolution, adaptive support sampling, energy cutoff and time refinement.

Reproduce in continuation-04/code:
python lowered_bose_axis_kernel.py --m 18
python lowered_bose_axis_lattice.py --NE 16 --K 8 --m 4 --tmax 0
python lowered_bose_projected_axis.py --NE 16 --K 16 --m 9 --tmax 0
python lowered_bose_projected_axis.py --NE 24 --K 16 --m 9 --tmax 0

These are local scattering calculations; the next physical requirement is conservative spatial transport on a compatible mesh, then coupled time-dependent Poisson, with exterior boundary and critical core convergence. No physical condensate or cosmological identification follows.
