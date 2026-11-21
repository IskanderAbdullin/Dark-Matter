# 20. A nonzero steady escape flux selects the **linear** energy-edge law in a minimal Bose Fokker–Planck model

**21 November 2026.** This is a conditional one-dimensional *open-system kinetic toy*, not a derived Uehling–Uhlenbeck collision kernel in a spherical galactic potential. The source model is Abdullin–Popov JCAP 2021 and our positive higher-order Bose energy truncations in notes 15–19. The established orbit-averaged Fokker–Planck treatment of tidal escape is relevant prior art, but the coefficients here are intentionally illustrative.

## A. Energy-space equation and an absorbing tidal boundary

Let \(E\) be specific particle energy, \(E_t\) a tidal escape boundary, and \(f(E,t)\ge0\) an occupation. Consider the local energy-space continuity law

\[
\partial_t f=-\partial_E J,
\qquad J=-D(E)\left[\partial_Ef+\beta f(1+f)\right],
\quad D(E)>0,\quad\beta=(k_BT)^{-1}>0.
\]

The \(\beta f(1+f)\) term is a **Bose-stimulated drift chosen for this toy**. An actual orbit-averaged equation conserves particles with energy density of states \(g(E)\): one should propagate \(N(E)=g(E)f(E)\), and the real \(C_{22}\) kernel is nonlocal and not generally expressible as the above simple Fokker–Planck law. This is a conceptual stress test of boundary exponents only.

Choose an absorbing boundary \(f(E_t)=0\) and a stationary **outward** number flux \(J(E_t)=J_*>0\). If \(D(E)\to D_t\in(0,\infty)\) continuously at the edge, then directly

\[
\boxed{f_E'(E_t^-)=-J_*/D_t<0,\qquad
f(E)=\frac{J_*}{D_t}(E_t-E)+o(E_t-E).}
\]

This yields a **linear** cutoff with \(k=1\) in the positive Taylor-remainder family. The C¹ quadratic cutoff \(k=2\) has \(f_E'(E_t^-)=0\), and hence \(J(E_t)=0\) for finite \(D_t\); it cannot carry a **nonzero stationary escape flux** with the toy coefficients. All \(k\ge2\) behave similarly. This makes the arbitrary smoothing order physically consequential.

At zero steady flux \(J\equiv0\), the absorbing ODE \(f_E'=-\beta f(1+f)\), \(f(E_t)=0\) has the unique solution \(f\equiv0\) for finite coefficients. A nontrivial smooth occupied steady state therefore needs an open source/sink, a flux, singular/degenerate coefficients, or a different boundary law. A real tidal steady **bound** distribution usually has escaping, replenishment, and orbital effects beyond this idealized closure.

## B. Exact solvable constant-coefficient example

Set \(s=E_t-E\ge0\), \(D\equiv D_0\) and \(j=J_*/D_0\). A steady solution obeys the Riccati equation

\[
\frac{df}{ds}=\beta f(1+f)+j,\qquad f(0)=0.
\]

For \(0<4j/\beta<1\), define \(d=\sqrt{1-4j/\beta}\), \(r_\pm=(-1\pm d)/2\), \(A_0=r_+/r_-\). Then

\[
\boxed{f(s)=\frac{r_+-A(s)r_-}{1-A(s)},\qquad A(s)=A_0e^{\beta ds}.}
\]

The solution is positive for \(s>0\) up to the first pole

\[
\boxed{s_{\rm pole}=\frac{\ln(r_-/r_+)}{\beta d}.}
\]

For \(\beta=D_0=1\), \(J_*=0.01\): \(s_{\rm pole}=4.6794065505\). An independent high-accuracy Runge–Kutta calculation agrees with this exact formula, e.g. at \(s=.01\): \(f=0.0001005017045102\), at \(s=2\): \(f=0.0663947606276\), at \(s=4\): \(f=1.0258280606\). Importantly, **constant positive escape flux and constant diffusion across arbitrarily deep binding energy do not furnish a globally regular thermal halo**: the simple toy diverges at finite \(s\). It must be used in a finite energy domain or supplied with real energy-dependent transport, a condensate sink, or other physics. This explicitly demonstrates why the slope matching is a **local edge result only**.

## C. Generalized local exponent

Let \(\Delta=E_t-E\to0^+\), and suppose physically derived coefficients scale \(D(\Delta)\sim D_0\Delta^a\), \(J(\Delta)\sim J_0\Delta^b\), with nonnegative exponents and \(f\to0\). If the derivative/diffusive term dominates the stimulated drift, then

\[
\boxed{f(\Delta)\sim\frac{J_0}{D_0(b-a+1)}\Delta^{b-a+1}}
\]

when \(b-a+1>0\) and the signs are compatible. Thus the edge exponent may differ from one if \(D\) vanishes or diverges, if the local flux vanishes, or if energy-dependent sources matter. The condition \(f\to0\) also rules out some exponents; when the drift is comparable one must solve the full balance. This prevents overclaiming that all physical tidal halos must have a King-like linear cutoff.

## D. Connection to our results and to prior studies

- **Note 15/16** showed a positive but overly smooth \(k=3\) cutoff can give \(10^9R\) outer radii and a \(k\ge4\) sequence with numerical noncompact tails. The toy reveals an *independent dynamical reason* why \(k>1\) needs additional physics: an ordinary finite diffusion coefficient with steady nonzero escape selects \(k=1\) locally.
- **Note 17** explained why bounded isotropic monotone \(f(E)\) can be covered by known **collisionless** linear-stability criteria. Those criteria do **not** mean any cutoff is in quantum collisional equilibrium; note 19 proved it fails closed-system detailed balance.
- The model's *steady escape* must eventually be coupled to Poisson, a tidal potential, true quantum or gravitational collisions, and the density of states. Fixing \(D,J\) by hand is not a microphysical prediction.

**References and provenance:** original JCAP 2021 https://doi.org/10.1088/1475-7516/2021/11/055; King and tidal Fokker–Planck evaporation https://doi.org/10.1111/j.1365-2966.2011.20183.x and King (1966) https://adsabs.harvard.edu/abs/1966AJ.....71...64K; quantum kinetic GP+Boltzmann formulation Proukakis–Rigopoulos–Soto https://doi.org/10.1103/PhysRevD.108.083513; established Bose quantum Boltzmann H theorem discussed in note 19. Only the local toy analytical and numeric statements above are actually established here.