# A boundary exponent selection rule for open Bose energy-space transport

**21 November 2026.** The result is an exact near-boundary asymptotic of a **specified one-dimensional toy Fokker–Planck model**, not a derivation of real self-gravitating quantum Boltzmann collision coefficients. The analysis connects continuation-05's energy-cutoff hierarchy with continuation-06's finite escape density of single-particle states.

## 1. An explicitly open, flux-carrying equation

Let \(s=(E_t-E)/(k_BT)>0\) be dimensionless binding energy. Assume a one-dimensional steady **energy-space** transport with an absorbing boundary \(f(0)=0\), outward/escaping current \(J>0\), effective diffusion \(D(s)>0\) and a Bose-stimulated friction/drift factor:

\[
\boxed{J=D(s)\left[\frac{df}{ds}-\beta f(1+f)\right],\qquad \beta\ge0.}
\tag{1}
\]

The constant-\(D\) case was studied in continuation-05 and yielded a linear edge. Here posit a local power law

\[
D(s)=D_0s^\alpha[1+o(1)],\quad D_0>0,
\]

with an asymptotically bounded \(\beta\). These are **model assumptions**; the actual orbit-averaged Bose collision operator, tidal escape dynamics, potential evolution and energy-state geometry may generate a different equation altogether.

## 2. Derive the exponent and the condition for a finite absorbing occupation

Rearranging the flux law,

\[
f'(s)=\frac{J}{D_0}s^{-\alpha}[1+o(1)]+\beta f(1+f).
\]

For a continuous boundary \(f(0)=0\), \(f\to0\), so the drift term is subleading compared with the first term if an integrable positive solution exists. Integrating from the boundary gives for \(\alpha<1\)

\[
\boxed{f(s)=\frac{J}{D_0(1-\alpha)}s^{1-\alpha}[1+o(1)]}.
\tag{2}
\]

Indeed, the ratio of the drift term to \(f'(s)\) is of order \(s\) with fixed \(\beta\), verifying self-consistency. For \(\alpha=1\), the compulsory diffusion integral is logarithmic; for \(\alpha>1\) it has a power divergence. Thus **no nonnegative continuous absorbing solution with finite \(f(0)=0\) and nonzero constant \(J\)** exists under these local conditions when \(\alpha\ge1\).

The energy-cutoff exponent is therefore **selected** by the effective diffusion law:

\[
\boxed{k=1-\alpha\quad (\alpha<1).}\tag{3}
\]

This does NOT mean every phenomenological Bose cutoff \(f_k\) is a solution of the full microscopic collision equation. It means the near-escape powers must match if that particular local flux equation holds.

## 3. Consequences for the positive Bose Taylor cutoffs

The earlier positive family has \(f_k(E)\propto(E_t-E)^k\) at the edge. For a nonzero \(J\), the simple transport (1) requires

| Effective \(D(s)\) near escape | \(\alpha\) | Selected occupation exponent \(k\) |
|---|---:|---:|
| Diffusion diverges \(\propto s^{-2}\) | −2 | 3 |
| Diffusion diverges \(\propto s^{-1}\) | −1 | 2 |
| Diffusion diverges \(\propto s^{-1/2}\) | −0.5 | 1.5 |
| Finite nonzero diffusion \(D(0)=D_0\) | 0 | **1** |
| Diffusion vanishes as \(\sqrt s\) | +0.5 | 0.5 |
| Diffusion vanishes as \(s^{0.8}\) | +0.8 | 0.2 |

Thus within (1) **a C¹ quadratic cutoff \(k=2\) with nonzero escape flux would demand \(D(s)\sim D_0/s\)**. A finite regular diffusion coefficient instead selects the C⁰ linear edge \(k=1\). Conversely diffusion that vanishes approaching escape selects less-regular edges, not smoother ones. This exposes a genuine tension between choosing an arbitrarily smooth cutoff and maintaining a finite outward particle flux.

For the finite-radius lowered-Bose gravitational equilibria in continuation-06 paper 01, the phase-space density of states \(g(E_t)\) is finite and positive. Therefore the energy population \(N(E)=g(E)f(E)\) carries the *same leading power* \(k\) as \(f\). There is no automatic change \(k\mapsto k+3/2\): that +3/2 is the **spatial density** velocity-integration factor. Still, \(g(E_t)\) is an ingredient of the true orbit-averaged diffusion coefficient, which was NOT computed here.

## 4. Direct numerical verification beyond the formal integral

Solve (1) with \(J/D_0=.01\), \(\beta=1\), starting at \(s=10^{-8}\) with the analytic positive boundary series and integrating the nonlinear Bose drift term. Logarithmic slopes fitted over \(10^{-6}\le s\le10^{-3}\):

| \(\alpha\) | Analytic \(k=1-\alpha\) | Numerical log slope |
|---:|---:|---:|
| −2.0 | 3.0 | 3.000034 |
| −1.0 | 2.0 | 2.000045 |
| −0.5 | 1.5 | 1.500054 |
| 0.0 | 1.0 | 1.000067 |
| +0.5 | 0.5 | 0.500089 |
| +0.8 | 0.2 | 0.200113 |

All integrations remain positive and monotone on the sampled interval. The ratio to the exact leading boundary law tends to unity; at \(s=10^{-6}\) the deviation stays below \(10^{-6}\) for the six examples. These confirm the exponent-selection formula within the **stated ODE**, not astrophysical flux rates.

## 5. What the result does and does not say

- It gives a **necessary local matching rule** for a simple steady diffusion+drift absorbing-boundary model with prescribed positive flux, and highlights how singular diffusion or nonlocal transport might be required for an exactly smoother edge.
- It does not determine \(D_0,\alpha,\beta,J\) for 86-eV dark matter. No boson self-interaction Lagrangian, collision integral, orbit averaging, tidal potential time history, Bose-star transfer or evaporation-rate prediction has been provided.
- With \(J=0\), or a sink/source term concentrated near the boundary, or different global boundary conditions, the conclusion can change. A static hard-cutoff Vlasov halo may exist without being an open kinetic steady state.
- The results do **not** imply observational preference for the linearly lowered halo; the latter still requires stability, galaxy fits and cosmological initial conditions.

**Prior mathematical/physical tools:** King (1966) escape-energy models; orbit-averaged FP transport https://doi.org/10.1088/1674-4527/21/5/108 and https://doi.org/10.1007/s41114-018-0013-8; our positive Bose cutoff sequence continuation-05 notes 15–20; the finite escape-state density calculation continuation-06 paper 01. The general power-law flux asymptotic is a straightforward elementary ODE result; no universal originality is claimed.
