# 09. Continuous family of **exactly critical** lowered Bose equilibria and the infinite-mass limit

**Date:** 2026-11-21. This is a one-parameter continuation of paper 04 with independently integrated finite-mass profiles. The lowered distribution is a phenomenological energy-cutoff ansatz with classical King-model antecedents; mathematical novelty beyond this program not established.

## 1. One parameter controls both the edge and the normalized central phase density

Recall the positive occupation \(f(q;\psi,\eta)=[e^{\eta-\psi+q^2}-1]^{-1}-[e^\eta-1]^{-1}\), for \(0\le q<\sqrt\psi\), and zero outside; \(\psi=(E_t-m\Phi)/kT\), \(\eta=(E_t-\mu)/kT\). With central **saturated fugacity** \(\psi_0=\eta\) one obtains the finite integral
\[
I_c(\eta)=\int_0^{\sqrt\eta}q^2
\left[\frac1{e^{q^2}-1}-\frac1{e^\eta-1}\right]dq.
\]
Its exact physical central degeneracy normalization is
\[
\boxed{\nu_{\rm lowered,c}(\eta)=n_0\lambda_{\rm th}^3
=\frac4{\sqrt\pi}I_c(\eta),}
\]
which is **strictly smaller** than the untruncated critical value \(\nu_c=\zeta(3/2)\) for every finite \(\eta\). Thus **critical central occupation \(z_0=1\) does not imply \(\nu=\zeta(3/2)\) after energy truncation**. This makes blindly transferring the original source's \(\nu_c\) boundary or 86-eV mass into the new model incorrect.

With density \(x(\psi)=I(\psi,\eta)/I_c(\eta)\), integrate \(\psi'=-Q/\xi^2\), \(Q'=\xi^2x(\psi)\), with \(\psi(0)=\eta\), \(Q(0)=0\). The outer radius \(\xi_t\) is first \(\psi=0\), and total mass is \(M_t=4\pi\rho_0R^3Q_t\).

## 2. Central slope and analytic endpoint limits

The ideal zero-momentum divergence is integrable in \(I_c\). As \(\psi=\eta-y\to\eta^-\),
\[
I(\eta-y;\eta)=I_c(\eta)-\frac{\pi}{2}\sqrt{y}+O(y),
\quad
x(\xi)=1-\frac{\pi}{2\sqrt6 I_c(\eta)}\xi+O(\xi^2).
\]
At small \(\eta\), use \(1/(e^x-1)=1/x-1/2+x/12+O(x^3)\):
\[
\boxed{I_c(\eta)=\frac23\sqrt\eta-\frac1{90}\eta^{5/2}+O(\eta^{9/2})}
\]
and hence \(\displaystyle -x'(0^+)\sim\frac{3\pi}{4\sqrt{6\eta}}\). This is a *formal central slope*, not a global galaxy cusp exponent.

At large \(\eta\), lowering vanishes in the interior of phase space and
\[
I_c(\infty)=\int_0^\infty\frac{q^2dq}{e^{q^2}-1}
=\frac{\sqrt\pi}{4}\zeta(3/2),
\quad
\boxed{-x'(0^+)\to
\frac{2\sqrt\pi}{\sqrt6\zeta(3/2)}
=0.553979545797\ldots}.
\]
This exactly recovers paper 1's critical untruncated Bose halo central slope, a strong analytic internal check on our energy-cutoff formulation.

## 3. Independent numerical scan, no observational fitting

For seven values \(0.5\le\eta\le6\), with \(\psi_0=\eta\) and 240-node Gauss–Legendre momentum quadrature and a monotone equation-of-state interpolation followed by Poisson integration:

| \(\eta\) | \(\nu_{\rm lowered,c}\) | Linear center slope | \(r_t/R\) | \(Q_t=M_t/(4\pi\rho_0R^3)\) |
|---:|---:|---:|---:|---:|
| 0.5 | 1.05944 | 1.36600 | 8.4536 | 1.0255 |
| 1 | 1.48011 | 0.97777 | 12.2445 | 2.9158 |
| 2 | 2.00004 | 0.72359 | 19.0712 | 8.4281 |
| 3 | 2.29326 | 0.63107 | 27.5266 | 16.1008 |
| 4 | 2.45368 | 0.58981 | 40.2988 | 26.3570 |
| 5 | 2.53670 | 0.57051 | 61.9187 | 40.3741 |
| 6 | 2.57751 | 0.56147 | 101.7407 | 60.7449 |
| \(\eta\to\infty\) | 2.612375… | 0.5539795… | divergent limit | divergent limit |

As the escape cutoff is loosened, \(r_t\) and \(Q_t\) increase monotonically over the sampled interval, while the central linear density coefficient decreases and the central normalized degeneracy rises toward \(\zeta(3/2)\). The last row is an analytic **limiting expectation** based on recovery of the untruncated infinite-mass isothermal tail, not a numerically integrated \(\eta=\infty\) finite-radius solution.

Our numerical assertions checked monotonic trends and \(\eta=4\) consistency with an independent 512-node critical computation. **The scan does not prove global monotonicity at all \(\eta>0\)**, uniqueness or stability of the entire solution family.

## 4. Why this extension matters

- A finite-mass near-critical **core-plus-envelope** exists as a one-parameter family within an explicitly nonnegative stationary \(f(E)\), not only as one isolated example.
- The original source's \(\nu=\zeta(3/2)\) surface is not the correct central critical constraint after truncation, even though the center's fugacity remains saturated.
- The untruncated ideal linear center coefficient is recovered **continuously** by removing the cutoff; this is stronger than matching a few rotation-curve points.
- Energy cutoff controls the whole density concentration and outer tidal scale; no physical inference should treat this extra parameter as free without a prior or formation history.
- Stationary positive DF + finite mass + virial identity (paper 07) remain insufficient to establish nonlinear Vlasov, thermodynamic, or cosmological stability.

## 5. Next step

For each \((\eta,\psi_0)\) calculate a calibrated sequence of conserved \(M,E,\mathcal C\), radial frequencies and a global linear Vlasov eigenvalue spectrum; seek turning points and stability transitions with the exact dynamically accessible perturbations. At the same time derive the evaporation/tidal mechanism that sets the cutoff instead of arbitrarily prescribing \(\eta\).

References: Abdullin & Popov JCAP 2021 055, https://doi.org/10.1088/1475-7516/2021/11/055; finite-extent models MNRAS https://doi.org/10.1093/mnras/stad117 and https://doi.org/10.1093/mnras/stag050; Guo–Rein King stability https://doi.org/10.1007/s00220-007-0212-8; Chavanis Bose halo review https://doi.org/10.3389/fspas.2025.1538434.


## Appendix A — exact monotonicity theorem for central phase density (beyond the numeric scan)

Let \(B(u)=(e^u-1)^{-1}\). The critical normalizing integral can be written in an exactly differentiated form:
\[
I_c(\eta)=\int_0^{\sqrt\eta}\!q^2 B(q^2)dq-
\frac{\eta^{3/2}}{3} B(\eta).
\]
The two moving-boundary terms cancel under differentiation, leaving
\[
\boxed{\frac{dI_c}{d\eta}=
\frac{\eta^{3/2}}{3}\frac{e^\eta}{(e^\eta-1)^2}>0
\quad\text{for every }\eta>0.}
\]
Thus **central quantum-cell normalization** \(\nu_{\rm lowered,c}=4 I_c/\sqrt\pi\) is *strictly increasing* and the one-sided central density slope \(A_c(\eta)=\pi/(2\sqrt6 I_c)\) is *strictly decreasing* for **all \(\eta>0\)**, not only the seven computed cases. As \(\eta\to\infty\), \(\nu_{\rm lowered,c}\nearrow\zeta(3/2)\) from below. This theorem does **not** establish global monotonicity of the total mass or cutoff radius; those remain numerically observed trends on the tested grid. A finite-difference check of \(dI_c/d\eta\) at \(\eta=0.5,1,2,4,6\) agreed with the exact derivative at relative precision \(\lesssim1.2\times10^{-5}\) using the existing quadrature implementation.
