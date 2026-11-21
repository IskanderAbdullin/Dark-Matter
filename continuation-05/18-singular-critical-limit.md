# 18. Why the regular stable-looking branch cannot be naively continued into the exact Bose critical center

**New asymptotic statement and numerical check, 21 November 2026.** This extends paper 17 and the previous critical-center investigations. The distinction is functional-analytic, not observational.

## A. Bounded subcritical vs unbounded critical f(E)

Consider the **C¹ quadratically lowered** Bose function \(F_2\) at fixed \(\eta=4\) with central potential depth \(\psi_0=\eta-\delta\), \(\delta>0\). The maximum occupation occurs at \(p=r=0\):

\[
F_{\max}(\delta)=\frac{1}{e^\delta-1}-\frac{1}{e^\eta-1}-\frac{e^\eta}{(e^\eta-1)^2}(\eta-\delta)
=\frac1\delta+O(1),\qquad\delta\downarrow0.
\]

For every \(\delta>0\) it is finite and the entire f(E) is smooth in occupied phase space and C¹ at the escape edge. For \(\delta=0\), \(F_2(r=0,p=0)\to\infty\). A classical Antonov/DFB/Pérez–Aly theorem applicable to the bounded family therefore **cannot be passed to the endpoint in sup norm**. Stability in a weaker norm might still hold, but requires a separate functional-analytic argument.

## B. The central density derivative changes discontinuously

Let \(I_2(\psi;\eta)=\int_0^{\sqrt\psi}q^2 F_2(\psi-q^2)dq\), and \(x=I_2(\psi)/I_2(\psi_0)\). At any \(\delta>0\), regular Poisson potential gives \(\psi(\xi)=\psi_0-\xi^2/6+O(\xi^4)\), so

\[
x(\xi)=1-\kappa(\delta)\xi^2+O(\xi^4),\qquad
\kappa(\delta)=\frac{I_2'(\eta-\delta)}{6I_2(\eta-\delta)}.
\]

For \(\delta\to0\), \(I_2(\eta-\delta)=I_{2,c}-(\pi/2)\sqrt\delta+O(\delta)\), hence \(I_2'\sim\pi/(4\sqrt\delta)\). Therefore

\[
\boxed{\kappa(\delta)\sim\frac{\pi}{24I_{2,c}\sqrt\delta}},\quad I_{2,c}(4)=1.0061675532.
\]

At **exact** \(\delta=0\), the central profile is instead \(x(\xi)=1-[\pi/(2\sqrt6I_{2,c})]\,\xi+O(\xi^2)\), with a nonzero one-sided first derivative. Thus pointwise differentiability at the center is *nonuniform* in the limit. Matching the subcritical quadratic and critical linear deviations suggests an inner transition scale of order \(\xi\sim\sqrt\delta\), despite the overall finite halo radius being continuous.

## C. Direct numerical test against the analytic singular coefficient

The positive 3D momentum integral and independent Poisson ODE were computed for \(\delta\in[0.2,10^{-5}]\), without replacing the profile by its exact-critical form. The analytic coefficient is \(\pi/(24I_{2,c})=0.1300973118\). We find:

| delta | Maximum occupation | \(\sqrt\delta\,\kappa(\delta)\) | First-zero radius \(r_t/R\) |
|---:|---:|---:|---:|
| 0.2 | 4.426 | 0.17749 | 78.425 |
| 0.01 | 99.406 | 0.13891 | 119.585 |
| 0.001 | 999.405 | 0.13280 | 126.927 |
| 0.0001 | 9999.405 | 0.13094 | 129.159 |
| 0.00001 | 99999.405 | 0.13036 | 129.857 |
| **Critical limit** | **infinite** | **0.130097** | **130.178** |

The density's *macroscopic* radius approaches a well-defined limit even while fmax and central curvature diverge. This is why just checking convergence of the rotation curve or virial residual cannot establish the validity of a critical-center **spectral/nonlinear** stability theorem.

## D. Priority consequences

For \(\delta>0\), verify the precise classical Antonov sufficient stability hypotheses and strong norms. For \(\delta=0\), investigate whether bounded equilibria converge in an appropriate \(L^p, p<3\) space and whether the gravitational linearized operator has a meaningful closed extension across an unbounded phase-space distribution. Do not infer the endpoint’s actual kinetic relaxation or coherent star formation from continuity of r_t/R alone.

References: original JCAP 2021 https://doi.org/10.1088/1475-7516/2021/11/055; Doremus–Feix–Baumann https://doi.org/10.1103/PhysRevLett.26.725; Pérez–Aly https://doi.org/10.1093/mnras/280.3.689; Guo–Rein https://doi.org/10.1007/s00220-007-0212-8.