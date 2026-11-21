# 11. Universal shallow-cutoff scaling limit of critically occupied finite Bose halos

**Theoretical note, 21 November 2026.** A new reduction and independent numerical experiment for our *phenomenological* positive lowered-Bose distribution. This is not a galaxy fit, new particle detection, dynamical stability theorem, or evidence that a cosmological tidal cutoff of this type forms.

## Set-up and rescaling

For a Bose occupation lowered at escape energy,
\[
f(q;\psi,\eta)=\frac{1}{e^{\eta-\psi+q^2}-1}
-\frac{1}{e^\eta-1}\quad (q^2<\psi),\qquad f=0\ {\rm else},
\]
take central saturation \(\psi_0=\eta\), density integral \(I(\psi;\eta)=\int_0^{\sqrt\psi}q^2f(q;\psi,\eta)dq\), and Poisson equation \(\psi'=-Q/\xi^2,\ Q'=\xi^2 I(\psi)/I_c\).

In the **shallow-cutoff limit** \(\eta\to0^+\), define
\[
U=\psi/\eta,\quad q=\sqrt\eta\,v,\quad s=\xi/\sqrt\eta,\quad {\cal Q}=Q/\eta^{3/2}.
\]
Then the phase-space integral has a universal limit
\[
I(\eta U;\eta)=\sqrt\eta\,j(U)+O(\eta^{5/2}),\qquad
\boxed{j(U)=U^{5/2}\int_0^1\frac{t^2(1-t^2)}{1-U(1-t^2)}dt.}
\]
Equivalently \(j(U)=\sqrt U-\sqrt{1-U}\arctan\sqrt{U/(1-U)}-U^{3/2}/3\) for \(0<U<1\), continuous at \(U=1\) with \(j(1)=2/3\). The integral form avoids severe numerical cancellation near \(U=0\).

Thus \(I_c=(2/3)\sqrt\eta-\eta^{5/2}/90+O(\eta^{9/2})\), and the entire normalized density tends to \(x(U)=(3/2)j(U)\). The flux system becomes **independent of \(\eta\)**:
\[
\boxed{U'=-{\cal Q}/s^2,\qquad {\cal Q}'=s^2(3/2)j(U),\qquad
U(0)=1,\ {\cal Q}(0)=0.}
\]
At its first zero \(s=s_t\), the physical dimensionless halo radius and mass satisfy
\[
\boxed{r_t/R\sim s_t\sqrt\eta,\qquad
Q_t\sim{\cal Q}_t\eta^{3/2}.}
\]

## Independent limiting numerical solution

Using Gauss–Legendre quadrature for \(j\) and a high-accuracy spherical flux ODE gives
\[
\boxed{s_t=11.8605981867,\qquad {\cal Q}_t=2.8955760808.}
\]
The exact critical central cusp slope rescales to
\[
\boxed{\sqrt\eta[-x_\xi'(0^+)]\to 3\pi/(4\sqrt6)=0.9619123726,}
\]
while at the **outer** edge \(j(U)\sim(2/15)U^{5/2}\) and \(x(U)\sim U^{5/2}/5\). This is not a single polytropic law throughout the halo.

| \(\eta\) | Direct \(r_t/(R\sqrt\eta)\) | Direct \(Q_t/\eta^{3/2}\) | \(\sqrt\eta[-x_\xi'(0^+)]\) |
|---:|---:|---:|---:|
| 0.5 | 11.955259 | 2.900580 | 0.965909 |
| 0.2 | 11.875694 | 2.896374 | 0.962553 |
| 0.1 | 11.864376 | 2.895775 | 0.962073 |
| 0.05 | 11.861552 | 2.895625 | 0.961952 |
| 0.02 | 11.860756 | 2.895583 | 0.961919 |
| 0.01 | 11.860642 | 2.895577 | 0.961914 |
| **Limit** | **11.860598** | **2.895576** | **0.961912** |

These six finite-\(\eta\) solutions came from the **separate original** lowered-Bose quadrature/Poisson code rather than from inserting the universal solution by hand. Three convergence checks pass.

## Connection to the fine-grained mixing constraint (paper 10)

The phase-density normalization \(\nu_\eta=4I_c/\sqrt\pi\sim(8/3\sqrt\pi)\sqrt\eta\), and \(Q_t\sim2.8955761\,\eta^{3/2}\). Therefore the initial thermal-patch-to-final phase-volume high-occupation margin in paper 10 follows the **independent analytic scaling**
\[
\boxed{s_\eta=\frac{\nu_\eta Q_t}{3\sqrt3\zeta(3)\sqrt{\pi/2}}
\sim0.5564958486\eta^2.}
\]
The strong \(\eta^2\) suppression explains why a shallow-cutoff finite halo cannot be the *all-mass-retained, purely collisionlessly mixed* descendant of that particular homogeneous relativistic thermal Bose relic. It is a conditional *phase-volume* obstruction, not a bound on Bose particle mass or on physically open tidal halos.

## Physical meaning and required next step

This one-parameter equilibrium family now has two controlled limiting directions: \(\eta\to0\), a universal finite compact rescaled profile; \(\eta\to\infty\), recovery of the untruncated infinite-mass isothermal envelope. Neither limit proves that a galaxy picks a specific \(\eta\), forms its critical Bose occupation, or is stable. The next task is to combine a kinetic escape/evaporation process with globally stable Vlasov–Poisson states, and then compare cosmological formation and observations.

Earlier work: Abdullin & Popov JCAP 2021 055 https://doi.org/10.1088/1475-7516/2021/11/055; finite-cutoff Vlasov background https://doi.org/10.1093/mnras/stad117; Guo–Rein King stability https://doi.org/10.1007/s00220-007-0212-8.
