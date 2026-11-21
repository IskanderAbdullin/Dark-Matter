# 16 — Why the smoothness threshold coincides with the Lane–Emden n=5 separatrix

**Analytic edge phase-plane reduction and independent numerical tests** for the self-consistent spherical Poisson equation of the positive F_k from note 15. Does not prove radius classification for every nonpolytropic inner profile.

The spherical equations are psi'=-Q/r² and Q'=r² x(psi), x=rho/rho0. Define tau=ln(r/R), u=-d ln(psi)/d tau=Q/(r psi), and w=r²x(psi)/psi. At the low-psi edge, x(psi)~C psi^n, with n=k+3/2. Hence the asymptotic autonomous system is

\[
\boxed{\dot u=u^2-u+w,\qquad
\dot w=w[2-(n-1)u].}
\]

For n>3 the singular-potential fixed point is u*=2/(n-1), w*=u*(1-u*). The Jacobian has trace T=(5-n)/(n-1) and determinant D=2(n-3)/(n-1). Therefore its complex eigenvalues have positive real part for n<5 and negative real part for n>5. The singular outer solution has

\[
\boxed{\psi(r)\propto r^{-2/(n-1)},\qquad Q(r)\propto r^{(n-3)/(n-1)}}.
\]

In particular k=3 means n=4.5, with eigenvalues +0.0714286 ± 0.923061i. Its tiny positive growth in *logarithmic radius* explains why a quasi-singular profile can persist over millions or billions of scale radii and eventually reach a finite edge. k=4 gives n=5.5 with eigenvalues −0.0555556 ± 1.05263i, making the singular outer phase-plane point attracting rather than repelling.

Independent log-radius integrations (eta=4, central psi0=4), measured on 10^7–10^10 R:

| k | n | Theoretical Q exponent | Numerical Q exponent |
|---:|---:|---:|---:|
| 4 | 5.5 | 0.55556 | 0.56752 |
| 5 | 6.5 | 0.63636 | 0.63813 |
| 6 | 7.5 | 0.69231 | 0.69296 |

For k=3 the radius is strongly nonmonotonic in central depth: at psi0=1 it is 45.81 R; at 2.5 it is 425.63R; at 2.75 it is 3658.73R; at 2.875 it is 4.45 million R; at 3.0 it is 1.26 million R; at 3.8 it is 2.95 billion R. This is one-parameter sensitivity of the mathematical equilibrium, not a fit to galaxy observations.

**Interpretation:** the edge-index crossing k=3.5 is a genuine analytic change in local Emden phase-space behavior. Our numerical k>=4 branches are consistent with the *known* n>5 infinite-extent / increasing-mass polytropic asymptotic. Yet deriving a strict all-radius, all-central-depth no-zero theorem for the *nonpolytropic Bose core* requires additional global comparison/phase-plane arguments; none is claimed.

Prior art: classic Lane–Emden threshold and perturbative solutions https://doi.org/10.1093/mnras/465.4.4265; Rein & Rendall (2000) https://doi.org/10.1017/S0305004199004193. The present connection to this Bose-cutoff family is our conditional application.
