# 15 — A positive hierarchy of higher-order Bose escape-energy truncations

**Source-derived setup:** Abdullin & Popov, JCAP 11 (2021) 055, eqs. distribBE, nonconddens, inhR. **New calculations for this research project:** positive Taylor-remainder family, edge power law and independent self-consistent radial Poisson profiles. No cosmological formation or thermodynamic equilibration asserted.

Let B(x)=(exp(x)-1)^(-1), x=(E-mu)/(kBT)>0, eta=(E_t-mu)/(kBT)>0, and t=(E_t-E)/(kBT)>0. For integer k>=1 set F_k=0 if E>=E_t, and for E<E_t let

\[
F_k(E)=B(\eta-t)-\sum_{j=0}^{k-1}\frac{\operatorname{Li}_{-j}(e^{-\eta})}{j!}\,t^j.
\]

Since B(eta-t)=sum_{n>=1}e^(-n eta)e^(n t), an exact positive representation is

\[
\boxed{F_k(t)=\sum_{n=1}^\infty e^{-n\eta}
\left(e^{nt}-\sum_{j=0}^{k-1}\frac{(nt)^j}{j!}\right)>0.}
\]

Each bracket is a positive exponential Taylor remainder for t>0, and its t derivative is also positive. Therefore df/dE<0 throughout the occupied region. F_k has C^(k-1) regularity at the energy escape boundary (not generally C^k); k=1 recovers the linearly lowered distribution, k=2 the previously introduced positive tangent subtraction.

For q=p/sqrt(2mkBT), psi=(E_t-mPhi)/(kBT), the positive density moment is I_k(psi,eta)=int_0^sqrt(psi) q²F_k(psi-q²)dq. Let a_k=Li_{-k}(exp(-eta))/k!>0. Then at the edge

\[
\boxed{I_k(\psi,\eta)\sim \frac{a_k}{2} B\!\left(\frac32,k+1\right)\psi^{k+3/2}.}
\]

The analogous J_k=int q⁴F_k dq yields dJ_k/dpsi=(3/2)I_k by integration by parts, proving hydrostatic pressure consistency with the same stationary collisionless Poisson DF. **Its *outer* effective polytropic index is n=k+3/2, not an exact full-halo polytrope.**

Numerical self-consistent spherical Poisson integrations for eta=4, central psi0=3.8:

| k | Boundary class | n | First zero r_t/R | Mass Q_t |
|---:|---|---:|---:|---:|
| 1 | C0 | 2.5 | 27.37249 | 18.75118 |
| 2 | C1 | 3.5 | 78.42542 | 17.90742 |
| 3 | C2 | 4.5 | 2.9517 × 10^9 | 3.3228 × 10^4 |
| 4 | C3 | 5.5 | No zero out to 10^10 R | Q(10^10 R)≈7.963×10^5 |

For exact central saturation psi0=eta=4, k=1 yields 40.2988R, k=2 130.1784R, k=3 approximately 1.6011×10^9 R; k>=4 up to k=6 gave **no zero within the computed range**.

This is a strong warning about adding smoothness to an arbitrary cutoff: it can generate enormous and apparently noncompact envelopes. The k>=4 numerical no-zero statement out to 10^10R is **not** a global theorem for every initial depth or every Bose microphysics. k=3 radii were checked with different ODE tolerances. At the originally printed F583-4 scale R=0.6kpc the k=3 profile would have illustrative radius ~960Gpc; this is not an observation or a new galaxy fit because R,V0 must be refitted for each new DF.

Related mathematics: Rein & Rendall compact-support theorem https://doi.org/10.1017/S0305004199004193; Lane–Emden n=5 threshold https://doi.org/10.1093/mnras/465.4.4265; caution about different edge-regularity damping results https://doi.org/10.1007/s00205-025-02114-y.

Run reproduction: reproduce_cutoff_ladder.py plus the more complete source scripts bundled with the downloadable research archive.
