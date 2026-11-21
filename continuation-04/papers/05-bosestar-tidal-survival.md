# 05. Microscopic Thomas–Fermi validity, halo tides and extended lensing impose radically different mass scales

**Status:** conditional analytic scales and independent SI arithmetic, not a derived star mass function or an observed bound. Starts strictly with JCAP 2021 GP/TF equations (g), (crosssec_via_g), (BECdens), (Rc), (RbsSIDM), and the factor-\(\pi\) audit in continuation-03.

## A. Fiducial condensate size and TF consistency

For a repulsive contact interaction with \(\sigma=8\pi a^2\) and \(g=4\pi\hbar^2a/m\), full first-zero radius in the \(n=1\) TF polytrope:
\[
R_c=\pi\sqrt{\hbar^2 a/(Gm^3)}.
\]
At *illustrative* \(m=86\) eV and \(\sigma/m=1.25\,\mathrm{cm^2/g}\): \(a=8.73210\times10^{-19}\) m and \(R_c=133.4461\) AU. The later source formula (RbsSIDM) instead gives \(42.4772\) AU because it omits the \(\pi\) between \(k^{-1}\) and first-zero radius. No claim that these parameters are measured.

A rough GP energy ratio in a spherical state is
\[
\frac{U_{\rm int}}{K}\sim\frac{4\pi Na}{R_c},\qquad N=M/m,
\]
up to profile constants. Therefore a **dimensional TF crossover** has a tiny nominal mass \(M_{\rm TF,scale}\sim mR_c/a\approx3.50\times10^{-3}\) kg. This small number does *not* mean a gram-scale object is physically self-gravitating or forms as a coherent star — it only says the continuum interaction term may dominate the gradient term at this assumed radius when the actual self-consistency, gravitational binding and formation constraints also hold. In particular treating \(R_c\) as fixed for masses below TF validity is invalid.

## B. Formation crossing-time sanity check

For an assumed self-gravitating star of fixed radius \(R_c\), \(t_{\rm dyn}=\sqrt{R_c^3/(GM)}\). Requiring \(t_{\rm dyn}<10\) Gyr gives a separate indicative scale
\[
M > M_{10{\rm Gyr}}=\frac{R_c^3}{G(10\,\mathrm{Gyr})^2}\approx1.20\times10^{15}\,\mathrm{kg}.
\]
This only compares orbital/dynamical times to an age; the actual condensation or nucleation time can be *much* larger. No precise lower mass from galaxy formation has been inferred.

## C. Galaxy tidal field may impose a much stronger floor

For a spherical host with approximately **flat** local circular speed \(V_c(r)\) and a compact object on a **circular orbit**, the Jacobi/tidal radius obeys \(r_t\simeq [GM/(2\Omega^2)]^{1/3}\), \(\Omega=V_c/r\), so \(R_c\lesssim r_t\) requires
\[
\boxed{M\gtrsim M_{\rm tidal}
=\frac{2V_c(r)^2R_c^3}{G r^2}.}\tag{1}
\]
At illustrative **host** \(V_c=44.1\) km/s at \(r=1\) kpc and \(R_c=133.4461\) AU:
\[
M_{\rm tidal}\simeq4.87\times10^{20}\,\mathrm{kg}
\simeq2.45\times10^{-10}M_\odot.
\]
This exceeds \(M_{10{\rm Gyr}}\) by ~\(4.07\times10^5\) and the formal \(mR_c/a\) scale by ~\(1.39\times10^{23}\)! These orders demonstrate that microscopic TF consistency is utterly insufficient to infer survivable Bose-star masses.

**Assumptions matter:** the factor 2 applies to a spherical locally flat rotation curve and circular orbit, not to all galactic potentials, eccentric orbits, or the actual inferred Vc of F583-4 at 1 kpc. Bose stars in a core may experience compressive tidal fields; cuspy/baryon-dominated regions differ. The true condition requires evaluating host potential curvature, orbital distribution, finite density gradients and stellar structure. We deliberately use a simple conditional scale, not a universal lower bound.

## D. Observational consequences depend on object size

With a point-mass reference Einstein radius \(R_E=\sqrt{4GM D_{\rm eff}/c^2}\), choose solely for illustration \(D_{\rm eff}=8\) kpc (e.g. \(D_L=10\) and \(D_S=50\) kpc). For \(M=M_\odot\), \(R_E=8.07173\) AU, hence \(R_c/R_E=16.533\). The crossover \(R_E=R_c\) occurs at \(M\approx273.3M_\odot\) if the TF radius remains constant. Objects below that mass can be strongly extended lenses rather than point masses. A full lensing calculation requires the projected \(\rho(r)\), finite source size, lens impact parameter, telescope cadence, the host spatial distribution and **a mass function \(dN_{\rm BS}/dM\)**. Extendedness does not automatically imply a detectable signal.

## E. A decisive research direction

Simulate GP/kinetic Bose-star populations or derive a formation/survival mass function first; then include tidal stripping, evaporation and mergers in the Galactic environment, and finally predict astrometric/photometric lensing by the **extended** \(n=1\) polytrope. A radius formula alone cannot explain dwarf galaxies, missing satellites or lensing event rates.

Prior art: Abdullin & Popov JCAP (2021); Guzmán et al. *Equilibrium halo solutions ... role of particle number* (2026) https://doi.org/10.1140/epjc/s10052-026-16373-w; Bhandari et al. (2026) https://doi.org/10.1088/1475-7516/2026/09/072; extended HSC lensing https://doi.org/10.1103/PhysRevD.102.083021; Chavanis review (2025) https://doi.org/10.3389/fspas.2025.1538434.
