# Next test: residual hot population accompanying Bose-stimulated cold burst

**Scope:** computational experiment, *not* a reproduction of a 1↔2 collision integral, and not an observational likelihood.

The relevant paper by Wen Yin, JHEP 05 (2023) 180, https://doi.org/10.1007/JHEP05(2023)180 (sections 3–4, especially equations 1.1, 4.10–4.12), studies χ1↔χ2+φ with Bose enhancement, inverse reactions and expansion. A highly occupied cold peak near p_burst/T≈ηM1²/(2T²) is accompanied by an ordinary hot φ component. The article explicitly warns that an O(1–10%) hot component can matter, and shows examples with ≈5% or ≈30% warm components depending on parameters. Kodai Sakurai & Wen Yin, JHEP 03 (2025) 202, https://doi.org/10.1007/JHEP03(2025)202, also study stimulated emission in 2→2 thermal scattering.

**Our earlier pure α-scan is not the distribution computed in either article.** It did not account for a residual hot component. It is therefore essential to test whether a small hot fraction qualitatively changes the transfer function.

We define a deliberately simple *single-species*, number-conserving, bimodal initial function:
\[
f_{\alpha,w}(p)=(1-w)\alpha^{-3} f_{\rm BE}(p/\alpha T_b)
  +w f_{\rm BE}(p/T_b),
\quad f_{\rm BE}(x)=1/(e^x-1).
\]
Here \(0\le w\le1\) is the *particle-number* fraction assigned to the hot component, not a separate fermionic WDM species or the late galactic noncondensed fraction. Because each term integrates to the same reference particle number, the full 86-eV real-scalar abundance is fixed. Its relativistic energy fraction compared with the reference thermal spectrum is \((1-w)\alpha+w\), so
\[
\Delta N_{\rm eff}(\alpha,w)
  =[(1-w)\alpha+w]\Delta N_{\rm eff}^{\rm thermal},
\quad
\frac{\rho_{{\rm hot,rel}}}{\rho_{{\rm all,rel}}}
  =\frac{w}{w+(1-w)\alpha}.
\]
At α=.02,w=.05, just 5% of particles carry ~72.5% of the *early relativistic energy* in this toy model. This does **not** mean a 72.5% hot dark matter mass fraction today.

The [code](burst_residual_hot.py) explicitly computes the mixed PSD in the CLASS convention with the mandatory \(1/(2\pi)^3\) cell factor, checks both number and radiation-energy moments, and calculates ten linear CLASS transfer spectra for α=.02,.01 and w=0,.01,.05,.1,.2. It uses a *single* noncold species rather than pretending the fractions are two separately conserved fundamental fields. With the collisionless linear equations, phase-space shape itself distinguishes the two populations. A subset is recomputed using an independent PSD grid.

The key observable diagnostics will be P(k,z=5)/P_CDM at k=1,5,10,20 inverse Mpc and the **first** power-half crossing \(P/P_{\rm CDM}=1/2\), when that crossing exists. Oscillations or numerical ringing at deeply suppressed k are reported, not removed or assumed physical. Our model's cold component shape remains artificial; the Yin paper includes dynamical momentum-dependent gain, loss, Bose stimulation, Pauli blocking, finite mother mass and hot-tail evolution. To genuinely connect microscopic physics to late Bose halos one must solve those **unintegrated** kinetic equations and feed the resulting actual f(q) into CLASS.

**Status:** executable hypothesis and public CI; no claim about the ten numerical results until actual successful workflow execution and convergence evaluation. No observational constraints imported or inferred automatically.
