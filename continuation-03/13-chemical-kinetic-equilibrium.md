# Chemical vs kinetic equilibrium and particle-number-dependent Bose stars

**Date:** 21 Nov 2026. Direct comparison of 2025/2026 primary publications and a conditional scaling analysis; not a complete rederivation of their coupled differential equations.

## 1. Two kinds of equilibration are different
*Chemical and kinetic equilibrium in cosmology: Facts and myths*, Phys. Rev. D 112 103509 (2025), DOI https://doi.org/10.1103/lfvd-bkk4, explicitly emphasizes that chemical equilibrium controls number densities while kinetic equilibrium controls distribution shapes, and either can cease before the other. This is particularly important when a complex scalar with approximate charge conservation, an ideal BE galactic \(\mu\), freeze-in production and Bose-star formation are discussed together:
- An elastic \(2\leftrightarrow2\) collision kernel can exchange momenta and, under suitable conditions, establish a Bose functional shape **while conserving particle number**. Chemical equilibrium with number-changing reactions would constrain or erase chemical potential depending on conserved charges.
- The original JCAP 2021 galactic BE DF presumes a global \(T\) and a \(\mu\), neither derived by virialization alone. Cosmological decoupling variables \(T_{b0}\) and galaxy virial temperature \(T\) describe completely different dynamical stages.
- A nonthermal Bose-stimulated cold production mechanism does not guarantee *late* kinetic equilibrium: the coupling responsible for production and that for elastic galactic scattering may differ.

## 2. GPP condensate with explicit particle number
Guzmán, Castellanos & Mastache, *Equilibrium halo solutions of the Gross–Pitaevskii–Poisson system: the role of the particle number*, EPJC 86 1103 (2026), DOI https://doi.org/10.1140/epjc/s10052-026-16373-w, retain mass, scattering length and total particle number \(N\) as independent inputs and map ground/excited/unbound branches. Their abstract notes an intermediate regime where quantum gradient pressure, gravity and repulsive interactions all contribute. This is an important qualification of the pure n=1 Thomas–Fermi radius law used in JCAP, which assumes a sufficiently strong interaction-dominated condensate. A fully condensed GPP galaxy is *not* the same as the original predominantly noncondensed mixed halo.

The original isolated TF Bose-star radius \(R_c=\pi\sqrt{\hbar^2a/(Gm^3)}\) is independent of star mass only after the kinetic gradient term is negligible. A dimensional criterion is obtained by comparing a kinetic energy scale \(K\sim N\hbar^2/(mR_c^2)\) with quartic interaction \(U_{\rm int}\sim gN^2/R_c^3\), so \(|U_{\rm int}|/K\sim 4\pi Na/R_c\), up to order-one profile factors. Therefore for fixed \(a,R_c\), the TF approximation may require \(Na/R_c\gg1\). This does not determine the astrophysical Bose-star mass function, stability or formation rate.

## 3. Very light hot relic abundance vs isolated bosonic mini-stars
Dror et al., *Indirect detection of hot dark matter*, JHEP05(2025)127, DOI https://doi.org/10.1007/JHEP05(2025)127, investigate **hot** light relics and their density profiles around CDM halos, with indirect-detection tests for axions mainly in \(\mathcal{O}(1-10)\) eV, not this model's 86 eV. This provides an independent dynamical treatment of a subdominant hot relic, but **not** evidence that a single 86 eV population makes up all dark matter and thermalizes into the original Bose halo. It also stresses that production and cosmic clustering must be studied before using galactic abundance as an input.

## 4. Required model cross-checks
- Keep chemical number density and kinetic shape separately in a comoving production solver; state exact conserved charge (if any).
- Derive star radius as a function of \(m,a,N\) *before* asserting a universal radius; check n=1 polytrope TF limit and the first-zero definition (factor-π audit).
- Solve excited and unbound GPP branches as possible alternatives; merely finding a stationary wavefunction does not establish cosmological production or long-lived dynamical stability.
- Distinguish thermal \(f(p)\), collisionless stationary \(f(E)\), and coherent \(\Psi\), including their normalizations and separate kinetic evolution.
- For any star mass function, compare dynamical crossing time and tidal disruption to halo age. No masses are predicted by the radius formula alone.

**Status:** conceptual map + conditional dimensional analysis, no observational constraints or microphysical Lagrangian selected.
