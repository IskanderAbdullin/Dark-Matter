# 17 — Existing Antonov / Doremus–Feix–Baumann theory resolves part of the purported stability gap

**Literature correction for the research program, 21 November 2026.** Previous notes called global stability of every lowered Bose halo an entirely open eigenmode problem. That is unnecessarily pessimistic for **bounded, isotropic, compact, strictly subcritical single-component collisionless equilibria**.

Established results:
- Doremus, Feix & Baumann (1971), https://doi.org/10.1103/PhysRevLett.26.725, prove sufficiency of df(E)/dE<0 for the spherical (radial) perturbation modes in appropriate regular classes.
- Pérez & Aly (1996), https://doi.org/10.1093/mnras/280.3.689, https://arxiv.org/abs/astro-ph/9511103, discuss nonradial and spherical perturbations, specifying preserving generators and spatial-density conditions, and establish stability of isotropic decreasing-energy models under stated hypotheses.
- The open textbook https://galaxiesbook.org/chapters/IV-01.-Formation-of-Dark-Matter-Halos_3-Virialization.html summarizes the Antonov nonradial and DFB radial criteria.
- Guo & Rein (2007), https://doi.org/10.1007/s00220-007-0212-8, prove nonlinear stability of the **classical King** distribution for a restricted class of dynamically accessible spherical perturbations; this is NOT a general nonlinear proof for arbitrary lowered Bose.
- Recent MNRAS mathematical modeling (2026), https://doi.org/10.1093/mnras/stag050, explicitly calls the monotone isotropic f(E) criterion a sufficient (not necessary) linear stability condition.

For the positive Bose family with t=(E_t−E)/kBT, one has

\[
\boxed{\frac{dF_k}{dE}=-\frac{1}{k_BT}
\sum_{n=1}^\infty n e^{-n\eta}
\left[e^{nt}-\sum_{j=0}^{k-2}\frac{(nt)^j}{j!}\right]<0}
\]

(for k=1, the bracket is simply e^{nt}). In addition the physical self-consistent density decreases outward because dI_k/dpsi>0 and dpsi/dr<0.

For psi0<eta, fmax is **finite** and the k=2 model is at least C1 at the energy cutoff; its eta4, psi0=3.8 spherical Poisson solution has finite outer radius ≈78.425R and total dimensionless mass ≈17.907. These are all natural candidates to apply the established **classical collisionless linear stability criterion**, after checking all regularity and perturbation-class hypotheses. No need to rediscover general radial Antonov theory from zero.

At exact saturation psi0=eta, however, F_k(0,0) diverges, and typical boundedness/regularity assumptions fail. Stability of the *singular critical endpoint* is not obtained merely by taking the limit of the bounded family. See note 18 for the quantitative nonuniform approach.

Moreover **linear Vlasov stability is not thermodynamic/kinetic stability**, does not entail damping, condensation, an H-theorem on a galactic timescale, Bose-star survival, nonlinear stability in every perturbation class or production from a cosmic relic. Nor does the theorem automatically apply to a two-component gas+stars system, Hartree–Fock self-field or a baryon-forced anisotropic halo.

**Revised priority:** carefully match regularity and perturbation assumptions of established theorems for the bounded k=2 Bose branch, then focus new mathematical/spectral work on the exact critical endpoint and dynamical two-component model. This is a literature-derived correction, not a new stability theorem.
