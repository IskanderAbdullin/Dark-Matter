# 02. Collisionless relic shape and the absence of automatic late Bose thermalization

**Research note, 2026-11-21.** Grounded in source Abdullin & Popov 2021 eqs. (distribBE), (DMabundance), (mUpper), (ncEnergy); modern kinetic-decoupling literature independently checked. The conditional statements below do not exclude asymmetric/BEC/nonthermal scenarios.

## A. A relativistically frozen Bose relic is NOT a late nonrelativistic Bose thermal distribution

For a single scalar with relativistic decoupling at \(a_d\), negligible chemical potential, and thereafter collisionless homogeneous FLRW transport, the comoving momentum \(q=ap\) is conserved. If its initial occupation is \(f_d(p_d)=[\exp(p_dc/kT_d)-1]^{-1}\), then at arbitrary later \(a\):
\[
\boxed{f(p,a)=
\left[\exp\!\left(\frac{pc}{kT_d a_d/a}\right)-1\right]^{-1}.}\tag{1}
\]
At \(p\ll mc\) this remains a function of **linear** \(p\) in the exponent, not a local nonrelativistic equilibrium BE distribution
\[
f_{\rm NR,eq}(p)=
\left[\exp\!\left(\frac{p^2}{2mkT_{\rm kin}}-\frac{\mu_{\rm kin}}{kT_{\rm kin}}\right)-1\right]^{-1}.
\tag{2}
\]
Equality for an open momentum interval would force \(a_1 p+b_1=a_2p^2+b_2\) as functions of \(p\), impossible with nonzero finite temperatures. This simple functional-shape contradiction holds even if the **actual particles become nonrelativistic** after the freeze-out. One must add elastic kinetic thermalization, collective processes, self-interaction, or choose a primordial nonthermal spectrum to obtain the source's late BE(kinetic energy) distribution. Virialization/phase mixing alone does not by itself establish quantum kinetic BE equilibrium or uniform chemical potential.

The exact number-weighted moments of distribution (1) in its nonrelativistic kinematic limit, in units \(T_p=kT_da_d/a\) (momentum energy scale), are
\[
\langle pc\rangle/T_p=\frac{3\zeta(4)}{\zeta(3)}=2.701178,\quad
\frac{\langle (pc)^2\rangle}{T_p^2}=\frac{12\zeta(5)}{\zeta(3)}=10.351534,
\]
\[
\frac{\sigma_{\rm 1D}}{cT_p/(mc^2)}
=\sqrt{\frac{4\zeta(5)}{\zeta(3)}}=1.857555,\quad
\gamma_2^{\rm frozen}
=\frac{9}{2}\frac{\zeta(7)\zeta(3)}{\zeta(5)^2}-3
=2.072844.
\]
This last kurtosis differs from the **nonrelativistic critical Bose equilibrium** value \(1.906876\) derived in continuation-03; a matching mass, average speed or momentum scale cannot make the entire spectra equal.

## B. A separate theorem: no condensation by homogeneous free expansion after *nonrelativistic* decoupling

Suppose instead a gas is already nonrelativistic and has \(f_d(p)=[z_d^{-1}\exp(p^2/(2mkT_d))-1]^{-1}\), \(0<z_d<1\), when kinetic collisions cease. Collisionless Hubble redshift \(p\mapsto p(a_d/a)\) keeps it *exactly* in the NR Bose family if parameterized as
\[
T(a)=T_d(a_d/a)^2,\qquad z(a)=z_d.
\]
Then \(n\propto a^{-3}\), \(\lambda_{\rm th}=\sqrt{2\pi}\hbar/\sqrt{mkT}\propto a\), and
\[
\boxed{n(a)\lambda_{\rm th}(a)^3
=\mathrm{Li}_{3/2}(z_d)=\text{constant}<\zeta(3/2).}
\tag{3}
\]
There is no spontaneous crossing into condensation from *only homogeneous free expansion*. An isentropic NR ideal Bose gas maintained in local equilibrium likewise conserves entropy per particle, \(s/(nk)=\frac52\,\mathrm{Li}_{5/2}/\mathrm{Li}_{3/2}-\ln z\), giving constant fugacity along reversible adiabats under standard assumptions. This statement does **not** prohibit condensation caused by inhomogeneous gravitational compression, particle-production asymmetry, entropy exchange, dissipative cooling, changes of effective degrees of freedom or mixed quantum field evolution.

## C. Cosmological abundance normalization: convention audit, not a categorical correction

Source eq. (DMabundance) writes, for the noncondensed component,
\[
\rho_{\rm DM} =\frac{2\zeta(3)m}{\pi^2 g_d^*}\left(\frac{kT_{\gamma0}}{\hbar c}\right)^3
\left(\frac{T_c}{T_d}\right)^3.
\tag{4}
\]
Under a DIFFERENT *fully specified* benchmark — **one real bosonic degree of freedom**, equilibrium \(\mu=0\) relativistic freeze-out, no condensate, no subsequent entropy injection to the bosons, standard SM photon+neutrino entropy bookkeeping — the late boson temperature satisfies \(T_{b0}/T_{\gamma0}=(g_{*S,0}/g_{*S,d})^{1/3}\), \(g_{*S,0}\simeq3.909\), so
\[
n_{b0}=\frac{\zeta(3)}{\pi^2}
 \left(\frac{kT_{\gamma0}}{\hbar c}\right)^3
 \frac{3.909}{g_{*S,d}}.\tag{5}
\]
If one *identifies* \(g_d^*=g_{*S,d}\) and sets \(T_c/T_d=1\) in (4), the number-density coefficients differ by \(3.909/2=1.9545\). At \(m=86\) eV, \(g_{*S,d}=61.75\), (5) gives \(n_b\simeq13.0003\;\mathrm{cm^{-3}}\), \(\Omega h^2\simeq0.10610\), while using (4) with no extra-condensate factor gives \(n\simeq6.6514\;\mathrm{cm^{-3}}\), \(\Omega h^2\simeq0.05428\). **Do not call this a proven error of (4):** the source may use a different \(g_d^*\) convention, an implicit species multiplicity, condensate history \(T_c/T_d\), conserved number asymmetry or entropy treatment, all of which require checking its cited Boyanovsky et al. (2008) reference and 2019/2025 author followups. The two formulas cannot be compared as identical physical scenarios without those details.

## D. Research consequences

1. For primordial 86 eV: write the **momentum-dependent** production rate, evolution of the collisionless distribution through relativistic/NR transition, and precise effective degrees of freedom; no late thermal-shape assumptions by fiat.
2. Distinguish \(g_b=1\) real scalar from \(g_b=2\) complex particle/antiparticle; a complex scalar with U(1) charge and chemical potential \(\mu\neq0\) requires a separate abundance calculation.
3. Compare initial \(f(q)\), galactic stationary \(f(E)\), and Bose-star coherent field \(\Psi\) **as three physically different objects**. Mathematically a fitted \(f(E)\) can be stationary under Vlasov without being a thermally populated DM gas.
4. Use the 6D Liouville/Casimir level-set test in paper 03, taking care not to confuse *fine-grained conserved occupancy* with a coarse-grained galaxy halo distribution.
5. A kinetic relaxation mechanism derived from the Proukakis–Rigopoulos–Soto GP/Boltzmann/Poisson framework is a credible bridge, but its validity and timing for 86 eV have not been demonstrated.

**Literature anchors:** Boyanovsky, de Vega & Sanchez (2008) https://doi.org/10.1103/PhysRevD.77.043518; *Chemical and kinetic equilibrium in cosmology: Facts and myths* (2025) https://doi.org/10.1103/lfvd-bkk4; Proukakis–Rigopoulos–Soto (2023) https://doi.org/10.1103/PhysRevD.108.083513; *Momentum distributions of cosmic relics* (2022) https://arxiv.org/abs/2201.06456; source JCAP 2021.
