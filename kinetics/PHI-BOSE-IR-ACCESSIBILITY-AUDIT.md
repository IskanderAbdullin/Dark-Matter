# Exact parent boost required to reach the massive 86-eV scalar infrared (2026-11-23)

## Decisive physical/numerical question

Earlier [parent-cutoff IR support](PHI-BOSE-PARENT-CUTOFF-KINEMATIC-IR-AUDIT.md) established a strict finite-parent-cutoff daughter energy floor. The missing inverse question is more useful for the physical model: **how energetic must a parent χ1 be for its two-body decay χ1→χ2+φ to reach a specified very low laboratory φ energy?** In particular, can an illustrative 1-GeV-temperature parent populate a scalar *at rest* through on-shell 1→2 decay at appreciable thermal rate?

This stage derives the exact threshold, including a critical subtlety: the minimum scalar energy **at a fixed parent momentum is nonmonotone** after the boost that produces a scalar at rest. The minimum over all parent momenta up to a cutoff is therefore not always attained at that cutoff. No finite-grid interpolation or invented zero mode is involved.

## Exact on-shell Lorentz derivation

Let \(m_1>m_2+m_\phi>0\). In the parent rest frame,

\[
E_\phi^*=\frac{m_1^2+m_\phi^2-m_2^2}{2m_1},\quad
q_\phi^*=\sqrt{(E_\phi^*)^2-m_\phi^2},\quad
\eta_*=\operatorname{arsinh}(q_\phi^*/m_\phi).
\]

Write the parent lab momentum \(p=m_1\sinh\eta\) and rapidity \(\eta\ge0\). At backward scalar emission \(\cos\theta^*=-1\),

\[
E_{\phi,\min}(p)=m_\phi\cosh(\eta-\eta_*).
\]

The exact **smallest** parent momentum allowing any scalar of laboratory energy \(E_\phi\le K\), \(m_\phi\le K<E_\phi^*\), is therefore

\[
\boxed{p_{\rm required}(K)=m_1\sinh\!\left[
\eta_*-\operatorname{arcosh}(K/m_\phi)\right].}
\]

For \(K\ge E_\phi^*\), \(p_{\rm required}=0\). For \(K<m_\phi\), there is no on-shell scalar. For a cutoff \(P\), the global daughter support floor is \(m_\phi\) if \(P\ge p_{\rm required}(m_\phi)\), and otherwise \(E_{\phi,\min}(P)\). This explicitly avoids incorrectly extending the fixed-p decreasing branch past its minimum. A cancellation-free positive evaluation is

\[
E_{\phi,\min}(p)=
\frac{m_\phi^2p^2+m_1^2(E_\phi^*)^2}
{m_1(E_1 E_\phi^*+p q_\phi^*)}.
\]

The physical daughter rest momentum \(q_\phi=0\) requires a parent with \(p_*=m_1 q_\phi^*/m_\phi\), not \(p\sim m_1\). For the intermediate regime \(m_\phi\ll K\ll E_\phi^*\), the leading ultrarelativistic inversion is \(p_{\rm required}(K)\simeq m_1(E_\phi^*+q_\phi^*)/(2K)\simeq m_1E_\phi^*/K\). The exact rapidity formula, not this approximation, is used for the table.

## Reproducible numeric benchmark

Use \(m_1=0.1\), \(m_2=0.02\), \(m_\phi=86\,\mathrm{eV}=8.6\times10^{-8}\) GeV. The parent-rest scalar energy is \(E_\phi^*=0.048000000000037\) GeV.

| Requested lab scalar energy K | Minimum parent p required |
|---:|---:|
| 86 eV (scalar at rest) | **55,813.9534883 GeV** |
| 172 eV (2mφ) | 14,955.3037634 GeV |
| 860 eV (10mφ) | 2,797.70952099 GeV |
| 1 keV | 2,404.45408546 GeV |
| 10 keV | 240.004427348 GeV |
| 100 keV | 23.999900271 GeV |
| 1 MeV | 2.39895833777 GeV |

For hard parent cutoffs \(P=2,8,16,24\) GeV, the respective **minimum** scalar energies are \(1.19925094\) MeV, \(299.988288\) keV, \(149.998548\) keV, and \(99.9995845\) keV. The corresponding scalar momenta near these values are practically equal to the energies because \(m_\phi\ll K\).

**Conditional thermal implication:** at parent temperature \(T_1=1\) GeV and fixed finite chemical potential, a Fermi-Dirac tail at \(p\ge2404\) GeV is exponentially suppressed at least by a factor of order \(e^{-2404}\) relative to a unit-scale parent occupation, while the 86-eV-rest threshold is at \(p\approx55814\) GeV. This is a kinematic accessibility and occupation-tail statement, **not a bound on the total Bose-stimulated collision rate near criticality** (Bose factors, phase space and actual cosmological history need separate control). The scalar rest energy is a single measure-zero state in a continuous excited distribution; this analysis does not derive or exclude a dynamically coupled delta-function condensate.

## Scientific consequence and negative control

A physical finite-volume Yukawa 1→2 collision solver with a parent cutoff of 2–24 GeV has **exactly zero direct event support** for the lower scalar energy interval \(m_\phi\le E_\phi<K_{\min}(P)\). This is true regardless of how accurate its energy, two dark charges and H theorem are. On this **restricted 1→2 event operator**, any perturbation of scalar occupations supported strictly below \(K_{\min}(P)\) is invisible to all its collision events: an infinite-dimensional spurious collision nullspace exists if no other processes are included. Thus a finite-parent-only calculation cannot demonstrate physical thermalization or determine a Bose zero-mode through this operator. An event graph which falsely deposits direct Yukawa daughters into that forbidden interval violates exact on-shell kinematics.

This **does not** mean those low scalar modes are truly decoupled in the full Universe: 2→2 scatterings, off-shell/thermal processes, cosmological redshift, interactions with a condensate and higher-energy parent populations can couple them. Nor does it imply the total 1→2 rate is small; most daughters in the thermal example are produced at much higher scalar energy.

## Reproduction and next test

- [Source](phi_bose_ir_accessibility.py), [focused regressions](../tests/test_phi_bose_ir_accessibility.py).
- \`python -m kinetics.phi_bose_ir_accessibility\`
- \`python -m pytest -q tests/test_phi_bose_ir_accessibility.py\`

Independent local direct evaluation of the Lorentz-inverted threshold and the rationalized boosted lower endpoint agrees to \(\lesssim 2.3\times10^{-15}\) relative at seven test energies; **public GitHub Scientific CI must be verified separately**.

**Next:** compare the exact threshold against the scalar bin edges of the existing discontinuous entropy-dual P1 solver; flag bins with no physical 1→2 event support as genuine cutoff nullspace, not proof of a BEC or a converged freeze-out. Quantify thermal tail versus Bose critical amplification before using any 86-eV scalar near-zero-mode claim. Preserve CLASS-ready=false until full physical coupled SM+FRW and observational validation.
