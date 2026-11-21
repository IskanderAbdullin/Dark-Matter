# Extended-object microlensing: corrected Bose-star radius changes the point-lens regime

**Source:** original JCAP 2021 formulas for repulsive TF Bose stars, with a distinct factor-π radius convention audit in note 01. Numbers here are conditional arithmetic for one illustrative σ/m upper-limit scale, not predicted object masses or event rates.

## 1. Lensing geometry
For lens distance \(D_L\) and source distance \(D_S>D_L\), define \(D_{\rm eff}=D_L(D_S-D_L)/D_S\). The physical Einstein radius in the lens plane of a **point mass** M is
\[
R_E=\sqrt{4GM D_{\rm eff}/c^2}.
\]
For a star with physical radius \(R_c\), if \(R_c\gtrsim R_E\) the object's projected mass profile affects magnification and centroid shift: point-lens exclusion curves cannot simply be read as bounds on this population. Large \(R_c/R_E\) alone does not predict a specific light curve or absence of lensing.

Use as an illustration \(D_L=10\) kpc and \(D_S=50\) kpc, giving \(D_{\rm eff}=8\) kpc. For 1 solar mass, \(R_E=8.071728\) AU.

## 2. JCAP conditional radius and threshold mass
At \(m=86\) eV, \(\sigma/m=1.25\) cm²/g and identical boson s-wave \(\sigma=8\pi a^2\), the *consistent first-zero* TF radius from source equations is \(R_c=133.4461\) AU, not 42.4772 AU (the latter is the length \(R_c/\pi\) in the later source equation).
The relevant dimensionless size ratio is
\[
\mathcal E\equiv R_c/R_E
=\frac{133.4461}{8.071728}\left(\frac{M}{M_\odot}\right)^{-1/2}
=16.5325(M/M_\odot)^{-1/2}.
\]
Threshold \(\mathcal E=1\):
\[
M_{\rm cross}=\frac{c^2R_c^2}{4GD_{\rm eff}}
=273.3\,M_\odot
\quad(D_{\rm eff}=8\text{ kpc}).
\]
Thus an **assumed** 1-Msolar Bose star of this radius is much more extended than the point-mass Einstein scale; a 300-Msolar object would be near crossover. Both masses are pure trial inputs. Changing the illustrative σ/m changes R_c as the fourth root and the crossover mass as the square root; a physically derived cross section and mass function are indispensable.

| Assumed compact-object mass | \(R_c/R_E\) |
|---:|---:|
| 0.1 M_sun | 52.28 |
| 1 M_sun | 16.53 |
| 10 M_sun | 5.228 |
| 100 M_sun | 1.653 |
| 300 M_sun | 0.955 |
| 1000 M_sun | 0.523 |

## 3. Independent gravitational timescale
The dynamical time proxy \(t_{\rm dyn}=\sqrt{R_c^3/(GM)}\) is a *necessary timescale scale*, not a proven nucleation time. With 133.446 AU:
- 1 M_sun: ~245 years.
- \(10^{20}\) kg: ~34.6 million years.
- \(10^{15}\) kg: ~10.9 billion years.
An age-scale comparison can flag candidates that cannot easily dynamically rearrange, but actual star formation can be faster/slower depending on environmental collapse, dissipative cooling, Bose enhancement and initial conditions.

## 4. Relevant observational theories and what they actually establish
- Bhandari et al., *Probing extended dark objects through astrometric microlensing from Gaia*, JCAP09(2026)072, https://doi.org/10.1088/1475-7516/2026/09/072. Modeled extended boson stars and Q-balls with Gaia astrometric signatures. The quoted rates/fraction constraints assume a chosen population, size/mass structure, survey selection; **do not** transfer directly to an unpredicted 86-eV Bose-star function.
- *Subaru-HSC through a different lens: Microlensing by extended dark matter structures*, PRD102 083021 (2020), https://doi.org/10.1103/PhysRevD.102.083021. Finite lens size weakens or changes point-mass searches.
- *Microlensing signatures of extended dark objects using machine learning*, PRD109 123004 (2024), https://doi.org/10.1103/PhysRevD.109.123004. Simulation/observational signature methodology.
- *Periodic Gravitational Lensing with Oscillating Boson Stars*, PRL137 031403 (2026), https://doi.org/10.1103/6w3p-8ym2. Predicts variability for **oscillating** boson-star models; not necessarily relevant to the static repulsively stabilized TF star of the source.

## 5. What is needed for real prediction
Specify from GP–kinetic evolution \(dN_{\rm BS}/dM\), spatial clustering and relative velocities; solve projected surface density \(\Sigma(R_\perp)\) of the n=1 sin profile with a first zero at R_c, compute deflection and centroid/photometric maps over lens/source distances, convolve with source angular size, sampling and noise, then integrate over a Milky Way halo distribution. With factor-π consistency and uncertainty in \(a,m\), construct an extended-lens exclusion/forecast surface in (M,R_c,f_BS). No actual microlensing likelihood has yet been evaluated.
