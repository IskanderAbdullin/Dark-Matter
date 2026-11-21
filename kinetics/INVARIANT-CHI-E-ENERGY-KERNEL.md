# Directly microreversible relativistic full-Pauli χ+e 2→2 energy kernel — NO post-hoc symmetry

**Date:** 2026-11-21. **Code:** [microreversible_chie_lattice.py](microreversible_chie_lattice.py), [independent unit tests](../tests/test_microreversible_chie_lattice.py), [runner](verify_microreversible_lattice.py), [compact measured scan](results/microreversible-chi-e-qkernel-2026-11-21.json).

This is an explicit answer to the earlier **negative finite-bin reversibility result**, not a physical 86-eV Bose dark matter relic, BEC, CLASS-ready calculation or a full thermal SM portal. We work with a Dirac χ of mass M=.1 GeV and massless e±, vector contact G=.001 GeV^-2, gχ=ge=4, physical T=1 GeV. Natural units ħ=c=kB=1. The independent continuum source is the **separate** old five-dimensional p,k,cosθ,cosθCM,φCM Gauss kernel, not our new q-parametrization.

## Derivation: exact invariant energy-differential 2→2 scattering kernel

For an isotropic fermionic χ and an electron/positron e, the integrated, spin-averaged, Lorentz-invariant collision event measure is

```text
dR = gχ ge dΠχ dΠe dΠχ' dΠe' (2π)^4 δ⁴(P+K−P'−K') |M(s,t)|²
         × fχ(E) fe(k) [1−fχ(E')] [1−fe(k')],

dΠ=d³p/[(2π)³2E], with |M|² averaged over the two incoming spins per particle.
```

Use energy variables E,k,E',k', three independent energy integrals after δ(E+k−E'−k'); physical momentum magnitudes p=√(E²−M²), p'=√(E'²−M²), and massless electron momenta k,k'. Introduce the physical momentum-transfer magnitude \(\ell=|\mathbf p-\mathbf p'|=|\mathbf k'-\mathbf k|\). At fixed positive energies the allowed interval is

```text
ℓ_min = max(|p−p'|, |k−k'|),
ℓ_max = min(p+p', k+k').
No physical 2→2 event if ℓ_max<=ℓ_min or E'<=M or k'<=0.
```

At a given ℓ, incoming χ/e cone cosines with respect to the transfer vector are

```text
cχ=(p²+ℓ²−p'²)/(2pℓ),
ce=(k'²−k²−ℓ²)/(2kℓ).
D0=2(Ek−pk cχ ce), t=(E−E')²−ℓ².
<D²>ψ = D0² + 2p²k²(1−cχ²)(1−ce²),
with D=s−M²=2(Ek−p·k), average over relative cone azimuth ψ.
```

The exact vector-contact spin-average

```text
|M|²=2G²[(s−M²)²+(u−M²)²+2M²t], u−M²=−(D+t)
```

has **analytic** cone-azimuth average

```text
<|M|²>ψ = 2G²[ 2<D²>ψ + 2 D0 t + t² + 2 M² t ].
```

Integrating the two independent spherical shell angular convolutions \(\int dΩ_p dΩ_{p'}δ³(\mathbf q-\mathbf p+\mathbf p')=2π/(p p' \ell)\) and likewise for e, with the invariant four-body measure, gives the positive energy-differential rate

```text
𝒦(E,k;E',k') = gχ ge/(256 π^5) ×
                ∫_(ℓmin)^(ℓmax) dℓ <|M(s,t)|²>ψ
and k'=E+k−E'.
Units: 𝒦 GeV, after dE dk dE' -> GeV^4 event-number rate per volume.
```

This expression is **time-reversal symmetric before any numerical symmetrization** because the invariant event measure and squared amplitude are unchanged by incoming↔outgoing exchange. The invariant q interval is identical after exchange (although χ/e cone labels change). A *second independent direct 2048-point azimuth average* of the raw vector amplitude at M=.1,E=1.1,k=2,E'=1.7,k'=1.4,ℓ=1.5 gives the same 8.564963736×10^-5 dimensionless squared amplitude as the analytic cone average, to floating precision. Four distinct physical incoming/outgoing energy configurations also pass forward/reverse \(\mathcal K\) equality to relative <10^-12, without a W=(Wraw+Wraw^T)/2 operation. Both tests are necessary but not sufficient for complete thermal plasma physics.

## A paired common-energy lattice without ad hoc W arithmetic averaging

Take χ energy centres \(E_i=M+(i+1/2)h\) and e energy centres \(k_a=(a+1/2)h\). For a given transition \(i,a\to j,b\), set \(b=i+a-j\), so **all four energy centres conserve total energy exactly**. Define an unordered positive quartet rate

```text
W_(ia,jb) = 𝒦(Ei,ka;Ej,kb) h³, i<j, b=i+a-j, and W_(jb,ia)=W_(ia,jb)
by the PHYSICAL kernel's proven reverse identity, NOT by post-hoc averaging.
```

The cell Fermi capacities for time evolution are exact within the chosen energy-bin shell volume, \(Vχ_i=gχ(p_{i+1}³−p_i³)/(6π²)\) and \(Ve_a=ge(k_{a+1}³−k_a³)/(6π²)\). The discrete number flow is \(J=W[F−R]\), with both χ/e initial AND final Pauli blocking. A direct sparse unordered event graph updates χ and e charges separately, and each quartet has \(E_i+k_a=E_j+k_b\). Thus \(\dot Nχ=\dot N_e=\dot Eχ+\dot E_e=0\) (for the sum-of-sources invariants), while \(\dot S=\sum W(F−R)\ln(F/R)\ge0\) for interior fermionic occupations. At common T and arbitrary independent constant χ/e chemical potentials, microscopic F=R for **each** quartet before summing. That is the crucial check which failed for the old bare interpolated 5D W.

The same sparse quartet coefficients feed a fully implicit four-fermion scalar-root update conserving both numbers and total χ/e energy with physical Pauli bounds. Sequential full-network time evolution is a fixed-volume quantum kinetic **energy-collocation** scheme; unlike previous work, there is *no physically artificial post-hoc repair* to W. Yet this discretization approximates the continuum energy integrals and does not automatically provide a fully physical Γ(T), cosmological g* or correct late relic.

## Independent numerical normalisation, convergence and important nonmonotonicity

The external-to-this-algorithm numerical reference is the earlier directly evaluated **five-dimensional CM scattering** (not q-transfer) with Gauss n=28 in *each* p,k,µ,cosθCM,φCM dimension, p/k cutoff 22 GeV and physical masses/couplings as above. It gives

```text
Qchi(Tchi=.9 GeV,Te=1 GeV)  = +1.9851492233308662e-8 GeV^5
Qchi(Tchi=1.1 GeV,Te=1 GeV) = -4.432408224385169e-8 GeV^5.
```

The new **independent one-dimensional ℓ integral kernel** is evaluated at nℓ=12 Gauss points (no CM-angle Monte Carlo) and the discrete energy convolution is evaluated on h-spaced χ/e lattice centres from 0 to ~22 GeV. Its paired Pauli FD two-temperature Qχ is:

| h (GeV), χ/e nodes | Qchi(Tχ=.9), GeV^5 | Qchi(Tχ=1.1), GeV^5 | Relative error to independent 5D reference (cold/hot) |
|---|---:|---:|---:|
| 1.0, 22/22 | +2.0066428148856442e-8 | −4.475030225116334e-8 | +1.083% / +0.962% |
| 0.5, 44/44 | +1.9854909200454386e-8 | −4.433145477989673e-8 | +0.0172% / +0.0166% |
| 0.25, 88/88 | +1.9833906119488718e-8 | −4.428968123117066e-8 | −0.0886% / −0.0776% |
| 0.125, 176/176 | +1.9851092070921427e-8 | −4.432106079842345e-8 | −0.00202% / −0.00682% |

**Important numerical negative control:** The h-convergence is NOT monotonic; coarsening from h=.5 to h=.25 increases the error magnitude before h=.125 improves it to <0.01%. This is a finite-energy-collocation/endpoint error, not a physical nonthermal reaction. Conserved charges and entropy are insufficient to infer q/cell accuracy. The nℓ=8 to 12 comparison at h=.5 shifts Qχ cold only ~1.45×10^-5 relative, indicating small transfer-quadrature error there; remaining variation is dominated by energy-cell representation/cutoffs. The last fine h=.125 values were measured independently in the sandbox and recorded in compact JSON; default CI tests use h=.5/.25 for size.

No finite scalar kinetic gamma is inferred from these numbers. A direct physical χ/e Pauli source for both cold and hot relative to SM, a discrete thermal FD fixed point and positive entropy are demonstrated for this **conditional toy EFT**. Massless electrons at T=1GeV are a simplifying approximation, while the true gauge-lepton UV portal and neutrino channels must be added; the contact EFT needs validity checking over the thermal momentum range.

## What remains falsifiable

1. The exact **continuum** Boltzmann collision operator is replaced by an energy-collocation quadrature: physically accurate only after q/energy cutoff, M/T, interpolation, dt and in-medium species convergence. The published 1D kernel is exact for the selected tree amplitude but its finite energy lattice is not.
2. The full χ/e network has Lie splitting and can show order/dt errors despite exact Pauli/N/E/entropy. Independently test time history and compare to 5D nonlinear strong Cχ(p), including IR/high-energy localized modes.
3. Real SM early-universe physics requires finite electron mass and e± dilution around MeV, gauge-required νe_L channel, in-medium corrections, mediator completion and number-changing channels, and realistic g*ρ/g*S(T).
4. Combine physically validated χ1 elastic scattering with Yukawa χ1↔χ2+φ in pinned BEST; independently converge the large-timescale Bose-stimulated IR and potential zero-mode physics before any physical 86eV ΩDM, freezeout, CMB/Lyman-α, or CLASS-ready statement.

Until those gates pass, **physical 86-eV dark matter is NOT validated and the draft PR must stay draft**.
