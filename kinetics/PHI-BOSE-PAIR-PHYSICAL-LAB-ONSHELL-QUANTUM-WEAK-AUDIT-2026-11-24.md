# First physical real-mass 5D lab-frame 2-to-2 quantum χ2+antiχ2↔φ+φ weak source (2026-11-24)

## Derivation and scope

The previous invariant Dirac Yukawa χ2+antiχ2→2 real φ tree cross section and exact Maxwell-Boltzmann equilibrium inverse ratio provide the amplitude, but not a fully quantum physical phase-space operator. Here both initial physical three-momenta, their relative angle and the CM final scattering polar angle/azimuth are explicitly integrated. The common positive event rate per physical volume is:

  dW = [gχ*ganti/(8π^4)] p1²dp1 p2²dp2 dcosθ_initial
       × v_Moller × (dσforward/dcosθ_CM) dcosθ_CM dψ_CM/(2π).

For two separate Dirac charge sectors, gχ=ganti=2. The spin-averaged vacuum dσforward/dcos already incorporates 1/2! for two identical real final scalars, so NO second identical factor is inserted in the common F/R kernel. The pair COM is Lorentz boosted to the laboratory plasma frame; each event carries exact Echi+Eanti=Ephi3+Ephi4, with massive on-shell daughter Ephi≥mphi. The outgoing CM azimuth matters for nonthermal isotropic LAB energy occupations.

Quantum forward and reverse factors use occupations at the ACTUAL FOUR LAB energies:

  F=fχ(E1) fanti(E2) [1+bφ(E3)][1+bφ(E4)]
  R=[1-fχ(E1)][1-fanti(E2)] bφ(E3)bφ(E4)
  A=ln(F/R)=ηφ(E3)+ηφ(E4)−ηχ(E1)−ηanti(E2)
  dH=dW(F−R)A≥0.

A single positive W for both directions prevents artificial thermal-equilibrium flux. Moment charges are dNχ=dNanti=−dR_events, dNφ=+2dR_events, so true net vector Dirac U1 is unchanged; energy Qχ=−E1 dR_events, Qanti=−E2 dR_events, Qphi=(E3+E4)dR_events and sums identically to zero for each event.

This is a physically normalized **instantaneous continuum-event quadrature / weak source** for a finite physical cutoff P, not a conservative full 5-field finite-time entropy-dual discretization, not an in-medium scattering amplitude, and not any real FRW relic or Bose condensate calculation.

## Independent rate normalization and negative grid results

At common T1GeV with classical fχ=fanti=exp(−E/T), no Bose/Pauli, the new LAB p,q,relative-angle,CM-polar/azimuth integration independently reproduces the invariant Gondolo-Gelmini thermal forward rate from the previous CM/Bessel integration:

  gamma_MB,GG(T=1GeV)=2.0886349101080e−16 GeV^4.

The lab angular and radial integrations show meaningful convergence:
| (nrad,nrelative,nCM,naz) | LAB gamma_MB(1GeV) GeV4 | Relative to independent GG |
|---|---:|---:|
| (12,12,32,4) | 2.10735030244e−16 | +0.8961% |
| (16,12,48,8) | 2.08976506317e−16 | +0.05411% |
| (20,16,64,8) | 2.08901536757e−16 | +0.01822% |
| (24,24,64,6) | 2.08880064525e−16 | +0.00794% |
| (24,32,96,1) | 2.08865651360e−16 | +0.001034% |

The nCM32 coarse grid is a material ~0.9% overestimate, despite exactly sharing forward/reverse event weights. Thus preservation of H/energy/detailed balance alone does NOT establish physical rate convergence; more realistic high-T and hard-cutoff physics must pass independent three-grid controls. The parent cutoff P=20GeV at T=1 excludes extremely cold φ near the actual 86eV mass, so no BEC/zero-mode conclusion is possible.

## Quantum FD/BE physical weak source at finite cutoff

For P20GeV, nrad16,nrelative12,nCM48,naz8 with initial occupations FDχ(1GeV,μ=.05), FDanti(1GeV,μ=−.05), BEφ(1GeV,μ=−.2), the direct+inverse quantum rate is:

  dR_events=+7.54787689038e−17 GeV4,
  quantum H=+3.01915075615e−17 GeV4,
  exact H/dR=0.4 (dimensionless chemical affinity).

Both physical number legs, U1 and energy ledger agree exactly at event level and energy residual is ≲1e−30 GeV5. A curved nonthermal trial FDχ(T=.9),FDanti(T=1.1),BEφ(T=1) with positive curvature ηφ+=.01*(1−cos E) yields dR_events≈+6.10894e−19 and H≈+1.20570e−17 GeV4. Common-T μχ=−μanti=.05 and μφ=0 is stationary to floating precision, in accord with the full real-scalar Yukawa network chemical affinities. These are instantaneous rates, not time-integrated fφ(q).

## Reproduction and missing physics

Source: kinetics/phi_bose_pair_physical_onshell_quantum.py; five focused tests tests/test_phi_bose_pair_physical_onshell_quantum.py. The code explicitly writes positive event W and all four massive E, so a third-party check can recalculate moments for arbitrary nonthermal η. The pre-existing independent CM/Bessel reference is kinetics/phi_bose_pair_independent_vectorized.py; the previous chemical-odds source is kinetics/phi_bose_pair_chemical_affinity.py.

Commands: python -m pytest -q tests/test_phi_bose_pair_physical_onshell_quantum.py ; python -m kinetics.phi_bose_pair_physical_onshell_quantum.

NEXT_ACTION: extend this exact common 2↔2 event measure into a positive/collision-invariant-preserving fully implicit finite-time five-field Dirac/anti quantum moment scheme, **together with both decays and antiparticle decays**. At the same time calculate in-medium χ and φ dispersion, Γ(T) and SM portal. Check parent radial/relative/CM/azimuth/IR/UV time convergence before physical freezeout/Ω/CLASS. No actual cosmological/galactic prediction follows yet.
