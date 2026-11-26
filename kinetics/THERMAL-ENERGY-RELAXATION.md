# Number-conserving FD thermal ENERGY-moment rate (not an elastic BGK or kinetic eigenvalue)

**Status:** independent numerical/analytic benchmark for the **restricted FD trial manifold** of χ1 scattering on a thermal massless e± bath. See [thermal_energy_relaxation.py](thermal_energy_relaxation.py), [unit tests](../tests/test_thermal_energy_relaxation.py), [real 5D microscopic collision model](relativistic_weak_form_2to2.py), and [microscopic portal/EW chiral-channel caveats](MICROPHYSICS-MODEL.md). It is not the actual full momentum-dependent γ_el(p,T), not a late freezeout or physical 86 eV relic.

## Exact thermodynamics at fixed χ number

Elastic χe contact preserves χ number. For an imposed FD(Tχ,μχ), demand `Nχ(Tχ,μχ)=Nχ(Te,μ0)` rather than incorrectly setting μχ=0 after heating. Using `E=√(p²+M²)`, `f=FD((E−μ)/T)`, `z=f(1−f)`, the six exact thermodynamic integrals are:

```text
Nχ = gχ/(2π²) ∫dp p² f
ρχ = gχ/(2π²) ∫dp p² E f
∂N/∂T|μ = gχ/(2π²) ∫dp p² z(E−μ)/T²
∂N/∂μ|T = gχ/(2π²) ∫dp p² z/T
∂ρ/∂T|μ = gχ/(2π²) ∫dp p² E z(E−μ)/T²
∂ρ/∂μ|T = gχ/(2π²) ∫dp p² E z/T
dμ/dT|N = −N_T/N_μ
c_N = dρ/dT|N = ρ_T − ρ_μ N_T/N_μ.
```

Number matching is performed with Brent, independently checking `N(T±ΔT,μ(T±ΔT))=N0`, and the analytical derivatives/c_N are cross-checked against a second finite-difference calculation using this root. These are thermal equation-of-state statements, not an assumed external reservoir reset.

## Kinetic energy exchange and a single conditional rate

The physical χe contact kernel provides the fully Pauli-blocked 5D instantaneous **weak energy transfer** Qχ(Tχ,μχ;Te). Define the linear-response energy conductance at *fixed χ number*:

```text
K_N = − (d Qχ/d Tχ)_{Nχ, Tχ=Te}  [GeV^4]
γ_T,moment = K_N / c_N                 [GeV].
```

Qχ>0 for a slightly colder χ FD state, Qχ<0 for a hotter one; at common T it is zero. This is an effective cooling rate for the **one-dimensional thermal-energy FD ansatz**, not a full eigenvalue of the linearized momentum-resolved Boltzmann operator. The true collision has infinitely many relaxation modes; a *single* BGK γ_el(p,T) cannot be inferred merely from this number. Any statement that the previously phenomenological γ_el=10^5 was microscopically calibrated is **false**.

**Sandbox numerical reference, not an external BEST run:** Mχ=.1 GeV, Te=1 GeV, μ0=0, gχ=ge=4, G=.001 GeV⁻², massless electrons, no neutrinos, thermal screening or SM heat evolution, `cutoff=32`. Number-matched c_N=0.9987065946 GeV³. From ΔT=.04,n=16 to ΔT=.005,n=28, `γ_T,moment` approaches ≈3.03183×10⁻⁷ GeV (ΔT=.01,n=24 gives 3.03184731×10⁻⁷). At ΔT=.01,n=24, `K_N=3.02792590×10⁻⁷ GeV⁴`; the number constraint is preserved numerically at the ≲10⁻¹² relative level. The numerical refinement of **both** ΔT and 5D mesh is shown in code/tests; do not mistake it for a full spectral or observational uncertainty budget.

For the **specific** equal-temperature massless SM lepton-doublet completion with G_L=G_R=G, the same χ energy moment adds the νe channel and scales by **6/4=1.5**, but this still presumes the same collision geometry per chiral lepton, bath temperature and neglects all realistic finite-me, neutrino, thermal/UV effects. Testing ge=6 is a *conditional model limit*, not automatically valid at arbitrary epochs.

Next independent physics gate: compute the **linearized nonthermal collision operator/eigenspectrum** from the generic Cχ(p), compare its slowest non-number mode(s) and heat mode against this FD-manifold rate, then integrate time dependence with a conservative, Pauli-preserving finite-volume scheme and finite bath. Do not assert `CLASS-ready`, Ωφ=ΩDM, BEC, or freezeout.
