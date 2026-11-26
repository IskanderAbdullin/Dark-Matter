# Microscopic finite-temperature moving Z′↔νν Pauli/Bose source: independent quantum inverse-production test (2026-11-22)

## Physical scope, hypothesis and predeclared tests

After the exact real-φ U(1)′ charge audit and daughter χ2–SM scattering comparison, a light Z′ is an **actual dynamical particle**, not merely a virtual elastic propagator. The prior reference quantified its hypothetical equilibrium energy and a Maxwell–Boltzmann neutrino decay clock. Neither calculation included full neutrino Pauli blocking, inverse ννbar→Z′ or Z′ Bose enhancement.

Here derive the **complete 1↔2 gain–loss source for one massive, unpolarized Z′ momentum mode**, interacting with two *externally thermalized* active neutrino flavor-antiflavor pairs of common temperature Tν and zero chemical potentials. The microscopic matrix element is normalized to the standard tree partial width Γ0=g′²mV/(12π). The on-shell two-body decay angles are explicitly integrated; the final formula is a **real microphysical momentum-resolved collision operator for Z′ in the restricted fixed-FD bath**, NOT a phenomenological χ elastic BGK Γel and NOT a self-consistent evolving finite SM reservoir.

Predeclared numerical gates: independent Gauss CM-angle average vs analytic logarithmic FD formula rel 2e−12 for mV/T=.03,.1,.3 and p/T∈[0,5]; p=0 analytic limit; exact Bose fixed point, correct gain-loss direction and entropy, energy and net ν-number ledgers; radial source n128→256 error <1e−6 in benchmark. Fail closed on any discrepancy; absence of finite ν dynamical backreaction and initial Z′ abundance must remain explicit.

## Exact derivation for moving-vector decay

For E_V=√(p²+mV²), β=p/E_V, unpolarized three-polarization Z′ decay is isotropic in its rest frame. Each massless final ν and antineutrino has thermal-frame energies

    ε1,2=(E_V ± p cosθ*)/2,   ε1+ε2=E_V.

With f1=f_FD(ε1/Tν), f2=f_FD(ε2/Tν) and arbitrary boson occupation bV(p)≥0, local quantum loss/gain are

    F=bV(1−f1)(1−f2),    R=(1+bV) f1 f2.

At common Tν and μν+μνbar=0,

    F/R = [bV/(1+bV)] exp(E_V/Tν),

so exactly F=R if bV=bBE(E_V/Tν), including for a moving vector.

The mode-resolved source is

    C_V(p) = −Γ0 (mV/E_V) <F−R>_(cosθ*)
           = Γ0 (mV/E_V) A(p,Tν)[bBE(E_V/Tν)−bV(p)],

where the **physical** FD blocking-minus-inverse coefficient is

    A=< (1−f1)(1−f2)−f1f2 >
     =<1−f1−f2>.

This is not an arbitrary relaxation constant: the coefficient comes directly from the neutrino quantum collision and two-body kinematics. Its analytic integral is

    A=1−(2Tν/p) ln[(1+exp(−(E_V−p)/(2Tν)))/
                       (1+exp(−(E_V+p)/(2Tν)))],  p>0,

    A(p→0)=tanh(mV/(4Tν)).

Thus the physical mode relaxation ΓV(p)=Γ0 (mV/E_V) A(p,Tν) is momentum dependent, suppressed by both time dilation and FD Pauli/inverse cancellation near mV/Tν≪1. **It is a neutrino↔vector 1↔2 rate, NOT the χ1/χ2 elastic Γel or its spectrum**.

An entirely independent 96-node Gauss angle integral of 1−f1−f2 reproduces analytic A within 2e−12 relative across the stated six cases. At mV/Tν=.03, A(p=0)=0.0074998593782; at p/Tν=1, A=0.24033296419; at p/Tν=3, A=0.57031607510. These nonconstant values invalidate copying the earlier MB K1/K2 as an exact FD/BE chemical rate.

## Number/energy and total quantum entropy

For massive vector spin degeneracy gV=3,

    NdotV = 3/(2π²) ∫dp p² C_V(p),
    QdotV = 3/(2π²) ∫dp p² E_V C_V(p),

    Ndot(ν+νbar) = −2 NdotV,
    Qdot_ν = −QdotV,

where the ν source counts both particles per physical 1↔2 event (two active flavors are included in Γ0). Net individual lepton-flavor asymmetries remain zero for zero ν chemical potentials. The energy identity is exact on each shell, with no free energy bath; **neutrino energy is not yet integrated dynamically**, so Qν is an instantaneous required ledger.

For bV>0, total instantaneous quantum entropy production is

    Sdot_total = 3/(2π²) ∫dp p² Γ0 mV/E_V
          <(F−R) ln(F/R)> >=0.

Independently compute SdotV=−∫C_V ln[bV/(1+bV)]
and Sdotν=−QdotV/Tν, and confirm Sdot_total=SdotV+Sdotν.
For bV=0, the vector entropy derivative has a logarithmic boundary
singularity; the runner reports **entropy=null**, not a fictitious finite number.

For a controlled family bV(p)=r bBE(E_V/Tν):
r=0 (empty), r=.5 (deficit) implies C_V>0 and inverse
production; r=1 has C_V=0 exactly; r=2 implies C_V<0
and decay, with nonnegative total entropy for both interior
nonthermal shapes.

## Direct finite-T inverse production and expansion-scale diagnostic

Use the **same arbitrarily chosen** zero-transfer contact coefficient
G m_e²=.001 and the convention |gχqχ|=√(4π), so
g′=G mV²/√(4π), for the three masses mV/me=.03,.1,.3.
Let m_e=0.51099895 MeV, Tν=m_e, assume a radiation reference
H=1.1640516263e−25 GeV at gstarρ=10.75, with no dark-sector
backreaction. With *initially empty Z′* and FD thermal ν,
the exact Pauli/Bose microscopic inverse production yields

    R_fill = NdotV(bV=0)/nV_BE(Tν),

where nV_BE is the **conditional thermal equilibrium** number
density of a three-polarization massive vector, used only as a
normalizing target. This is not an assumed initial Z′ density.

| mV/me | g′ from the exploratory contact G | R_fill/H from quantum inverse source | Previous unblocked MB Γ/H |
|---:|---:|---:|---:|
| .03 | 2.53885313e−7 | **1.34973037e3** | 3.37206100e3 |
| .10 | 2.82094792e−6 | **1.84903984e6** | 4.57678983e6 |
| .30 | 2.53885313e−5 | **1.33809254e9** | 3.16439783e9 |

This *new quantum inverse result* is approximately 40–42%
of the older MB time-dilated decay clock, showing explicitly why
a vacuum width or unblocked rest scattering frequency is not a
precise finite-T production rate.

At mV=.03me, n128→n256 changes R_fill/H by
about 7.36e−9 relative; n256→512 by <4e−13.
For mV=.1,.3, the converged n256/n512 differences
are at roundoff. These are 1D radial convergence tests
**for a fixed thermal ν reservoir**, not a full cosmological
time-integration, not an actual mediator abundance or ΔNeff.

The corresponding quantum energy-transfer ledger for mV=.03me
at bV=0 is QdotV=+9.28950708536e−36 GeV^5,
Qdotν=−QdotV; in r=.5 case QdotV halves,
Sdot_total=+2.23340381e−33 GeV^4. At r=2,
QdotV reverses and Sdot_total=+4.03273497e−33.
At r=1, all number, energy and entropy collision moments
vanish to the implemented exact detailed balance.

**Scoped negative physical result:** in this particular
illustrative gauge-coupling convention, even the correct full
moving-vector quantum inverse ν source is extremely rapid
compared with the assumed Hubble clock. Neglecting a real Z′
population is therefore **not justified by the
virtual-exchange χ–lepton amplitude alone**. This is
*not* a terminal exclusion or proof of cosmological
thermalization: once produced, the mediator transfers
heat and particle number back to neutrinos, changing
Tν, expansion and inverse rates; other χ1/χ2 channels,
finite-temperature masses, constraints, gauge symmetry
restoration and actual initial conditions must be evolved.

## Reproducibility and required next physical gate

[Executable source](vector_neutrino_quantum_1to2.py),
[10 regression tests](../tests/test_vector_neutrino_quantum_1to2.py),
[compact results](results/vector-neutrino-quantum-1to2-2026-11-22.json).

    python -m pytest -q tests/test_vector_neutrino_quantum_1to2.py
    python -m kinetics.vector_neutrino_quantum_1to2 --mV-over-me 0.03

Physical context: [light gauge-boson impact on neutrino
decoupling and N_eff, PRD 110, 075032 (2024)]
(https://journals.aps.org/prd/abstract/10.1103/PhysRevD.110.075032).
The microphysical shell source here is independently derived,
not a published numerical cosmology.

Next P0: **remove the fixed ν thermostat** by evolving the
neutrino/antineutrino phase-space or at least a finite
energy/entropy-conserving ν bath with self-consistent Tν
and number/chemical charges; integrate actual Z′
production and decay with H and ν/e±γ coupling.
Add open Z′↔χ2χbar2 and χ1χbar1↔Z′Z′
processes, gauge-breaking S and light mediator bounds.
Then couple **both** χ1/χ2–SM full relativistic
Pauli 2→2 operators and χ1↔χ2+φ Bose
1↔2, establish full q/dt/cutoff/time
convergence and a late frozen 86-eV spectrum.
No physical χ Γel, ΩDM, BEC or CLASS-ready
model has yet been demonstrated.
