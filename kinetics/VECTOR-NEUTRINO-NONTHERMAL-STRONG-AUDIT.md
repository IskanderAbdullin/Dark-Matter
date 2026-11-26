# Arbitrary nonthermal ν strong mode source vs independent Z′ weak collision (2026-11-22)

## Falsifiable objective and model scope

Extend the physically normalized moving-vector Z′↔ννbar quantum Pauli/Bose 1↔2 source beyond an assumed thermal FD neutrino bath. Construct a **mode-resolved strong neutrino occupation derivative** Cν(k) for prescribed arbitrary isotropic 0<fν(k)=fνbar(k)<1 and positive bosonic bV(p), then independently reproduce the *full weak* vector number and energy source and total entropy from two different integration coordinates.

This is **not** the full evolving neutrino distribution, Hubble/weak SM kinetics, χ1/χ2 elastic Γel or 86-eV φ relic. The vector unpolarized rest decay is isotropic; only the two active neutrino flavors and antiflavors are modeled, with gν,total=4, gV=3 and vacuum combined Γν0=g′²mV/(12π). The physical interacting neutrino/antineutrino distributions are **assumed equal**; extension to distinct µ/asymmetries and finite ν self-scattering is still required.

Predeclared tests: full-weak p+angle vs strong k+p integrated N,E and entropy relative below 1e−9 for nonthermal sign and amplitude cases at mV=.3, n64; exact kinematic inverse mapping; and explicit *negative* unsplit-k quadrature producing a biased number source even at moderately large order. Expand mV and cutoff scans, preserving failure cases.

## Derivation and exact prefactor

For a vector momentum p, E=√(p²+mV²), an on-shell decay has neutrino daughter energy k=(E+p cosθ*)/2, hence dk=(p/2)dcosθ*, with allowed k∈[(E−p)/2,(E+p)/2] and dP/dk=1/p for an unpolarized vector. Conversely for fixed k>0, the parent integration starts at

    p_min(k) = | k − mV²/(4k) |.

This has a **cusp at k=mV/2** even though the underlying phase space is smooth in p,cosθ*. A naive single Gauss rule across this cusp can generate an artificial strong–weak mismatch.

At any on-shell event define

    F = bV(p) [1−fν(k)] [1−fνbar(E−k)],
    R = [1+bV(p)] fν(k) fνbar(E−k),
    J = F−R.

Vector weak source and energy:

    NdotV = −gV/(2π²) ∫dp p² Γν0 mV/E <J>θ,
    QdotV = −gV/(2π²) ∫dp p² E Γν0 mV/E <J>θ.

For equal fν=fνbar, the **strong** occupation source per one helicity mode is

    Cν(k) = [2gV/gν] Γν0/k²
            ∫_{p_min(k)}^{pmax} dp (p mV/E) J(p,k),

with gV=3 and gν=4: the factor 2 explicitly counts BOTH massless daughter particles emitted per event and the g ratios. This is not a BGK relaxation ansatz: it is the Lorentz-invariant two-body differential decay measure with full nonthermal Bose/Pauli weights. It satisfies

    gν/(2π²)∫dk k² Cν(k) = −2 NdotV,
    gν/(2π²)∫dk k³ Cν(k) = −QdotV.

The entropy identity is independently audited:
Sdotν = gν/(2π²)∫dk k² Cν log[(1−fν)/fν],
SdotV = gV/(2π²)∫dp p² C_V log[(1+bV)/bV],
Sdotν+SdotV = ∫dν Γν0 (F−R) log(F/R) ≥0.
This tests the quantum H theorem under nonthermal occupations, not only at a thermal fixed point. Equality holds only at detailed balance in all sampled modes. There is no external thermal reset in these instantaneous collision moments.

## Independent numerical checks

Use units mV=.3, momentum cutoff pmax=12, Γν0=1 for rate normalization, and prescribe

    fν(k)=fνbar(k)=expit[−k+ε exp(−k)],
    bV(p)=r/[exp(√(p²+mV²))−1],

with (ε,r)=(+.15,.75),(−.15,1.25),(+.15,1),(0,.75). These are prescribed test shapes, not solutions of the cosmological model.

At (ε,r)=(+.15,.75), Gauss n64 for vector p, its rest angle, neutrino strong p and **each side** of the k=mV/2 split gives:

    NdotV=+8.45748524863e−3,
    −2NdotV=−1.69149704973e−2,
    Sdot_total=+3.14257134000e−3,

in the chosen dimensionless energy units. Independent strong vs weak relative discrepancies at n64 are ≲8×10⁻¹⁴ for number, energy and total entropy. The three other nonthermal controls likewise reproduce N, E and entropy below 1e−9 relative; signs of integrated flow vary with the test occupations, while entropy stays positive.

At mV=0.03, the forward/pair inverse kernels develop a stronger scale hierarchy, and strong/weak consistency converges **more slowly**:

| mV | Gauss n32 relative strong N error | n48 | n64 | n96 |
|---:|---:|---:|---:|---:|
| 0.03 | +7.553e−4 | −1.866e−5 | −2.750e−7 | −5.464e−8 |
| 0.30 | −2.500e−7 | −5.390e−11 | ~7e−14 | ~3e−14 |
| 1.00 | ~−7e−12 | ~2e−14 | ~−8e−15 | ~−2e−14 |

**Negative quadrature result:** at mV=.3, using a *single unsplit* 32-node Gauss k integration over the entire support gives relative ν number mismatch **−3.1083×10⁻4**, whereas splitting at k=mV/2 (32 nodes per side) reduces that mismatch below 1e−6; at n64 per side it is below 1e−12. A false source imbalance can therefore emerge even from the correct analytic collision equation solely because the strong integration ignores the support cusp.

The **physical cutoff remains independently uncertain**. Even when weak and strong agree at the same cutoff, raising pmax from 8 to 12 to 16 changes the actual integrated number rate (mV=.3, same ε,r): −0.01685285268, −0.01691497050, −0.01691665823, respectively. Agreement between two discretizations at a shared cutoff does NOT certify the true continuum rate or a physical neutrino thermalization time.

## Reproducibility and remaining gates

Source: [vector_neutrino_nonthermal_strong.py](vector_neutrino_nonthermal_strong.py);
tests: [test_vector_neutrino_nonthermal_strong.py](../tests/test_vector_neutrino_nonthermal_strong.py);
compact measurements: [vector-neutrino-nonthermal-strong-2026-11-22.json](results/vector-neutrino-nonthermal-strong-2026-11-22.json).

    python -m pytest -q tests/test_vector_neutrino_nonthermal_strong.py
    python -m kinetics.vector_neutrino_nonthermal_strong

The tests use TWO independent integration parameterizations, not an arbitrary symmetric rate closure. Next integrate **actual** bV(p,t), fν(k,t), fνbar(k,t) with positivity/fermion Pauli constraints and exact finite-volume number/energy/entropy, compare time and grid errors, and include physical ν elastic and weak e±γ reactions, χ1/χ2 vector channels, original χ1↔χ2+φ Bose source, mediator mass/UV bounds and H. No physical χ Γel, 86-eV ΩDM, BEC, frozen fφ or CLASS result.
