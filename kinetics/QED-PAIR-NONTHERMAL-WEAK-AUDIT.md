# First arbitrary-isotropic nonthermal QED e−e+↔γγ collision weak moments (2026-11-22)

## Scientific claim and boundary

The earlier full five-dimensional quantum QED e−e+↔γγ
integral assumed **thermal** Fermi electron/positron spectra and a
Bose photon spectrum, at a small set of temperatures/chemical
potentials. Here drop that restriction: accept arbitrary
isotropic bounded positive occupation functions f_e−(E), f_e+(E)
and bγ(E), and evaluate the **same on-shell quantum forward and
inverse QED 5D weak moments** without a thermal detailed-balance
identity. The thermodynamic H theorem, electric charge and total
energy are checked for genuinely nonthermal spectra.

**Still NOT solved:** mode-resolved strong collision operators
C_e−(p), C_e+(k), Cγ(k), a quantum Pauli/Bose preserving
finite-volume time integrator, physical weak ν and
charge-neutral baryon plasma, UV-valid χ1 portal, actual
physical Γel, or 86-eV production/freeze-out/observables.

Numerical gates: thermally parameterized arbitrary f,b reproduce
the independently implemented stabilized thermal FD/BE 5D
reference at Te/me=.1 with pair affinity .02 to <2×10⁻⁷
relative; when exact common thermal equilibrium is supplied,
net R/Rforward<10⁻¹⁴ and entropy/Rforward<10⁻²².
A positive/negative nonthermal electron logit perturbation
changes the net number sign accordingly and preserves
nonnegative entropy; perturbing a photon occupation bump
in opposite directions reverses net reaction but keeps
entropy positive. Total energy to relative 10⁻¹², charge
exact, photon number changes in pairs.

## Formula and physical role of weak moments

With arbitrary physical 0<f_e−,f_e+<1 and
0<bγ<∞, define (on the exact relativistic 2→2 shell)

    F = f_e−(E−) f_e+(E+)
        [1+bγ(ω₁)] [1+bγ(ω₂)],
    B = [1−f_e−(E−)] [1−f_e+(E+)]
         bγ(ω₁) bγ(ω₂),
    J = F−B.

The positive on-shell event measure dν consists of
the exact QED spin-averaged differential amplitude,
the two-spin electron/positron incoming weights, the
identical-photon 1/2, and full p,k,μ,CM θ*,φ*
phase-space; see [previous 5D derivation](QED-PAIR-FULL5D-AUDIT.md).
No phenomenological reaction rate or posthoc symmetrization.
The measured sources are

    Ndot_e−=Ndot_e+=−∫dν J,
    Ndot_γ=+2∫dν J,
    Qdot_e−=−∫dν E− J,
    Qdot_e+=−∫dν E+ J,
    Qdot_γ=+∫dν (ω₁+ω₂)J,
    σ = ∫dν (F−B)ln(F/B) ≥0.

For every event ω₁+ω₂=E−+E+; therefore
Q_total=0 exactly up to floating-point, independent
of the shape of f,b. Electric charge source
Ndot_e−−Ndot_e+=0 even though electron, positron,
photon number individually vary. The H theorem is
**not** restricted to thermal occupations: the integrand
has sign σevent≥0 for any F>0 and B>0.
The routine deliberately requires strictly interior
Fermi and strictly positive Bose occupations; degenerate
endpoints 0/1 and zero-photon/vacuum limits are covered
by the prior separate vacuum benchmark but not
numerically stable in this unrestricted log ratio yet.

## Independent thermal compatibility checks

At Te/me=Tγ/me=.1 and μe−/Te=μe+/Te=.01,
the new code numerically produces R_net=
+1.51124370419×10⁻¹⁹ me⁴ versus the previous
**independent thermal FD/BE-specific implementation**
+1.51124361519×10⁻¹⁹ at nrad24, nangle16,
nφ12 (relative difference +5.89×10⁻⁸).
Entropy from both codes agrees to within 2×10⁻⁷.
At μsum=0, the new code's direct F−B subtraction
produces a roundoff residual R_net≈1.68×10⁻³⁴,
~2.25×10⁻¹⁷ of the positive forward rate,
and tiny nonnegative σ≈1.15×10⁻⁴⁷.
This is numerical cancellation in the strictly
pointwise balanced FD/BE gas, not a false physical
chemical source.

**Negative finite-precision check:** at Te/me=1 and
μsum/Te=.02, the direct F−B vs analytically stabilized
thermal-affinity integrand differs by ~3.02×10⁻⁶
relative at nrad24/nangle16. This failure of
uniform high-T precision prevents declaring a
universal 10⁻⁷ nonthermal weak operator. Improve
event-level log ratio/scaled gain/loss and recheck.

## Nonthermal synthetic perturbations

Choose me=1, common nominal T=.1 me and
f_eq(E)=1/(exp(E/T)+1), b_eq(E)=1/(exp(E/T)−1).
The *electron* logit packet

    f_e−(E)=expit[−E/T+ε exp(−E/T)],
    f_e+(E)=f_eq(E), bγ(E)=b_eq(E)

is not FD thermal at any single common chemical potential
when ε≠0. With ε=+0.15,
R_net=+1.69053734456×10⁻²³ me⁴,
Qγ=+3.78888806086×10⁻²³ me⁵,
σ=+6.12490965689×10⁻²⁹ me⁴.
With ε=−0.15, R_net=−1.69053122009×10⁻²³
and σ=+6.12488288452×10⁻²⁹.
Charge and energy cancel eventwise. This is a
**prescribed nonthermal perturbation**, not one derived
from χ production or cosmological freezeout.

A distinct photon nonthermal perturbation

    bγ(E)=b_eq(E){1+ε exp[−((E−1.1me)/(.2me))²]}

with thermal e± gives, at nrad24/nangle20/nφ12,
ε=+.2: R_net=−1.57586871116×10⁻¹⁸,
σ=+4.01488072146×10⁻¹⁹, and
ε=−.2: R_net=+1.39592094768×10⁻¹⁸,
σ=+4.11054392776×10⁻¹⁹.
A photon occupation excess drives pair creation
but still produces entropy; a deficit drives
annihilation. A loss-only annihilation model
could not reproduce the sign reversal.

## Precision and negative controls

Electron logit bump ε=+.15, nangle20/nφ12,
incoming nrad12,16,24,32,40 yields
R_net [me⁴] =
(1.69159953405,1.69054680534,1.69053734456,
1.69053734621,1.69053734624)×10⁻²³.
At nrad12 the error relative to n40 is ~6.28×10⁻⁴;
nrad24→40 is ~1×10⁻⁹ for this smooth tiny packet.

But photon bump ε=+.2 at nrad24→32 changes
R_net from −1.57586871116×10⁻¹⁸ to
−1.57588673140×10⁻¹⁸ (relative ~1.14×10⁻⁵);
its entropy changes too. **There is no universal
grid accuracy from one perturbation, and no
proof that a discretized strong collision preserves
a nonthermal equilibrium/Pauli automatically.**

Code: [qed_pair_nonthermal_weak.py](qed_pair_nonthermal_weak.py),
[regression](../tests/test_qed_pair_nonthermal_weak.py),
[compact data](results/qed-pair-nonthermal-weak-2026-11-22.json).

    python -m pytest -q tests/test_qed_pair_nonthermal_weak.py
    python -m kinetics.qed_pair_nonthermal_weak --T-over-me 0.1 --electron-logit-bump 0.15

## Required next step

Derive the **mode-resolved strong operator** and a Pauli/Bose
positive conservative collisional time evolution for arbitrary
e±/γ occupations, with net energy and electric charge
conservation, photon pair number changes, entropy growth
and q/dt/cutoff/angle convergence tested independently
of this weak reference. Avoid interpolating distributions
or summing bins in a way that creates a false equilibrium,
as already found for χe projection. Couple to realistic ν,
charge-neutral baryon plasma, and explicit gauge/UV-valid
χe/χν mediator. The quantum-kinetic χ1↔χ2+φ
Bose/Pauli source, FRW and observations remain open;
no physical Γel, 86-eV ΩDM, macroscopic BEC,
frozen f(q) or CLASS readiness follows.
