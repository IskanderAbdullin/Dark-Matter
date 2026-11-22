# First finite-volume nonthermal ν, antiν and Z′ event-by-event quantum time evolution (2026-11-22)

## Falsifiable numerical question and assumptions

The earlier exact moving-vector Z′↔ννbar 1↔2 source established
independent **strong ν** versus **weak vector** number, energy and entropy
moments, but **no nonthermal time evolution**. A reduced finite ν
Fermi-Dirac reservoir time pilot assumed instantaneous neutrino
thermal re-equilibration. This experiment tests a stricter
alternative: independently propagate **distinct arbitrary
nonthermal** fν(k), fantiν(k) and bV(p), with NO thermal reset.

Predeclared gates: microscopic eventwise positive Bose/Fermi
bounds; exact discrete conserved reduced charge
Nν+Nantiν+2 NV, lepton asymmetry Nν−Nantiν and energy; nondecreasing
quantum entropy; exact chemical equilibrium fixed point μV=μν+μantiν;
demonstrate finite-grid cutoff and dt errors against a separately
integrated continuous angular weak source. This is a **numerical
conservative lattice pilot**, not a continuum-converged strong
operator or physical cosmological Γ/H.

Reference framework:
[Froustey, Pitrou & Volpe 2020, full neutrino quantum kinetics]
(https://arxiv.org/abs/2008.01074) and
[Discrete Quantum Kinetic Equation (2023)]
(https://doi.org/10.1007/s44007-023-00070-1).
Our implementation is independently derived, no source code copied;
licenses of third-party source implementations do not apply.

## Exact lattice shell and on-shell event weights

Choose natural units, unpolarized massive vector with gV=3,
massless two-active-flavor neutrinos with gν=2, gantiν=2.
Use a common energy lattice of spacing h:

    E_V,j=j h,   p_V,j=sqrt((j h)^2−mV^2),
    kν,i=i h,   kantiν,l=l h,   j=i+l.

Require integer mV/h and Emax/h; vector j>mV/h excludes its
zero-measure p=0 endpoint. Exact daughter kinematics requires

    |kν,i−kantiν,l| <= p_V,j.

The discretized positive physical phase-space weights are

    wV,j=3/(2π²) p_V,j E_V,j h,
    wν,i=wantiν,i=2/(2π²) k_i² h,
    W_jil=wV,j Γ0(mV/E_V,j) h/p_V,j >0,

from the exact isotropic moving-vector rest decay
d cosθ*/2=dkν/p. Here Γ0=1 is an **arbitrary unit of time**,
not a viable measured coupling. Each allowed j,i,l denotes a
Z′↔ν+antiν reaction with quantum statistical flux

    F=bj(1−fi)(1−al),
    R=(1+bj) fi al,
    J=F−R.

For event extent ξ=dt×W×J evaluated *implicitly at the
updated occupations*, set

    bj' = bj − ξ/wV,j,
    fi' = fi + ξ/wν,i,
    al' = al + ξ/wantiν,l.

The solver finds the **unique** root of
ξ−dt W[F(ξ)−R(ξ)]=0 by Brent on the exact physical
Pauli/Bose interval. No after-the-fact clipping, mass-energy
remapping or an independently imposed thermal spectrum.
Inside the interval F decreases and R increases with ξ,
so the function is monotone and a unique bound-preserving
reaction extent exists. Each event preserves reduced charge,
neutrino-antineutrino difference and energy **exactly**, because
E_V,j=kν,i+kantiν,l. The quantum entropy derivative along ξ
is ln(F/R), so this event extent moves in a direction of
nonnegative entropy change. All events share positive W.

**Limit:** a finite lattice with energy-collocated daughters
only approximates continuous angular phase space; angular
endpoint inclusion and Emax truncation generate genuine
discretization bias. Per-event invariants are necessary, not
proof of a physically accurate continuum collision rate.

## Numerical source and nonthermal initial fixture

Reference continuum integrates incoming physical vector
momentum p and CM cosine with independent Gauss–Legendre,
not the discrete event graph. Both use the identical physical
fermion/Bose shell factor for a synthetic (not observational)
nonthermal fixture at T=1:

    fν(k)=expit[−k+0.15 exp(−k)],
    fantiν(k)=expit[−k−0.10 exp(−k/2)],
    bV(p)=0.75/[exp(sqrt(p²+mV²))−1].

Set mV=.4 and Emax=4 in common energy units, Γ0=1.
Independent reference uses 144 incoming-p and 144 angle
Gauss points with pmax=sqrt(Emax²−mV²).

| Energy step h | Positive event edges | dNV relative to continuum | dEV relative | entropy source relative |
|---:|---:|---:|---:|---:|
| .4 | 45 | −21.0550% | −12.6612% | −20.3440% |
| .2 | 189 | −8.0845% | −4.8044% | −8.0491% |
| .1 | 774 | −0.8873% | −0.8330% | −1.3797% |
| .05 | 3114 | **+0.3778%** | **+0.4095%** | **+0.3305%** |

The continuum finite-cutoff reference is
dNV=+0.00580880057694, dEV=+0.01157798731686,
dS=+0.00144991988903 in consistent natural time/volume units.
These numbers are signed *vector* production and
heat; neutrino+antineutrino changes have opposite energy.

**Negative:** exact discrete conservation, positivity and H theorem
do NOT guarantee accurate continuum moments. Even h=.05 is still
off by ~0.3–0.4%, with nonmonotonic bias as h is refined.
The continuum **cutoff** is independently important: moving
Emax from 4 to 8 changes the continuum number source by
~14% and energy by ~29% relative to the Emax=8
reference; Emax=6 remains ~2.48% low in number and
~6.88% low in energy. No global cutoff-converged claim.

## Finite nonthermal evolution and temporal error

With mV=.4, h=.2, Emax=4, 189 events, begin from the
nonthermal fixture above and evolve to Γ0t=4.
The initial ledger is

    C=Nν+Nantiν+2NV = 0.677688967485,
    L=Nν−Nantiν    = 0.00629871876221,
    E_total        = 1.087744100397,
    S_total        = 1.617320525670.

Every tested step preserves C, L and E to absolute
floating-point residual O(10^−16), while all fermion
occupations remain in [0,1], bosons nonnegative, and
the entropy grows after each step. The endpoint at
a common T, μν, μantiν, μV=μν+μantiν is stationary
to floating-point tolerance, also for μν=μantiν=.1.

| dt in Γ0^-1 | steps to t=4 | bV at lowest-energy grid node | ΔS total |
|---:|---:|---:|---:|
| 1 | 4 | 0.8887712453 | +0.0006936347 |
| .5 | 8 | 0.8722947720 | +0.0007487274 |
| .25 | 16 | 0.8615340320 | +0.0007886689 |
| .125 | 32 | 0.8578753698 | +0.0008123882 |

**Negative:** dt=1 vs dt=.125 changes that vector
occupation by ~3.6%. Exact eventwise energy conservation,
Fermi/Bose positivity and monotonic entropy are **not**
a time-accuracy certificate. Serial first-order reaction
splitting requires a dt extrapolation or controlled
higher-order/implicit coupled solve.

An initially EMPTY vector spectrum bV=0 is also
valid: inverse ννbar production populates boson
modes without violating Pauli/charges/energy and
increases total entropy. This pilot can handle
distinct ν/antiν occupations; it does not assume
zero lepton asymmetry.

## Reproducibility and scientific limitations

- [Executable event graph](vector_neutrino_reversible_event_graph.py)
- [Six focused regressions](../tests/test_vector_neutrino_reversible_event_graph.py)
- [Compact JSON](results/vector-neutrino-reversible-graph-2026-11-22.json)
- [Permanent stage archive](../archive/stages/2026/11/2026-11-22-vector-neutrino-reversible-event-graph.md)

    python -m pytest -q tests/test_vector_neutrino_reversible_event_graph.py
    python -m kinetics.vector_neutrino_reversible_event_graph

Numeric units are scaled to T=1, Γ0=1 and h=0.4/0.2/0.1/0.05;
no random seed or external runtime datasets. Python+NumPy,
SciPy; use the Git SHA and exact CLI in the stage archive.
The on-shell vector matrix element and spins are assumed
as in the preceding [full Z′↔ν source]
(VECTOR-NEUTRINO-NONTHERMAL-STRONG-AUDIT.md).

This is a **first numerical event graph**, NOT a validated
nonthermal continuum time solver; uniform physical k/p
grids, angular adaptivity, radial cutoff, Z′→χi pairs,
weak e±γν processes, true time-dependent mediator
mass, H/entropy, χ1/χ2→real φ, BEST and its 86-eV
f(q)/ΩDM/BEC/CLASS all remain open. A coupled
gauge/SM model must replace Γ0=1 and the synthetic
initial distribution before any cosmological rate claim.

**Next step:** investigate continuous-source
finite-volume accuracy and the observed nonmonotonic
h and dt errors; use a conservative higher-order
implicit/IMEX solver and independently verify
strong ν/antiν and vector mode shapes with
multiple angular and cutoff refinements. Only
then introduce the evolving weak SM/χ plasma
and a gauge-valid microscopic production model.
