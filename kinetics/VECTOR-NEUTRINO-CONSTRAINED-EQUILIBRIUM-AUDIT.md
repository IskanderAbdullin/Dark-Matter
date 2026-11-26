# Independent maximum-entropy equilibrium versus nonthermal Z′+ν+antiν event evolution (2026-11-22)

## Research question and prespecified falsification

An exact quantum 1↔2 event-network audit found two **spurious conserved** highest-energy daughter occupations in the old integer ν grid. The new half-offset ν/antiν energy grid has stoichiometric left-nullspace of dimension exactly three, spanning the real reduced charge C=Nν+Nantiν+2NV, asymmetry L=Nν−Nantiν and total energy E, at the tested finite domains. **Does its actual finite-time nonlinear Pauli/Bose evolution approach the independently derived entropy maximum under precisely those three discrete constraints?**

Predeclared gates: solve maximum-entropy occupations without calling the collision/time code; all three discrete invariants must match initial state to 2e−11 absolute; positive Bose energies and 0<f<1; the stationary entropy-maximizing state must be unchanged by an explicit microscopic implicit event sweep. Long-time arbitrary nonthermal ν and antiν must approach the **same** solution, not a hidden FD thermostat. Preserve residual entropy gap, spectral errors and dt-bias, even if all conservation tests pass. This tests **discrete** kinetics only, not full physical continuous transport, SM weak equilibration or cosmology.

## Independent constrained thermodynamic solution

For the restricted broken-phase Z′↔ν+antiν network with physical degeneracies gV=3, gν=gantiν=2 and exact on-shell E_j=kν,i+kantiν,l, maximize

    S = Σ_j wV,j[(1+bj)ln(1+bj)−bj ln bj]
       −Σ_i wν,i[fi ln fi+(1−fi)ln(1−fi)]
       −Σ_l wa,l[al ln al+(1−al)ln(1−al)]

at fixed physical *discrete* C,L,E. The Bose/FD entropy is strictly concave inside the physical domain. Three Lagrange multipliers yield

    fi = expit[(μν−ki)/T],
    al = expit[(μanti−kl)/T],
    bj = 1/{exp[(Ej−μν−μanti)/T]−1},

subject to μV=μν+μanti<min_j Ej. Because the selected reaction graph has no extra exact left-null invariants in these three finite test cases, there are no additional discrete chemical charges to fit. The thermodynamic solver independently solves three nonlinear constraints for (log T, μν/T, μanti/T) using SciPy root, with NO calls to event source or time integration. The bath temperature/chemical potentials emerge from conserved densities and are **not imposed**.

The same initial synthetic but genuinely unequal arbitrary nonthermal spectra were used:

    fν(k)=expit[−k+0.15 exp(−k)],
    fantiν(k)=expit[−k−0.10 exp(−k/2)],
    bV(E)=0.75/[exp(E)−1].

The relative vector phase-space cell weights use the earlier independently audited **clipped exact wV integral**; daughter shell uses exact positive angular overlap and half-offset k bins. Mass mV=.4, vector Emax=4, energy spacing h=.4,.2,.1. Energies and T in the synthetic unit Tinitial=1 and Γ0=1 sets arbitrary time: NOT a viable mediator coupling, Γ/H or a physical 86-eV model.

## Measured independent equilibrium solutions

| h | Constrained T/Tinitial | μν/T | μanti/T | μV/T | Maximum-entropy gain |
|---:|---:|---:|---:|---:|---:|
| .4 | 1.03058018 | −0.10959360 | −0.16466833 | −0.27426193 | +0.00108105814 |
| .2 | 1.02842113 | −0.10581032 | −0.16095289 | −0.26676321 | +0.00113072711 |
| .1 | 1.02762325 | −0.10440679 | −0.15957578 | −0.26398257 | +0.00115047138 |

All three final C,L,E equal their respective initial discrete values to absolute 4.5e−16 or better. The nonzero **lepton asymmetry** demands μν≠μanti; forcing μν=μanti would violate the real L constraint. At h=.2 initial invariants were C=0.66573105994, L=0.00625489214, E=1.05132922830, initial S=1.57107159701. Their independent equilibrium yields S=1.57220232412; the finite entropy gain is not a numerical arbitrary normalization.

After computing the above root **independently**, passing these equilibrium FD/BE occupations into the on-shell microreversible event sweep leaves every state unchanged to ~1e−13 at dt=3 (tested h=.4,.2,.1), and the microscopic quantum vector weak source vanishes. This establishes a scoped, nontrivial agreement of independent equilibrium thermodynamics with the event graph's detailed balance, not a time-evolved cosmological plasma.

## Actual nonthermal time pilot without neutrino thermal reset

Keep all individual ν, antiν and Z′ occupancies dynamically independent, and advance on-shell positive microreversible events using bounded implicit reaction extents. At h=.2, Γ0t=500 with dt=1, the maximum absolute occupation-space difference from the **precomputed independent constrained FD/BE solution** is ~0.00103721 and the remaining total entropy gap is +1.8630×10⁻⁶. At Γ0t=1000, dt=1, the difference falls to ~0.00033923 and entropy gap +3.9929×10⁻⁷; the simulation continues to conserve discrete C,L,E to floating precision, with nondecreasing entropy throughout.

The h=.4 test also approaches its thermodynamic maximum: at Γ0t=100, dt=1, maximum occupation difference ~0.00261739, entropy production reaches ~99.12% of its theoretical maximum; at Γ0t=1000, dt=1, difference ~4.71×10⁻⁵ and entropy nearly saturates. These are finite-domain, long-time numerical validations at **specific** Γ0=1, not claims of a continuum relaxation eigenvalue or physical γelastic.

## Negative results that cannot be ignored

First-order sequential event-splitting retains a **real temporal bias** even though every elementary event is bounded and conserves energy. At h=.2 and τ=100, lowest-energy bV is 0.75305223 for dt=2, 0.74397501 for dt=1, 0.74667793 for dt=.5 and 0.74874261 for dt=.25. The behavior is **nonmonotonic** in dt and does not justify a Richardson extrapolation based on just two resolutions. Reducing dt is necessary but not sufficient for a certified physically relevant source.

The equilibrium values vary with h: h=.4/.2/.1 yields T/Tinitial=1.03058/1.02842/1.02762 and μV/T=−.27426/−.26676/−.26398. Numerical convergence of *discrete* conserved quantities is not evidence of continuum physical freeze-out. The independent 5D weak/strong and finite-cutoff errors remain from earlier stages: h=.05 number weak bias ~−0.142%, Emax tail and missing threshold bin; the finite lattice's extremely high-energy and angular continuum is not certified. A small entropy gap also does not certify correct physical relaxation rate; a rate-normalization error may change evolution time while leaving the equilibrium invariant.

## Sources and reproduction

- [Independent equilibrium and nonthermal time comparison](vector_neutrino_constrained_equilibrium.py)
- [Focused regressions](../tests/test_vector_neutrino_constrained_equilibrium.py)
- [Compact numeric data](results/vector-neutrino-constrained-equilibrium-2026-11-22.json)
- [Permanent stage record](../archive/stages/2026/11/2026-11-22-vector-neutrino-constrained-equilibrium.md)

Run:

    python -m pytest -q tests/test_vector_neutrino_constrained_equilibrium.py
    python -m kinetics.vector_neutrino_constrained_equilibrium --h 0.2 --tau 500 --dt 1

No random seed. Python NumPy/SciPy/pytest; deterministic event ordering; no external state or temp Action artifact needed to reproduce. Context: [neutrino quantum kinetic collisional evolution](https://arxiv.org/abs/2008.01074), [prior exact stoichiometric audit](VECTOR-NEUTRINO-EVENT-INVARIANTS-AUDIT.md). No external implementation copied.

**Next P0:** verify real Scientific CI for constrained equilibrium and previous topological audit; then investigate opposite reaction order/Strang or coupled implicit integrator to quantify and reduce temporal bias, matched 3+ h/dt/cutoff full spectra and finite ν weak chemistry. Only after nonthermal strong source and gauge valid χ1/χ2/φ+SM+Z′ plasma integration may one assess physical Γel/H, 86-eV Bose freeze-out, ΩDM, BEC or CLASS/P(k). This stage is **not terminal**.
