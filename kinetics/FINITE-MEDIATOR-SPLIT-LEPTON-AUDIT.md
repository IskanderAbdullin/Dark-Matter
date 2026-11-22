# Split electron–neutrino finite-mediator 5D heat: error budget (2026-11-22)

## Conditional physical question

A previously examined anomaly-free U(1)_(L_e−L_μ) **candidate**
requires vector χ scattering off both e± and active ν_e/ν_μ.
If photons/electrons are hotter than χ but neutrinos are colder,
their heat-transfer moments can have opposite signs. A single
electron BGK or undifferentiated SM reservoir can consequently
produce the wrong sign. Here the finite-vector t-dependent
tree amplitude is included in an independently adaptive
five-dimensional full-Pauli weak heat integration.

**Crucial limits:** the mediator couplings/mass and original
χ1↔χ2+φ gauge charges remain unvalidated. Two active ν flavors
are counted with g=4 total, and e−+e+ with g=4. We assign
the illustrative instantaneous post-e± reheating ratio
Tν/Tγ=(4/11)^(1/3)=0.7137658555; this is NOT a computed
neutrino-decoupling history, nor an assertion that χ1 actually
couples with the assumed common G. The benchmark Mχ=.1m_e
(~51 keV), not the 86-eV φ. No nonthermal strong C[f] or
Γ_el eigenvalue is inferred.

## Numerical setup

Natural benchmark mass unit me=1, Mχ=.1, fixed zero-transfer
G=.001 m_e^−2, Tγ=Te=1, Tχ=.9, Tν=.7137658555.
Each χ+lepton contact channel uses the full massive-vector
propagator [mV²/(mV²−t)]² with t≤0,
all χ/lepton FD forward–reverse Pauli factors, exact
COM-to-thermal-frame outgoing boosts, and finite
physical energy transfer. Incoming p,k,μ Gauss
nrad=nmu=8, outgoing 16 azimuths, kinetic cutoff 14
times the corresponding lepton T. The polar integral
is independently adaptive SciPy QUAD or 32-node transfer-u
Legendre at otherwise identical incoming grids.

## Measured physical weak moments (m_e^5 in chosen units)

| mV/me | Qχ←e, adaptive | Qχ←two ν, adaptive | Total, adaptive | Total, transfer32/adaptive − 1 |
|---:|---:|---:|---:|---:|
| .03 | +6.98337798520e−16 | −9.04684963424e−16 | −2.06347164905e−16 | +2.9624% |
| .1 | +5.81526778486e−14 | −7.21469485358e−14 | −1.39942706872e−14 | +1.2707% |
| .3 | +2.73727704625e−12 | −3.15980864410e−12 | −4.22531597843e−13 | +0.07743% |

For mV=.03, separate n32 transfer angle errors are:
electron −10.7065%, ν −7.5888%; the *net* heat ratio
happens to shift by only +2.9624% due to partial
cancellation. **Do not conclude** the individual kernels
are 3% accurate. At mV=.1 the separate biases are
−1.0021% e and −0.5612% ν; the net bias is +1.2707%.
Even an apparently stable signed total may mask much
larger unresolved channel errors.

On an independently refined incoming nrad=nmu=10 grid
for mV=.03 (same adaptive polar/nphi16/cutoff), the
electron heat is +6.99162148219e−16,
neutrino −9.04441370849e−16,
total −2.05279222630e−16. The total differs
by roughly **0.52%** from incoming order8.
This is an explicit *negative convergence warning*:
even the adaptive angular **total** is not yet a
precision cosmological source, because radial and
thermal-tail errors remain.

With equal temperatures Tν=Te=1, both channels
would instead heat the cooler χ and total Qχ>0
in the same conditional model. Thus the
assumed split ν/e thermal state can reverse the
sign. This is a counterexample to **single
SM bath as a universal closure**, not a
physical prediction of the sign in the real
early universe.

## Reproduction and next research stage

Source: kinetics/finite_mediator_split_lepton_reference.py
Tests: tests/test_finite_mediator_split_lepton_reference.py
Compact data: kinetics/results/finite-mediator-split-lepton-2026-11-22.json

    python -m pytest -q tests/test_finite_mediator_split_lepton_reference.py
    python -m kinetics.finite_mediator_split_lepton_reference --mV 0.1 --nr 8 --nmu 8

Next priority: independently converge nrad/nmu/nphi/polar
and both incoming/outgoing tails at finite electron mass,
massless left-helicity neutrinos, and genuinely evolved
Tν(T),Tγ(T), electron chemical potentials and
e±↔γγ + weak ν collisions, under model-consistent
mediator constraints. Propagate **absolute** errors in
Qe and Qν through the net when near cancellation.
Then test gauge-consistent χ1,χ2,φ production
and full nonthermal Bose/Pauli FRW.
No physical Γ_el, 86-eV ΩDM, BEC, frozen f(q)
or CLASS-ready result has been established.
