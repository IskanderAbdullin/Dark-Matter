# Full 5D χ–lepton finite-vector angular convergence pilot (2026-11-22)

## Why this calculation matters

The independent adaptive polar QUAD control at a *single incoming
event* already refuted 48-point angular convergence for a very
light vector mediator. The next physics gate is to check whether
the error persists after integrating the full thermal incoming
phase space, rather than treating one unusual momentum configuration
as representative.

This new stand-alone 5D implementation evaluates physical Pauli
weak χ energy transfer with **no fixed BGK rate**. The incoming
p,k and cosine μ use Gauss quadrature; at *every* incoming phase-space
node the outgoing CM polar angle is integrated independently by
adaptive SciPy QUAD, by Gauss cosine, or by the rational transfer-u
map, while outgoing azimuths reconstruct boosted χ/e energies.
All channels use full vector Dirac spin trace, t-dependent
[mV²/(mV²−t)]² factor, and full FD forward/reverse. A full-domain
1/2 avoids double-counting forward and reverse reactions.

The integrand at fixed incoming p,k,μ and relative azimuth is
shared across polar methods, so this is a **quadrature**
cross-check, *not* an independent matrix element nor a proof of
nonthermal strong collision dynamics. One regression also compares
the conventional 5D cosine normalization with the pre-existing
independent_massive_chie_5d.py implementation.

## Conditions and numerical protocol

All masses and energies are in consistent GeV benchmark units:
Mχ=.1, mℓ=0 (two-flavor active neutrino massless proxy),
Te=1, Tχ=.9, G=.001 GeV^-2,
incoming kinetic cutoff=14 Te, gχ=gℓ=4.
Mediator mV=.03,.1,.3. Polar adaptive epsabs=1e−25,
epsrel=2e−10, limit300, 16 evenly spaced CM azimuths.
Direct/transfer Gauss polar order32. Incoming orders
(nrad,nμ)=(8,8) and (10,10) were both measured.

## Measured numeric heat (GeV5; same radial grid per comparison)

| mV | incoming n | adaptive polar Qχ | cosine32/adaptive−1 | transfer32/adaptive−1 |
|---:|---:|---:|---:|---:|
| .03 | 8 | +1.0032391470910e−15 | −7.831% | −9.531% |
| .03 | 10 | +1.0032302358139e−15 | −7.876% | −9.559% |
| .10 | 8 | +8.1920040675003e−14 | +0.8032% | −0.9619% |
| .10 | 10 | +8.1926942054035e−14 | +0.8124% | −0.9102% |
| .30 | 8 | +3.7414271378324e−12 | +0.1783% | −0.04141% |
| .30 | 10 | +3.7444785932879e−12 | +0.1775% | −0.04482% |

**Quantitative negative:** even after the full incoming
(p,k,μ) integration, the *transformed* polar Gauss32 heat
can be about **9.6% too low** for the mV=.03 mediator,
relative to adaptive polar at the *same incoming grid*.
The cosine32 estimate is ~7.9% too low; the rational
transfer transform is therefore **not uniformly better**
than direct cosine for signed Pauli energy-transfer moments.
At mV=.1, the two 32-point methods err in opposite
directions by ~1%. A correct t Jacobian and equal-T
detailed balance never imply convergence of Qχ.

Radial order 8→10 changes the adaptive reference by
about −8.9e−6 relative (mV=.03), +8.4e−5 (mV=.1)
and +8.2e−4 (mV=.3). These are **not yet
fully converged** radial/cutoff studies; they
indicate that the polar discrepancy is appreciably larger
than this specific radial refinement for the light mediator.
The test isolates the angular issue, not the final
cosmological uncertainty.

Independent equal-temperature check (mV=.1, incoming
nr=nμ=5, massless neutrino or massive mℓ=.5):
Qχ(Tχ=.9)>0, Qχ(Tχ=1.1)<0, and Qχ(Tχ=Te)
consistent with zero to floating roundoff. No Fermi
equilibrium heat should be inferred from a raw positive
forward collision rate. This is a valuable thermal weak
moment reference, not the full nonthermal kinetic spectrum.

## Reproduction and next gates

Source: kinetics/finite_mediator_5d_adaptive_reference.py
Regression: tests/test_finite_mediator_5d_adaptive_reference.py
Compact data: kinetics/results/finite-mediator-5d-adaptive-2026-11-22.json

    python -m pytest -q tests/test_finite_mediator_5d_adaptive_reference.py
    python -m kinetics.finite_mediator_5d_adaptive_reference --mV 0.03 --nr 8 --nmu 8 --nang 32

**Next scientific P0:** extend incoming nrad/nμ, azimuth,
polar adaptation and both incoming/outgoing tails at
finite T/m_e for physical electrons and chiral ν, including
split e/ν temperatures. Construct an absolute error
budget for near-canceling e and ν heat, not a ratio of
two unconverged positive numbers. Then constrain the
gauge portal mediator, χ1χ2φ charge consistency,
e±γν SM evolution and the full nonthermal Pauli/Bose
collision + FRW. No physical Γ_el, macroscopic
BEC, 86-eV ΩDM, frozen f(q) or CLASS readiness.
