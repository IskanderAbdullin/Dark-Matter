# Finite-mediator contact-EFT validity: exact kinematic bound (2026-11-21)

## Question and scope

The existing full-Pauli χ–electron/neutrino collision benchmarks assume a momentum-independent vector contact coefficient G. This note derives a **strict pointwise tree-level spacelike-propagator bound**, without pretending to compute a physical thermal rate. No mediator mass or coupling has been fixed by experiment; this is a necessary diagnostic, not a sufficient UV validation.

For a massive vector of mass mV, t-channel amplitude has factor 1/(mV²−t). In elastic 2→2 scattering in the CM frame, t=−2p_*²(1−cosθ), hence −4p_*²≤t≤0. Relative to the contact amplitude 1/mV², the exact squared propagator factor is

    R(t) = [mV²/(mV²−t)]² = (1+|t|/mV²)^(−2).

Therefore the contact approximation *overestimates* the differential squared amplitude at every nonforward angle for this isolated t-channel contribution. Its maximum **relative error measured against the contact result** is

    δ_contact,max = 1 − (1+4p_*²/mV²)^(−2).

Requiring δ_contact,max≤ε at all angles gives

    mV/p_* ≥ 2 / sqrt((1−ε)^(−1/2)−1).

For ε=.1 the threshold is mV/p_*≥8.7178 (approximately); for ε=.01 it is ≥28.195 (approximately). These are exact algebraic bounds for the isolated propagator, not cross-section errors in a gauge-completed plasma.

**Alternative relative error convention:** relative to the exact result, the maximum contact overestimate is (1+4p_*²/mV²)²−1. Its ε threshold is mV/p_*≥2/sqrt(sqrt(1+ε)−1), i.e. ≈9.053 for ε=.1. Never mix these conventions.

For elastic masses M,m and invariant s, p_*²=λ(s,M²,m²)/(4s), λ(a,b,c)=a²+b²+c²−2ab−2ac−2bc. A thermal bath has no sharp maximum p_*; any finite contact validity claim must state a phase-space tail tolerance and separately bound the tail contribution. A single benchmark temperature cannot imply uniform validity of G for the full collision operator. The propagator factor does not establish mediator gauge completion, thermal screening, resonance behavior, chemical kinetics, or an 86-eV relic.

## Gauge-charge consistency gate

For a Yukawa interaction y φ \\barχ₂ χ₁+h.c. under U(1)' charges q₁,q₂,qφ, gauge invariance requires qφ−q₂+q₁=0 (for the displayed field ordering). Thus qφ=q₂−q₁. A **diagonal** vector current for two distinct charge eigenstates is not automatically an **off-diagonal** χ₁↔χ₂ coupling. If q₁=q₂, φ may be neutral, but this by itself establishes neither cosmological stability nor Bose production. If q₁≠q₂, a charged φ requires specifying symmetry breaking, spectrum and stability. A vectorlike singlet χ used in the current χe/χν pilot is not yet a UV embedding of the original three-field decay model.

## Falsification / next executable test

Implement R(t) in the existing independent CM 5D full-Pauli χe and χν heat integrals, holding the *same* incoming distributions, spin factors and kinematic cuts. Scan mV/p_thermal and compare against the analytic angular envelope above; independently verify with a direct t quadrature. Track electron and neutrino sectors separately and repeat with split temperatures. Quantify tails outside the contact-safe domain before asserting any physical gamma_el. Include real mediator experimental bounds and the charge assignment above before cosmological extrapolation.

Status: **analytic conditional bound only**; no new 5D numerical verification, no finite-mediator thermal rate, no claim of model viability. CLASS-ready=false.
