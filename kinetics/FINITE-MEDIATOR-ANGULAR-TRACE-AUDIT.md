# Nonconstant vector spin-trace angular reference (2026-11-22)

## Predeclared scope and criterion
For tree-level elastic Dirac χ+lepton scattering at fixed center-of-mass invariant s, with vector-current numerator
`N(t)=D²+U²+2(M²+m²)t`, `D=s−M²−m²`, `U=−D−t`,
test an independent analytic integration of the **nonconstant numerator**
times the massive spacelike propagator factor against direct adaptive numerical
quadrature at massive and massless lepton points. Relative tolerance 3×10⁻¹¹
for seven specified kinematic points; test independently in public CI.

With `y=1−cosθ*`, `a=2p*²`, `b=a/mV²`, `t=−ay`, and
`B=M²+m²`, the numerator is
`N(y)=2D²−2a(D+B)y+a²y²`.
Therefore the exact angular moment (without the overall `2G²`) is
`A=2D²J0−2a(D+B)J1+a²J2`, where
`Jn=∫₀² yⁿ/(1+by)² dy` and
`J0=2/(1+2b)`,
`J1=[ln(1+2b)+(1+2b)⁻¹−1]/b²`,
`J2=[2b−2ln(1+2b)−(1+2b)⁻¹+1]/b³`.
At small b, direct subtraction is unstable, so the implementation uses
`Jn=Σₖ(−1)ᵏ(k+1)bᵏ2^(n+k+1)/(n+k+1)`, convergent for `2b<1`.

## Independent verification and limitations
The new source is `finite_mediator_trace_angular.py` and the independent
adaptive `scipy.integrate.quad` tests are
`tests/test_finite_mediator_trace_angular.py`. In a separate local
calculation, the three moment expressions at b=0.1,1,10 matched numerical
quadrature to relative error no greater than 3.5×10⁻¹⁴. The committed
seven-case physical trace test is **pending actual GitHub CI**, not yet
measured. The previous finite-mediator 5D commit `64d6321` passed
[Scientific verification #44968592762](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/44968592762).

**Crucial non-factorization:** this A(s) is not the angular integral of
finite-temperature heat. Boosted outgoing χ and lepton energies depend on
scattering polar **and azimuthal** angles; final Pauli factors and
`E′−E` prevent replacing the full 5D integrand with A(s). This
diagnostic checks the amplitude/propagator part only. It does not establish
χe/χν heat convergence, a physical relaxation rate, gauge-consistent
thermal mediator, 86-eV relic, BEC, ΩDM, or CLASS readiness.

## Executable next step
Run `python -m pytest -q tests/test_finite_mediator_trace_angular.py tests/test_finite_mediator_5d.py`
in public Scientific verification and record the real run ID. Next compare
full 5D finite-mV heat against an independently coded direct-t integration
**retaining boosted final-state Pauli and energy transfer**; scan
n=12/16/24/32, cut=18/26/36 T, mV/T=1/3/10/100,
massive e and chiral ν, and split Tν. Do not replace full thermal
kernel by the analytic vacuum angular moment.
