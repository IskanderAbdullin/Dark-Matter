# Finite-mediator angular control — 2026-11-21

This is a deliberately restricted, **independently quadrature-checked propagator-only** benchmark, not a physical thermal collision integral.

For CM momentum p*, u=(1-cosθ)/2 in [0,1] and x=4p*²/mV², the isolated t-channel vector propagator squared divided by its contact limit is R(u)=(1+xu)^(-2). If and only if the remaining angular numerator is constant, its isotropic mean is

I0(x)=∫₀¹du R(u)=1/(1+x).

The normalized transport-weighted mean (weight 1-cosθ=2u) is

I1(x)=2∫₀¹du u R(u)=2[ln(1+x)-x/(1+x)]/x²,

with I0(0)=I1(0)=1. The code evaluates I1 with a convergent small-x expansion to avoid cancellation. These ratios are **not** valid as total χe or χν cross-section ratios when spin traces, Pauli occupations, CM boosts and energy-transfer weights depend on scattering angle.

Independent composite Simpson integration (32768 intervals, standard-library only) over mV/p*=0.1,0.3,1,3,8.59926345,10,30,100 gives maximum absolute analytic/quadrature difference below 2×10^-9 (actual maximum ≈7.4×10^-12). At the pointwise 10%-contact-envelope threshold mV/p*=8.59926345, the angle-integrated ratios are less suppressed than the maximal pointwise defect. This is a **negative control** against interpreting the 10% pointwise criterion as a measured 10% thermal-rate error.

Reproduce: `python kinetics/finite_mediator_angular_control.py`. The script asserts all listed cases and prints JSON to stdout, with no committed large artifacts.

Next required physical step: replace the contact coefficient in the existing 5D χe/χν full-Pauli heat integrands by 1/(mV²-t), retain exact angular spin traces and split e/ν temperatures, perform independent t/angular and thermal-tail convergence, and check the actual χ1χ2φ charge/mediator model. No Γel, freezeout, ΩDM, BEC or CLASS-ready result follows from this control.
