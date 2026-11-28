# A coherent limiting branch has a distinct self-resonance

Instead of interpreting the snapshots as quasiparticles, make the
opposite **conditional** assumption: put the same instantaneous
scalar energy into a homogeneous real field in the potential
m₀²φ²/2+λφ⁴/24. At its turning point the energy determines
amplitude A. Its conservative local solution is A cn(Ωt|k²),
with Ω²=m₀²+λA²/6 and k²=λA²/(12Ω²).
Linear perturbations see the periodic curvature of the potential.
The one-period transfer matrix has unit determinant. Its
eigenvalue above unit magnitude defines the Floquet exponent.

The pure quartic band reproduces the Lamé result
3/2 < (p/Ω)² < √3 and a peak amplitude exponent 0.03598 Ω.
See [Greene, Kofman, Linde and Starobinsky](https://arxiv.org/abs/hep-ph/9705347).
For the .06 snapshot energy at x=0.05 on the third grid, the
conditional unstable comoving momentum interval q=p/T
is 0.001066–0.001146 and the peak linear amplitude rate is
about 1.30 × 10¹⁰ in units of radiation-reference H for the
illustrative 86 GeV parent mass. At x=0.10 the band maps to
q≈0.00193–0.00207 with the same dimensionless exponent.
The exact small bare mass is included, and the monodromy and
solver tolerances are checked.

This calculation starts with all the scalar energy coherent, not
with a measured coherent fraction from the spectra. It does not
determine the actual phase correlations, onset of coherence,
backreaction, nonlinear saturation or fragmentation. The
conditional instability is not the same quantity as the bare
particle 4→2 rate, cannot be added to it, and does not predict
the final abundance. Very small q alone does not justify
discarding field evolution if coherence actually develops.

~~~sh
python -W error early_pole_coherent_quartic_floquet_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Known-band and snapshot-energy controls](results/early-pole-coherent-quartic-floquet-2026-11-26.json).
