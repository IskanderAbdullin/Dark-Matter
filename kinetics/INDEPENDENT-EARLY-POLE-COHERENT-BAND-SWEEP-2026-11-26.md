# Fast instantaneous coherent growth need not survive a narrow-band crossing

As a sensitivity check on the
[homogeneous-field limit](INDEPENDENT-EARLY-POLE-COHERENT-QUARTIC-FLOQUET-2026-11-26.md),
assign a chosen fraction c of the snapshot scalar energy to a
coherent oscillator. Give the remaining sector a static mass
squared contribution (1−c) times the original Hartree shift.
This interpolation is an organizing assumption, not a derivation
of the coherent fraction or a coupled field equation.

The exact Lamé unstable interval narrows at small elliptic
parameter k²; its width in (p/Ω)² scales as (k²)².
Use a rationalized width and the narrow-band profile where
machine precision makes direct trace subtraction unreliable.
At elliptic parameter 0.002 the numerical peak and band
integral agree with the asymptotic expressions within 0.8%.
Integrate the frozen exponent across the band for an illustrative
linear sweep. The chosen sweep is the local value implied if
Ω/T scales as N^(1/3) while the stored net scalar source changes N.
It is not a time-dependent calculation of phase locking.

For the .06 third-grid input at x=0.05, coherent energy fraction
10⁻⁴ gives a local peak rate about 22 Href but only 1.9 × 10⁻⁹
estimated logarithmic amplitude gain in one passage.
For 1% coherent energy the estimate is about 0.20, versus 0.098
on the previous grid. These illustrative one-passage figures
are very sensitive to the assumed sweep. They do not exclude
growth for a trapped or repeated passage, do not show that the
sector actually has coherent domains, and do not predict the
eventual energy transferred by a resonance.

~~~sh
python -W error early_pole_coherent_band_sweep_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Raw sweep controls](results/early-pole-coherent-band-sweep-2026-11-26.json.gz).
