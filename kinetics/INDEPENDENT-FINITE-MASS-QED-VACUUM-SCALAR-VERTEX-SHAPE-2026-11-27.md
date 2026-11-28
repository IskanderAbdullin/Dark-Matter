# Vacuum baseline: the relative scalar vertex and unresolved real recoil

For equal gauge charges and unequal on-shell fermion masses, compute
the one-photon vertex correction of the neutral scalar current in
vacuum Feynman gauge. Only compare the normalized scalar-current
amplitude at different invariant transfers q², using q²=0 as
the reference (this is not zero four-momentum). At fixed external
masses this difference removes the q²-independent UV term,
external on-shell wavefunction factor, and scalar current
normalization. It does not fix their absolute finite values.

For light photon regulator mass lambda, parametrize the two
fermion propagators by x=u*t, y=u*(1−t). With d=P·b, write

    A = M² t² + m²(1−t)² + 2d t(1−t)
    Delta = u² A + (1−u) lambda².

In four dimensions the projected loop numerator is
`4 L² + 4d + u N1 + 4u² A` for shifted virtual momentum L,
with N1 from the mass and momentum contractions in the code.
An explicit spinor trace checks that coefficient to 8.1e−16
relative. After integrating over L, the UV-independent
momentum comparison uses the numerator divided by Delta
and the logarithmic part from L². The remaining vertex soft
log coefficient is `4d∫dt/A = 4χ coth χ`.
The common wavefunction terms cancel in a transfer difference,
so its logarithmic rate coefficient agrees with the difference
of the [real angular soft factors](INDEPENDENT-FINITE-MASS-QED-SOFT-SUBTRACTION-2026-11-27.md).
This is the relative version of the usual
[vacuum soft real/virtual cancellation](https://inspirehep.net/literature/2313).

Subtract the logarithm at scale M, then combine with the
massive-photon regulated real eikonal integral and with the
finite exact massless real recoil remainder. The latter is
already integrable and does not need the auxiliary photon mass.
The physical real-photon cap is in the parent rest frame;
the regulator is sent to zero at fixed cap. A separate
finite real-band integral reproduces every cap difference
to 2.1e−14 absolute, without using a vertex or regulator.
The regulator ladder also reproduces the analytic subtracted
relative vertex for transfers well beyond the archive poles.

At the two archived scalar poles and rest cap omega/M=.01,
the combined relative shape coefficient (relative to q²=0,
before the common g²/(4 pi²) factor in the rate) is +.003211
and +.003175. The endpoint regulator defects near the archive
poles are below 2.5e−10. These numbers describe transfer
dependence in a specified *vacuum* comparison only. They are
not the magnitude or sign of the total one-loop rate correction,
not an in-medium vertex, and not a finite-temperature
renormalization or hard/soft matching prescription.

Replay from kinetics with the archived spectra:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_vacuum_scalar_vertex_shape_reference.py
```
[Raw receipt](results/finite-mass-qed-vacuum-scalar-vertex-shape-2026-11-27.json).
