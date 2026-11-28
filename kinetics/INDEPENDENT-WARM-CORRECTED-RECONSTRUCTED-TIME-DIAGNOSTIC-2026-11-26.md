# Recover admissible moments after first-order time cancellation

Use the completed [matched q32 time ladder](INDEPENDENT-WARM-CORRECTED-TIME-AND-SUPPORT-2026-11-26.md).
At each of its 13 common checkpoints form two first-order diagnostic
combinations, `2*U24-U12` and `2*U48-U24`. Here `U` comprises
comoving number and **kinetic** energy in each shell. Do not treat
an algebraic linear combination of pressures or distribution functions
as an evolved state. Instead reconstruct a regular Bose or Pauli-admissible
Fermi phase in each shell, with no clipping of inputs to the closure.
Integrate its pressure and recheck number and kinetic energy on a
separate 192-point radial rule.

All extrapolated states pass moment and per-mode admissibility checks.
The worst independent relative shell moment defect is 6.24e-14.
Aggregate new exterior particle counts extrapolated separately
remain nonnegative and close the combined count ledger within
5.54e-13, without constructing an exterior phase distribution.

At the common 7.394049 MeV endpoint, the scalar/light number
differences between the two extrapolants are 9.93e-6 / 7.27e-6
in shell L1 over their initial moments. The corresponding
kinetic-energy differences are 9.57e-6 / 7.24e-6, and the
reconstructed pressure differences are 9.57e-6 / 6.98e-6.
These endpoint differences are about 50–60 times smaller than
the raw 24/48 differences. Across all 12 noninitial checkpoints
the two-moment improvement is at least a factor of 6.99.
The two extrapolated cold occupations (4.09167394e-5 and
4.09153893e-5 at the endpoint) continue to decrease at every
common frame. This is a consistency check on the apparent
first-order error, not an accepted 96-step trajectory or proof
of accuracy at a claimed tolerance.

With the exact `ladder/` output from the linked time comparison:

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  warm_multishell_q32_extrapolated_moments_reference.py EOS.dat \
  ladder/seed32.json ladder/q32-12 ladder/q32-24 ladder/q32-48
~~~

[Thirteen checked frames](results/warm-corrected-q32-extrapolated-moments-2026-11-26.json.gz).
