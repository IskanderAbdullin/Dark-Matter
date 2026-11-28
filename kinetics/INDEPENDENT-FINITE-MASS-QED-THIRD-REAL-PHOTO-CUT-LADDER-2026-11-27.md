# The third archive reproduces the bare real photoabsorption cutoff dependence

Complete the [fine archive photon cutoff ladder](INDEPENDENT-FINITE-MASS-QED-REAL-PHOTO-CUT-LADDER-2026-11-27.md)
under the identical kinematic and thermal conventions for the third
charged and scalar spectra. This remains one real absorption sign,
one charge assignment, a fixed scalar mode and a declared lab
photon energy cutoff.

At the lowest cut 4e−5 the bare signed source is +.05964387
of the positive Born source at g=.005, versus +.05860619
on the fine fixture. The third fixture's maximum two-mesh
relative drift is .001421; its restricted-window overlap
with the older reference is within 2.84e−6 after its charge
factor is removed. The common-temperature net divided by
turnover is −1.32e−16. The audit keeps gain and inverse
loss separate and compares their expected Born soft scales,
not the difference of two independently rounded large rates.

The cutoff ladder is **not** a physical positive correction.
It omits the equal local leading virtual subtraction and the
finite dressed terms. The [separate photon pole and cut gates](INDEPENDENT-FINITE-MASS-QED-FROZEN-TRANSVERSE-PLASMON-GATE-2026-11-27.md)
also explain why this vacuum photon basis cannot be carried
unchanged into the frozen plasma.

Replay in kinetics with the pinned spectra:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_hard_photo_complement_reference.py third
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_hard_photo_cut_audit.py fine third
~~~
[Third receipt](results/finite-mass-qed-hard-photo-angular-cut-third-2026-11-27.json),
[both-archive audit](results/finite-mass-qed-hard-photo-cut-both-audit-2026-11-27.json).
