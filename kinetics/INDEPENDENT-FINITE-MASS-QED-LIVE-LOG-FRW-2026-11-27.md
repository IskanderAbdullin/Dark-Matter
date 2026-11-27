# A live conservative update with the finite-mass log energy

The [hard-shell log reference](INDEPENDENT-FINITE-MASS-QED-SHELL-LOG-2026-11-27.md)
tested fixed occupations. Here the same symmetric energy changes
with live two-moment Fermi reconstructions and is tested with
[the original comoving decay/pair packet fixture](INDEPENDENT-EARLY-POLE-LIVE-JOINT-HARD-PACKETS-FRW-2026-11-26.md).
The photon term, Hartree scalar energy, tree packet weights, and
imposed expansion/interaction clock are unchanged. This is a
finite active-domain numerical control; it adds no exterior term.

Integrate the quadratic log on unordered cell pairs, using one
triangle for the same cell and continuous corner transformations
for adjacent cells. Differentiate the *same* quadrature with
respect to the two entropy fit coefficients on both endpoints of
each pair. Transform that gradient back to the evolved moments
with the exact susceptibility of the Fermi moment map. The
interaction contribution to the comoving stress-energy trace is
the same integral weighted by
`2−(a m/E_p)^2−(a m/E_k)^2` at fixed cell moments. The
derivatives go to the same packet energy remap and the FRW
scale-work integral, not to a separate energy adjustment.

At two scale factors, direct log-only finite differences agree
with the explicit energy and trace derivatives within 3.4e−25
absolute. Changing the log quadrature from order 14 to 24 changes
the log moment gradient and trace by less than 2.1e−9 relative
on the tested state. Across both packet quadratures and the
static-scale control the complete redshift-subtracted energy
drift stays below 1.7e−16, with positive entropy gain. In the
prescribed expanding runs, adding the real log reduces the
final scalar gain by 3.48e−7 and 3.50e−7 relative; the two order
variants have different base gains, so this is a paired response,
not a continuum error bound.

This does not match a finite-mass collision operator, quasiparticle
residues, or dielectric widths. It does not establish a realistic
cosmological rate or a Friedmann solution. The earlier dressed
energy, packet truncation, and tagged-domain assumptions still
apply; in particular the [archived exterior term](INDEPENDENT-FINITE-MASS-QED-LOG-EXTERIOR-2026-11-27.md)
has to be included before using the same calculation as an archived
live/exterior comparison.

Replay in kinetics at both packet quadratures:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_live_log_frw_reference.py --level 0
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_live_log_frw_reference.py --level 1
```
Receipts: [level 0](results/finite-mass-qed-live-log-frw-l0-2026-11-27.json),
[level 1](results/finite-mass-qed-live-log-frw-l1-2026-11-27.json).
