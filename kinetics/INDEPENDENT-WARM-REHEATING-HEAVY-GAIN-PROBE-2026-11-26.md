# Test whether heavy regeneration is perturbatively small on the prompt lepton spectra

The [weighted-count upper](INDEPENDENT-WARM-DIRECT-RESTRICTED-COUNT-UPPER-2026-11-26.md)
explicitly omits dark creation of new heavy pairs. A
separate calculation asks whether that omission is
numerically justified as a small correction on the
**imposed direct-lepton** spectra at 50 and 100 MeV
reheating. It is not a bound on the full signed rate.

Produce the stable fermion via direct vacuum FD
lepton pairs and prompt heavy decay; produce the
scalar from the same decays. Redshift the particles
between production and each labeled bath temperature.
Reconstruct positive momentum counts on two independent
histogram sizes, then fold cross sections over available
lab incident angles. Use the spin-antisymmetrized
same-charge channel, the opposite-charge channel,
and the scalar-pair channel separately. Ignore outgoing
heavy Pauli blocking to isolate the gross kinematic
opportunity for regeneration. Each regenerated pair
adds two units of the weighted count, because the
two incoming light particles each carry weight one
and each heavy carries weight two. Prompt subsequent
decay redistributes, but does not further increase,
the ledger under that convention.

| Reheat → measured T (MeV) | Gross extra Q per Hubble time, divided by visible entropy | Stable input processed per Hubble time, relative to passive stable count | Scalar input processed, relative to passive scalar count |
|---|---:|---:|---:|
| 50 → 35 | 1.60e−5 | 1.27 | .199 |
| 100 → 75 | .00393 | 33.20 | 23.50 |
| 100 → 60 | .00512 | 34.50 | 16.80 |
| 100 → 50 | .00439 | 28.11 | 9.61 |

The processing numbers are rates per Hubble time on
fixed distributions, **not** claims that a fixed
population can actually be used repeatedly. For
example, a regenerated stable pair may return two
stable daughters, but their momenta are different.
With isotropic unpolarized prompt decay one scalar
takes the mean energy share E*_scalar/M of its parent.
At 100 → 60 MeV, the two stable-incoming channels
then transfer more energy to scalars than the
scalar-incoming channel removes. This net rate is
about 45 times the passive scalar energy per Hubble
time. It is balanced by the opposite sign in stable
energy. Other reactions or inverse rates can change
the result once the distributions respond.

The most active labeled epochs change the gross
count by less than 0.2% between the sampled grids;
the already tiny late tail is less well converged.
Recover both stable charge rates against an independent
existing fold to 5e−13 relative on the same atoms,
and check the scalar rate with a separate lab-angle
projection to better than 0.1%. The resulting rates
are too fast to justify evolving these fixed spectra
as a first-order approximation at 100 MeV. The test
does not estimate a final yield, the sign of all
coupled temperature changes, or whether regeneration
closes the missing abundance. It identifies the need
for population and momentum feedback before that claim.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_heavy_gain_probe_reference.py EOS.late.dat local-spectra
~~~

[All rates, energy ledgers and numerical controls](results/warm-prompt-reheating-heavy-gain-probe-2026-11-26.json).
