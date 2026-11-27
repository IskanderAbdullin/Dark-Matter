# Bound the rate from existing boosted parents and test hard reactivation

At the [7.394 MeV cold recoil gate](INDEPENDENT-WARM-LATE-COLD-SCALAR-PROMPT-GATE-2026-11-26.md)
a parent can emit a scalar below p/T=.1 only if its energy
exceeds 3.24662 GeV. The [one-use branching bound](INDEPENDENT-WARM-LATE-COLD-SCALAR-SINGLE-USE-2026-11-26.md)
is not an instantaneous rate bound when the heavy population
already exists. In an open unpolarized decay at boost gamma,
multiply its vacuum cold angular fraction by the lab rate
Gamma/gamma. The Bose factor in the cold band is held at its
archived maximum. Pauli blocking only reduces this partial
rate, so no lower bound on the total blocked width is needed.

The minimum daughter energy satisfies
E_scalar,min ≥ E_scalar,rest/(2 gamma), and the full interval
width is 2 p* sqrt(gamma²−1) ≥ 2p* beta_min gamma on the
allowed parent domain. For rate per parent (r=1), or rate
per unit parent energy/M (r=2), bound the vacuum coefficient
by (C−E_scalar,rest/(2 gamma))_+ divided by
2 p* beta_min gamma^(r+1). Its unrestricted positive maximum is
at gamma=(r+2) E_scalar,rest/[2(r+1) C]. A separate refined
rapidity grid checks the bound against the exact lab endpoints;
at p/T=.1 it lies less than 0.061% above the sampled vacuum
maxima. After the archived Bose factor, the rate coefficients
are at most 1.0841e−6 times the rest width per eligible parent,
or 2.3777e−8 times the rest width per unit of eligible energy/M.
They are coefficients, not constraints on the abundance or
time history of high-energy parents.

## Reactivate moment-preserving hard packets

At the corrected q64 archive, reintroduce only the recorded
exterior cohorts with their same individual count and mean birth
energy, using admissible two-moment shapes. All tests use a
7.5 GeV lab cap, half-width .02 GeV, and unchanged cold cells.
The reconstructed exterior energy is about 2.826e−11 of the
dark energy in each branch. These deliberately hard packets
select possible stored moments, not measured exterior shapes.

| Packet center (GeV) | Gross prompt gain below .1T | Gross loss below .1T | Net |
| --- | ---: | ---: | ---: |
| .713 | 0 | 0 | 0 |
| 1.55 | 0 | 0 | 0 |
| 3.3 | +5.391e−14 | 2.177e−11 | −2.172e−11 |
| 6.5 | +9.923e−14 | 8.612e−11 | −8.602e−11 |

Rates are divided by represented scalar count times H.
Closed and open comparisons agree at this precision for
these phases. The hard cases open the channel but lose
rather than gain cold scalar number. An independent high
parent fold measures the eligible stationary heavy energy
as 1.632e−19 and 8.409e−20 of the dark energy, respectively.
The rate inequalities bound the corresponding positive
gains by 3.267e−13 and 1.683e−13 in the table units, enclosing
the direct folded gains. The upper is applied to this
gross gain, not the signed result.

No finite duration is assumed for either frozen current.
Rapid decay and finite numbers make a momentary rate
times one Hubble time unsuitable as a population yield.
Decay of an untracked heavy reservoir, repeated use with
changing shapes, other channels and total later production
remain separate questions.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/warm-corrected-q32-time48-step48-2026-11-26.json.gz
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
python -W error warm_late_cold_parent_time_dilation_bound_reference.py EOS.dat "$a" "$b"
python -W error warm_late_cold_scalar_reactivated_prompt_reference.py EOS.dat "$b"
~~~

[Rate limits](results/warm-late-cold-parent-time-dilation-bound-2026-11-26.json);
[reactivated currents and heavy-energy checks](results/warm-late-cold-scalar-reactivated-prompt-2026-11-26.json).
