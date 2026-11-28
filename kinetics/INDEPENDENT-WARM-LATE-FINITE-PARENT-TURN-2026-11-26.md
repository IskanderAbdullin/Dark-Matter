# Finite heavy-parent response on the late .68 GeV turning trajectory

The earlier [turning clock](INDEPENDENT-WARM-LATE-TAIL-TIME-REFINEMENT-2026-11-26.md)
enforced the stationary heavy occupation at every energy in every
collision evaluation. On the same native closed decay/inverse-decay
chart, let the heavy occupation evolve at all 7,730 parent nodes.
For each implicit step the parent equation is linear when the trial
daughter phase is fixed. Eliminate it algebraically:
h_new = (h_old + Δt λ A) / [1 + Δt λ (A+B)],
where A and B retain the full occupation factors and λ includes
the time-dilated width and the angular normalization. The daughter
equations are solved with that updated parent, not its stationary
limit. The stationary archived endpoint is only a Newton guess;
the accepted previous finite-parent phase provides the next target.

We enforce the conservation of each daughter count plus the parent
count and of the total three-species energy. Check entropy with the
Fermi entropy of the finite parent included; the daughters alone
need not gain entropy when forming a parent out of vacuum. The
existing two-species solver keeps its original entropy check when no
additional guard is supplied. Verify the condensed Jacobian by
finite differences and the parent step independently against its
direct number equation with cancellation-aware tolerance. Read all
accepted snapshots back, recompute W independently from phase space,
and record separate SHA-256 fingerprints for daughter and parent states.

## Finite parent with daughter feedback

At the original six-step clock beginning at Hdt = 1e-8, initializing
the parent at its instantaneous balance changes W by at most
2.90e-7 relative to the same discrete stationary evolution. The
positive-to-negative current bracket stays [7.0e-8, 8.5e-8].
Starting instead from an empty parent on this closed chart changes W
by up to 3.05e-6, but leaves the same bracket. The empty parent is
only a sensitivity control. In the section of the 24-step clock
starting at Hdt = 5.125e-8, a separate phase-matched replay gives
the same earlier bracket [6.25e-8, 6.625e-8]. The two corrected
signed currents at that bracket are +6.83213e-13 and −7.98342e-13
in units of initial dark energy per Hubble time. These are discrete
values; the underlying time discretization still moves the bracket.
Replaying all 24 microsteps from the earlier common Hdt = 1e-8
phase independently keeps [6.25e-8, 6.625e-8]. Its two endpoint
currents are +6.83195e-13 and −7.98357e-13, and its largest relative
W displacement is 2.93e-7 over all 24 saved samples.
Counts and energy in the six-step and short finer replays are
preserved to better than 2e-15 relative, including the parent.

## Independent driven angular comparison

As a distinct diagnostic, prescribe piecewise linear entropy
coefficients between the original saved daughter phases and use the
physical angle normalization independently of the native chart.
Integrate the scalar parent lag with linear forcing and linear W
weight, with analytic small-rate limits and a separate Radau check.
There is no feedback in this comparison. For the six-step .68
example the integrated correction is +2.87e-24 of the initial dark
energy for the equilibrium start, and −6.38e-24 for the vacuum
start. Near the finer turning bracket the lag contribution to the
instantaneous signed rate is about +3.1e-17, smaller than either
signed current. For the coarsely sampled .80 GeV daughter path the
driven integrated correction is +7.43e-21, also small relative to
its late W. The driven comparisons are not error bounds on the
fully coupled solve or statements about every initial parent state.

All results hold only for this finite, prescribed background and the
closed decay/inverse-decay channel. They do not cover visible or
other dark sources, unknown returning exterior populations, a
different heavy initialization, or a cosmologically converged tail.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
# Prepare turn6 and turn24 using the linked time-refinement note.
run() { python -W error warm_late_finite_parent_prompt_reference.py EOS.dat "$@"; }
run finite6 turn6/step-{0,1,2,3,4,5,6}.json
run finite6vac turn6/step-{0,1,2,3,4,5,6}.json --initial vacuum
run finite24local turn24/step-{11,12,13,14,15,16,17}.json
run finite24full turn24/step-{0..24}.json
for pair in 'turn6/step-0.json finite6' 'turn6/step-0.json finite6vac' 'turn24/step-11.json finite24local' 'turn24/step-0.json finite24full'; do
  python -W error warm_late_finite_parent_audit_reference.py $pair
done
drive() { python -W error warm_late_tail_driven_heavy_lag_reference.py EOS.dat "$@"; }
drive drive068low.json turn6/step-{0,1,2,3,4,5,6}.json
drive drive068high.json turn6/step-{0,1,2,3,4,5,6}.json --order 8 --pieces 2
drive drive068vac.json turn6/step-{0,1,2,3,4,5,6}.json --order 8 --pieces 2 --initial vacuum
drive drive068fine.json turn24/step-{11,12,13,14,15,16,17}.json --order 8 --pieces 2
# Prepare competing080/tapes using the .80 competing-flux note.
drive drive080low.json competing080/tapes/step-{0,4,5,6}.json
drive drive080high.json competing080/tapes/step-{0,4,5,6}.json --order 8 --pieces 2
~~~

[Finite replays and independent audits](results/warm-late-finite-parent-turn-summary-2026-11-26.json.gz),
[accepted daughter phases from the three shorter replays with heavy-state hashes](results/warm-late-finite-parent-turn-daughters-2026-11-26.json.gz),
[driven parent diagnostics](results/warm-late-tail-driven-heavy-lag-2026-11-26.json.gz).

## Matched finite-parent control on the .80 GeV declining path

On the actual .80 GeV archived phase at Hdt = 1e-7, regenerate all
six original stationary microsteps through Hdt = 1e-6. The regenerated
terminal discrete moments are exactly equal at recorded precision.
The decoded phases seed Newton but do not set the finite parent's
accepted targets. Each of the six finite-parent currents remains
negative. The relative difference in W increases only to 1.01e-8
at the end of this selected fixed-background interval; terminal
currents are about −2.88415e-8 in both runs. This is a feedback control
for the same discrete prompt channel, not the earlier expanding
pair/elastic comparison.

~~~sh
t=competing080/tapes
python -W error warm_late_competing_locked_guides_reference.py EOS.dat "$t/step-5.json" "$t/step-6.json" guides080 --steps 6
run finite080 guides080/step-{0,1,2,3,4,5,6}.json
python -W error warm_late_finite_parent_audit_reference.py guides080/step-0.json finite080
~~~

[Matched .80 audits and regenerated guide phases](results/warm-late-hard080-finite-parent-2026-11-26.json.gz).
