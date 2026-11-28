# Direct FD histories through both local-clock temperature labels

Continue the exact .55 restart to .91, then to log(a/a_start)=.96, using direct vacuum charged-lepton timelike and elastic terms. The .91 temperatures (about 0.38624 GeV) remain above the earlier fast-clock labels, so .91 alone does not establish coverage. Both reach about 0.36932 GeV at .96. The first attempted BDF tail at the earlier tolerance exceeded the entropy closure gate at .96; rerunning the identical tail with tighter tolerances and bounded steps preserves the original gates. The final intensive BDF/Radau paths agree within 6.26e-11 over that tail. An independent extensive-number-and-energy Radau integration resolves the intensive state coordinates within 5.02e-13. All these checks retain the existing pressure-work and entropy ledgers.

Locate each earlier linearized rate/H equality merely as a temperature label on the actual, continued trajectory. For the combined case the slow label is at log(a/a_start)≈0.519466, where live eta is -0.034813 rather than the instantaneous algebraic -0.19641. The fast label is at log(a/a_start)≈0.921075, where live eta is -0.075516 rather than -0.90402. Both frozen values are outside a small-force approximation. These are not times inferred from an equilibrium condition changing the evolution.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results; e=eos-range.dat
before="$p/fd-portal-direct-expansion-tracking-continuation-$d.json"
clocks="$p/fd-portal-local-mode-clock-crossings-$d.json"
first="$p/fd-portal-direct-expansion-055-091-$d.json"
final="$p/fd-portal-direct-expansion-055-096-$d.json"
ext="$p/fd-portal-fast-clock-tail-conserved-coordinates-$d.json"
python -W error fd_portal_direct_expansion_beyond_local_clocks.py "$e" "$before" "$clocks"
python -W error fd_portal_direct_expansion_fast_clock_tail.py "$e" "$first" "$clocks"
python -W error fd_portal_fast_clock_tail_conserved_coordinates.py "$e" "$first"
python -W error fd_portal_fast_clock_tail_coordinate_closure.py "$final" "$ext"
python -W error fd_portal_live_states_at_local_clock_temperatures.py "$e" "$before" "$final" "$clocks"
```

Use the earlier local EOS normalization. This is a controlled direct numerical continuation, not a complete cosmological implementation or a freeze-out threshold or abundance prediction.
