# Cold-window integral of the frozen-background parent relaxation

Integrate the [tagged finite-time relaxation](INDEPENDENT-EARLY-POLE-SCALAR-PARENT-MOMENTUM-SCAN-2026-11-27.md) from q=1e−6 to 1e−5. Weight each tagged fractional gain by the **initial** archived scalar measure q² f_s(q) dq/(2π²). Freeze the daughter and scalar distributions in the parent equation and use the same frozen distributions to evaluate the scalar diagnostic; only the parent can relax within this one channel. The integration panels are split at the known scalar cell edges and at crossings of the hard parent and daughter cell edges. Integrating the weight alone reproduces the archived number in the window within 5e−14 relative.

| Archive | Integrated gain with all occupations frozen | Integrated gain with explicit parent relaxation | Difference |
| --- | ---: | ---: | ---: |
| Fine | 5.57373252e−6 | 5.54530454e−6 | −0.51003% |
| Third | 4.58662524e−6 | 4.56303981e−6 | −0.51422% |

Exposure is 7.2e−5. The comparison at outer momentum order 3→6 and inner radial order 14→24 changes the displayed reduction by less than 0.00030 percentage points. This is the accumulated source evaluated on fixed background occupations, not the evolved abundance in the five-channel calculation. The held daughter and scalar inputs do not receive reciprocal work in this diagnostic, so a closed energy or number evolution is not implied.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
for g in fine third; do
  python -W error early_pole_archive_scalar_parent_window_reference.py "$g" \
    --radial-order 14 --q-order 3
  python -W error early_pole_archive_scalar_parent_window_reference.py "$g" \
    --radial-order 24 --q-order 6
done
```
Raw: `results/early-pole-archive-scalar-parent-window-2026-11-27.json.gz`.
