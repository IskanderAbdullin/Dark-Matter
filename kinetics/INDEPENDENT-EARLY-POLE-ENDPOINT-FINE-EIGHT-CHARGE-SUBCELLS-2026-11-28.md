# Eight charge subcells reduce the next number step but leave the scalar shape step comparable

Repeat the full fine charge-only trajectory with eight charged subcells per original cell and the same 67-bin scalar grid. The archived startup and all fixed clock joins use the same step counts as the previous two- and four-subcell histories, with the full paired number, Pauli, and continuity checks. This is an independently evolved history, not a projection of the final d4 state.

| x | Fine d8 scalar yield | d8 relative to d4 | Last fixed-start time comparison |
| --- | ---: | ---: | ---: |
| .01 | 3.47415512e−7 | −.00459% | 1.30e−7 |
| .02 | 1.01901173e−5 | −.08442% | 1.69e−6 |
| .04 | 1.68971140e−4 | −.82412% | 1.28e−4 |
| .05 | 3.55909633e−4 | −1.37587% | 3.26e−5 |

At x=.05 the latest absolute yield decrement is .4757 times the previous decrement; that previous decrement was 1.464 times the first. These ratios are not evidence of a common asymptotic order. The remaining charge-only yield step is still larger than the measured last fixed-start time difference; the latter is not a global time bound. On the common **unchanged** scalar grid, normalized full histogram TV is .8931% from d2 to d4 and .8650% from d4 to d8. Cold conditional W₁ is respectively 1.3834e−8 and 1.3031e−8. Shape contrasts are therefore still of comparable size even as the yield decrement shrinks. At a common fixed mass for moment postprocessing, d4→d8 absolute total count falls by 1.376%, while the pressure numerator rises by .0432%; their different warm/cold weights matter. None of these tests certifies a physical spectrum.

The d8 pole mass is .0008606946, down 1.303% from the original fine archived value. At the upper cold edge and cosine zero the kinematic interval is E/T=21.9742–21.9800. The largest threshold over the whole cold q/angle band is 29.947, still below the retained parent ceiling near 32. This threshold measures access to the archive's support, not occupancy or a rate.

Reproduce a fine `--div 8 --steps 2048 --one` startup with the endpoint cell driver to x=.0075, then exact continuations with coarse 512, 1024, 1024, 1024 at the four later epochs. See the [continuation](results/early-pole-cell-endpoint-fine-d8-to-x05-2026-11-28.json), [full four-level ladder](results/early-pole-cell-endpoint-fine-charge-ladder-d8-2026-11-28.json), [matched scalar shape](results/early-pole-cell-endpoint-fine-d8-shape-contrasts-2026-11-28.json), and [threshold audit](results/early-pole-cell-endpoint-fine-d8-threshold-2026-11-28.json). The shape script takes d2, d4, and d8 saved x=.05 states on the same fine chart. The independent third d8 history is a separate test.
