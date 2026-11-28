# Where the retained parent-energy cap clips the cold angular domain

The accepted conditional cold-band integrals hold (k=.0015). Change **only** k in a kinematic pilot, refitting the transverse pole at each k to the same archived charge state. Keep the fitted scalar mass fixed for each state. On the tested past-spacelike branch let (u=omega-sqrt{q^2+m_s^2}>0), (v=|mathbf{k}-mathbf{q}|), and (p_c=sqrt{E_c^2-M^2}) at the highest archived parent energy. A parent at or below the cap can solve the past-spacelike radial equation only if
`v p_c - u E_c - (M^2 - m_f^2 + u^2 - v^2)/2 >= 0`.
The expression is strictly increasing in (v). Its positive root is computed with the rationalized quadratic; at every partially open case an independent parent-energy inversion reproduces (E_c). A 4,096-point direct cosine classification checks the normalized geometric measure, and independent fixtures cover fully closed, partially open and fully open caps.

The reported fraction is of **cosine measure**, not of rate or population. It is only the permission from this parent cap on this branch; other cuts and occupations can still remove weight.

| Chart | k=.00140 cap fractions at q=(1,3.5,15)e−6 | k=.00145 cap fractions | k=.00150 |
| --- | --- | --- | --- |
| Fine D1 | 1, 1, .6760 | 1, 1, 1 | 1, 1, 1 |
| Fine D2 | .2605, .4331, .4910 | 1, 1, 1 | 1, 1, 1 |
| Fine D8 | 0, 0, .0810 | 1, .8853, .5965 | 1, 1, 1 |
| Third D8 | 0, 0, .0199 | .9431, .6281, .5365 | 1, 1, 1 |

The sampled transition to all cosine angles at q=1.5e−5 occurs near k=.001490853 on fine D8 and .001496962 on third D8 (common parent cap (E_c=32.000039)). This is close to the previously fixed k=.00150; assuming that all angles remain available when moving to smaller k is not justified. Detected crossings are sign changes on the sampled .00135–.00165 scan, not a claim of uniqueness or exhaustive phase-space support. No cold scalar production coefficient is integrated here.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_early_pole_endpoint_past_k_cap_aperture
p=results/early-pole-cell-endpoint; d=2026-11-28
python -W error early_pole_endpoint_past_k_cap_aperture.py \
  "$p-fine"{,-d2,-d8}"-x05-state-$d.json.gz" "$p-third-d8-x05-state-$d.json.gz"
```
