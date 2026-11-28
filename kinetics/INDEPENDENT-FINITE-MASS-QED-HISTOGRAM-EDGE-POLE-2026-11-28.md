# A histogram edge adds a logarithm to the local spacelike pole limit

The [smooth two-pole slices](INDEPENDENT-FINITE-MASS-QED-CUT-LOCAL-SEQUENTIAL-SUBTRACTION-2026-11-27.md) explicitly assumed no occupation edge at a root. Here place the first parent-side charged propagator root exactly at a genuine archived step by varying the parent azimuth. On the fine archive it crosses a daughter-occupation edge while the parent cell is unchanged; on the third archive it crosses a parent-occupation edge while the daughter cell is unchanged. The second propagator zero stays outside the diagnostic interval.

Work in x=D1/(2 q·K) so the other propagator is proportional to 1−x. Let w⁻ and w⁺ be the one-sided forward Pauli products divided by their mean, f the smooth first diagonal weight (including dE/dx and a fixed normalization), and C the smooth mixed spinor numerator. For both equal diagnostic imaginary widths ε, the subtraction of the positive diagonal Lorentzian area is
```text
π f(0) [w⁻ + w⁺] / (2 ε).
```
This alone is not sufficient. The linear part of the diagonal and the simple real interference pole produce a logarithm with coefficient
```text
(w⁺ − w⁻) [f'(0) + 2 Re C(0)].
```
Indeed x/(x²+ε²) has opposite logarithmic half-integrals. The constant diagonal tails, finite mixed term and any imaginary numerator δ term are retained by the analytical moment formulas. The coefficient is defined with this explicit x variable, window, width convention and chosen Pauli weight.

| Archived step | w⁻, w⁺ | Log coefficient | Local value after area and log removal |
| --- | --- | ---: | ---: |
| fine daughter | .969923, 1.030077 | .1203223 | −1.9291535 |
| third parent | .765695, 1.234305 | .9373309 | −2.3522306 |

Two spinor interpolation orders reproduce the smooth numerator with the archived step held outside the polynomial fit. Independent direct finite-width integration agrees to 4.71e−12 in the tested moderate widths. The limiting logarithm is predominantly coherent interference in these examples; removing diagonal poles alone misses it. For both samples the simultaneous area-and-log removal approaches the declared local finite constant, not a cutoff-independent physical rate.

As a separate regularization diagnostic, shift the occupation boundary from x=0 to x=λε while keeping the same fitted numerators. The leading scaled area εI changes continuously from the left to the right weight as λ varies, as
`f(0)[w⁻(π/2+atan λ) + w⁺(π/2−atan λ)]`.
It agrees with the retained finite-width checks to 2.98e−6 absolute. Therefore taking the root/edge coincidence first is different from holding a displaced edge on one side as the width is removed. A real collision model needs a physically resolved distribution and matching; a step histogram does not supply either.

## Smooth occupation profiles do not select one common offset

Replace the sharp step by an explicitly artificial symmetric ramp of halfwidth σ, with linear or endpoint-flat cubic profile on |x|≤σ. Outside it the two original one-sided rates are unchanged. At fixed r=σ/ε the change for the linear ramp has the exact leading form
```text
log_coefficient × [1 − atan(r)/r − log(1+r²)/2].
```
It follows by integrating y(y/r−1)/(1+y²) in the scaled coordinate y=x/ε. The alternative cubic profile replaces y/r by 1.5y/r−0.5(y/r)³. Both independently integrated spinor-based controls reproduce their predicted joint limits to 8.20e−12 absolute. For fixed smoothing σ=.002 the area-subtracted integral saturates as ε shrinks, unlike the logarithmic sharp edge. As σ is then removed, the two profiles share the same logarithm but differ by a finite offset: cubic minus linear tends to one third of the log coefficient after that common logarithm is removed. At the actual finite ratio the independently computed offset agrees with its finite-r prediction within 5.00e−7.

Neither ramp is proposed as the physical distribution. This example shows why a histogram plus a free choice of smoothing and a constant imaginary width cannot by themselves provide a unique overlap prescription.

Run `python -W error finite_mass_qed_histogram_edge_pole_reference.py` from `kinetics/`, and compare the [machine receipt](results/finite-mass-qed-histogram-edge-pole-2026-11-28.json). No in-medium width, fully regulated Ward prescription, or matched scalar source is claimed.
