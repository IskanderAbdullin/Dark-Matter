# Initial net number loss can converge while hard-state feedback does not

The expanding [angular threshold comparison](INDEPENDENT-EARLY-POLE-ANGULAR-CHEMICAL-CONTROL-2026-11-26.md)
puts particles into a narrow band near q=√3 m/T. Its numerical
width controls how many available states share the same number
of products. At fixed cold populations and pole, linearize the
two hard Bose factors and their reverse reaction with respect
to the hard cell occupations. For four high-occupancy cold inputs,
the leading hard quadratic terms cancel, but the linear terms
still increase production as hard states fill.

Compare the representative-node Jacobian with an independent
finite-volume integration. The latter splits the exact outgoing
energy interval at all receiving shell faces and their reflected
partner faces, uses stepwise hard occupations, and integrates the
outgoing angular measure analytically. The physical outgoing
pair energies sum to the incoming energy exactly. This second
operator is a frozen feedback diagnostic; it does not evolve
representative-node energies.

For the .06 third-grid input at x=0.05 the finite-volume **initial
formal loss** stays about 69.15 per log a as the hard interval
is resolved. The largest local formal bare hard Jacobian
eigenvalue increases from 1.42 × 10¹⁰ to 9.36 × 10¹¹, and then
1.04 × 10¹² per log a at the tested 300 and 700 inserted-cell
partitions. Changing the incoming angular quadrature from 32
to 128 at 300 cells changes the latter comparison by about
2%, still smaller than the observed partition effect. These
values describe the formal collision graph only; they are not
a measured inverse lifetime or an in-medium rate.

A separate positive representative-node stencil splits outgoing
intervals at the **nodes** instead of faces, then adds Gauss
points within each segment. This keeps exact projected energy,
stoichiometry and thermal balance while reducing sampling of
only a few outgoing angles. On the frozen input, one to four
Gauss points give the same initial rate at the tested resolution;
that does not independently establish the later hard feedback.

As a separate width sensitivity, add an illustrative energy-sharing
half-width in quadrature to the on-shell angular width while keeping
the sum of outgoing energies fixed. This is not a spectral-function
calculation and need not satisfy exact on-shell momentum kinematics.
At 700 inserted hot cells, increasing the extra width from zero to
0.01 m, 0.03 m and 0.10 m changes the hard-feedback eigenvalue from
1.04 × 10¹² to 2.72 × 10¹¹, 9.56 × 10¹⁰ and 2.91 × 10¹⁰ per log a.
Initial formal net loss does not change on this snapshot because
the broadened products still sample essentially the same initial
hard occupancy. These stable nonzero-width trends do not justify
assigning any of these widths to the actual scalar medium.
[Width diagnostic data](results/early-pole-chemical-hot-width-sensitivity-2026-11-26.json).

~~~sh
python -W error early_pole_chemical_hot_feedback_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
python -W error early_pole_chemical_hard_shell_feedback_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
python -W error early_pole_scalar_number_change_interval_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
python -W error early_pole_chemical_hard_shell_feedback_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz --broadening-study
~~~

[Raw feedback and interval diagnostics](results/early-pole-chemical-hard-feedback-2026-11-26.json.gz).
