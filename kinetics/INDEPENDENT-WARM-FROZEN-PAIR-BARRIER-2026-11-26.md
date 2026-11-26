# Uniform frozen scalar-pair sign at a half-filled tagged light mode

Use the actual corrected q32 endpoint at 7.394049 MeV and the
accepted [common-chart q64 replay](INDEPENDENT-WARM-PAIR-COMMON-CHART-2026-11-26.md)
at 7.754805 MeV. Hold every bath mode fixed except a single tagged
light mode, whose occupation is tested at f_L = 1/2. This is a
directional test on the shape reconstruction in each archive,
not an assumption about the evolution of the bath.

Write z_phi = log[(1+B)/B] and z_L = log[(1−F)/F].
For a represented fermion partner and both scalars represented,
the eventwise gain/loss ratio at the tag test value is

    G/L = exp[−(z_phi(E₁) + z_phi(E₂) − z_L(E_p))].

If either scalar is in vacuum, the gain on that event is zero.
Each stored entropy is linear within its shell, so the bracket is
affine on every three-energy box. Minimize it on the box intersected
with the tag-energy slab m_L ≤ E₁ + E₂ − E_p ≤ E_L,max.
The extrema are box vertices or edges intersected by either plane.
This drops angular restrictions, which can only lower the gap.
The two minimum gaps are 10.8535 and 10.7063, giving inward-partner
ratios at most 1.94e-5 and 2.25e-5, respectively.

Outside the represented partner grid, a vacuum partner can still
be produced. Both occupied scalar inputs must then have
E_i ≥ m_L + E_L,max − E_phi,max and their sum must
exceed m_L + E_L,max. Use a separate two-dimensional affine
box bound on the scalar entropy sum and a constant upper amplitude
bound for this restricted gain. For the downward term, keep only
the incoming-cosine half-range [−1, 0]. A uniform positive
lower bound on the longitudinal squared amplitude there is
8y⁴ p_*² M²/A², with A = M² − m_L² − m_phi² + s/2.
Take each shell's most conservative (p_*, A), flux and partner
occupation, and add these disjoint positive loss rectangles.
The resulting upper exterior gain divided by the lower inward loss
is 5.28e-9 on q32 and 2.93e-20 on q64. Thus the *scalar-pair* term
has a strictly downward bound over the entire represented tag
momentum range on these two frozen reconstructions, without
excluding the positive vacuum-final channel.

An independent tagged quadrature verifies the sign at 26 q32 and
31 q64 momenta. On q32 the open stationary-parent prompt surrogate
can point upward by 3.44e-8 per Hubble time at rest; the
represented prompt is exactly absent there, and the pair term is
negative and much larger (−0.229 per Hubble time). The sampled
pair-plus-prompt signs are downward. This sampled statement does not
supply an all-momentum inequality for prompt, test the elastic term,
or control untracked incoming occupations. In particular, it does
not establish a persistent light/parent cap or validate a bound that
assumes such a cap over a live cosmological trajectory.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_pair_frozen_continuum_cap_bound_reference.py \
 results/warm-multishell-trajectory-lifted24-q32-2026-11-26.json.gz \
 results/warm-pair-common-chart-q64-terminal-2026-11-26.json.gz
python -W error warm_frozen_light_cap_components_reference.py EOS.dat \
 results/warm-multishell-trajectory-lifted24-q32-2026-11-26.json.gz \
 results/warm-pair-common-chart-q64-terminal-2026-11-26.json.gz
~~~

[Uniform pair bounds](results/warm-pair-frozen-continuum-cap-2026-11-26.json)
and [sampled prompt/pair directional rates](results/warm-frozen-light-cap-components-2026-11-26.json.gz).
