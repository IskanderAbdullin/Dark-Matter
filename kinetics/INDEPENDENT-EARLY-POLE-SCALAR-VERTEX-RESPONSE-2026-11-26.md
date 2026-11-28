# An infrared vertex-resummation sensitivity check

The [reversible elastic control](INDEPENDENT-EARLY-POLE-SCALAR-ELASTIC-STENCIL-2026-11-26.md)
evaluates bare scattering on a highly occupied scalar distribution.
As a sensitivity analysis, compute the isotropic nonrelativistic
retarded response from [Walz, Boguslavski and Berges](https://arxiv.org/abs/1710.11146),
Eq. (26), for the infrared portion of the same snapshots. The
nonrelativistic test uses q<10⁻⁵, well below the Hartree pole. It
uses g=λ/(8m²) from matching the real-scalar threshold quartic
to the nonrelativistic contact convention. The literature
large-N expansion is **not** controlled for this real N=1 model.

For piecewise constant occupations, integrate the logarithmic
real part per shell with its elementary antiderivative, splitting
the imaginary part at its two exact bounds. An independent
adaptive quadrature checks the real part at the selected angles,
and the response at opposite frequency is the complex conjugate.
Over transfers 0.25–4 times the peak in q²f and frequencies
0–2 times P²/(2m), the .06 third-grid x=0.05 response magnitude
is roughly 1.4 × 10³ to 2.6 × 10⁴. At P equal to that peak and
ω=P²/(2m), the squared vertex ratio to the bare value is
1.47 × 10⁻⁹ (2.09 × 10⁻⁹ on the previous grid).

As a distinct frozen control, integrate one transfer channel
with its actual incoming/outgoing kinetic-energy difference and
azimuth over the same closed IR graph in both directions.
Do not carry a single constant vertex ratio across the operator.
The third-grid cold shape flux at x=0.05 for the .06 case falls
from 4.14 × 10¹⁵ to 9.37 × 10⁶ per reference H, a reduction of
about 4.4 × 10⁸. The intermediate momentum grid gives
3.04 × 10⁷ after the same change; two angular orders on the
third grid differ by 7% in the screened result. These are
**formal** shape fluxes, not scalar number production, a
particle-width measurement, or a reliable bound on relaxation.

The cold kernel is a transfer-dependent large-N-inspired proxy,
not the complete in-medium real-scalar scattering amplitude.
It omits the hot sector, other diagrams, self-consistent response
during redistribution, possible condensate, off-shell effects and
number-changing collisions. The comparison demonstrates that
medium treatment matters greatly; its residual flux magnitude
must not be read as a physical timescale. The graph keeps its
conservation and entropy checks after this weight change.

~~~sh
python -W error early_pole_scalar_nr_vertex_response_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
python -W error early_pole_scalar_cold_screened_flux_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Polarization samples](results/early-pole-scalar-nr-vertex-response-2026-11-26.json.gz);
[cold projected comparison](results/early-pole-scalar-screened-cold-flux-2026-11-26.json).
