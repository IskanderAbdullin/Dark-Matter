# Preserve charge steps without independently evolving thin face volumes

The [two-width replay](INDEPENDENT-EARLY-POLE-BOUNDED-FACE-LONG-WIDTH-2026-11-28.md) found that making separately evolved faces thinner greatly changed the live result, even when the initial static source and the time steps were accurate. Here use two **cell-owned endpoint nodes** instead. An original charged cell can have a different occupation on either side of an interface even though the two touching endpoints have the same momentum. Interpolation is only between the two nodes *owned by the containing cell*, never across its neighbour.

At the archive reference time give the two nodes positive number volumes whose sum is the physical cell volume, and whose energy-weighted sum reproduces the exact first energy moment of that cell. Hold those node volumes fixed during comoving evolution, while their physical energies change with the scale factor. This avoids arbitrarily small face volumes and keeps the number invariant well defined. In each reaction use the exact allowed continuum energy aperture, with two-node deposition and geometric occupation/hole factors for both charged legs. A uniform initial occupation within the original cell reproduces its exact step. Node-wise Pauli, paired number and projected energy constraints are still checked; this is a quadrature representation, not a claim that finite probability resides at a physical interface.

An independent all-parent-bin × all-daughter-bin step integral uses the expanded long-double aperture endpoints and thus does not reuse the sorted quadrature chart. At the initial snapshot its scalar mode and gross coefficients agree with the projected chart to 2.83e−14. The initial growth clocks are 4.41826567 and 3.67048175 in the two archives and agree across the nested charge subdivisions; the scalar cells are unchanged. Maximum normalized event-energy and common-temperature detailed-balance deviations are below 5e−15.

For a local evolving interval Δlog(a)=5e−5, independently refine only the charge intervals by splitting each original cell into 1, 2, 4, 8 or 16 parts:

| Archive | Scalar increment, no subdivision | Scalar increment, 8 | Scalar increment, 16 | Last relative spatial change |
| --- | ---: | ---: | ---: | ---: |
| Fine | 3.8225372e−6 | 3.8060737e−6 | 3.8053432e−6 | −1.92e−4 |
| Third | 3.1336713e−6 | 3.1225396e−6 | 3.1212762e−6 | −4.05e−4 |

The last 32→64 time change is at most 1.19e−7. The independent normalized energy continuity residual stays below 1.93e−16. The large artificial face-width parameter is absent; actual spatial resolution is still a distinct error source. These local tests constrain neither earlier cumulative production nor a full new archive evolution, and the default driver is not changed.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_endpoint_direct_aperture_audit.py
for g in fine third; do
  python -W error early_pole_cell_endpoint_collision_pilot.py "$g" --levels 1 2 4 8 16
done
```
Compare the [all-pairs direct check](results/early-pole-endpoint-direct-aperture-2026-11-28.json) and [local time/charge mesh receipt](results/early-pole-cell-endpoint-collision-2026-11-28.json).
