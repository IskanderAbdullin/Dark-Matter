# 28. Orbit-invariant conservative phase remap
2026-11-25; resumed and independently rerun at 11:20 UTC, base 11321f144. Frozen potential, bound orbits only. This recovers the unsaved last experiment from the preceding launch.

For fixed binding s=ψ−q² and squared angular momentum L=2r²(q²−p²), both invariants are exactly constant along frozen characteristics. The normalized spherical measure r²dr dE dp/4 (with E=q²,p=q_r) equals ds dL dt/(8√2) on each radial branch. Full-period phase θ∈[0,1) has measure ds dL T/(8√2)dθ. The code uses a common factor four larger, which cancels in all normalized ratios. Integrate turning-point flight time with r=rp+(ra−rp)sin²z and map by monotone PCHIP. Streaming is translation by Δt/T. Exact overlap of positivity-limited linear reconstruction on nonuniform periodic phase cells conserves each orbit's mass, hence total N, binding and L. No artificial orbital mixing or exterior sink. Unbound orbits and changes of potential not included. Exponential phase mesh k=8 resolves the narrow inner portion of long weak-binding orbits. Roundoff clipping monitored separately.

Paper 26 weak-binding Gaussian patch, independent first-crossing CDF .8901858457653883 at t=10:
| nL | phase cells | flight table nz | steps | outside fraction |
|---:|---:|---:|---:|---:|
|128|128|512|1|.8102499978404|
|128|256|512|1|.8505816717088|
|128|512|512|1|.8847723653216|
|256|1024|1024|1|.8896395585468|
|256|2048|1024|1|.8901262676864|
|256|4096|1024|1|.8901790330851|
|256|2048|1024|20|.8895600138460|
|256|4096|1024|20|.8901155500008|
|256|4096|1024|100|.8900038642445|
Fine one-step error 0.000765% relative; 100 steps error 0.0204%, showing repeated remapping diffusion despite exact invariant ledgers. Max aggregate N and L drift 2.3e−14; minimum before roundoff clipping zero. Uniform phase 512 cells fails (.6487582) since weak-binding orbital period is ~3503 dimensionless, much longer than the inner transit. Nonuniform grid is essential here, not a universal guarantee.

Independent oscillator, 20 remaps to t=1, phase n=128/512/2048 gives weighted L1 .00569210/8.91797e−5/1.84434e−6; N error ≤5.6e−16. No new claim about physical halo evolution, central singularity, nonlinear self-gravity or condensates. CI not run. Reproduce: python lowered_bose_orbit_remap.py --nL 256 --nphase 4096 --nz 1024 --steps 100; --oscillator --nphase 2048 --t 1 --steps 20.

Next: conservative projection between local collision cells and orbit cells without remapping the whole distribution at each collision step; quantify momentum and entropy limitations. Preserve energy binding exactly on frozen potential, handle unbound open characteristics separately.
