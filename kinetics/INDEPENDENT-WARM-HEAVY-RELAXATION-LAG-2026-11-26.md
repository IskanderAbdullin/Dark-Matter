# Check the heavy decay lag on a prescribed interpolation of saved states

The [frozen heavy energy diagnostic](INDEPENDENT-WARM-FROZEN-HEAVY-ENERGY-2026-11-26.md)
does not measure how a heavy mode follows a changing daughter distribution.
Use the [actual four saved warm blocks](INDEPENDENT-WARM-LIGHT-OCCUPATION-NEXT-MATCHED-STAGE-2026-11-26.md)
to define a separate driven model of that lag.

At each saved point and fixed heavy comoving momentum, compute
the stationary source (f_*=A/(A+B)) and linear relaxation in
log expansion (lambda=Gamma M(A+B)/(E_H H)).
Interpolate (f_*) linearly between the two endpoints of
each cell. Set its (lambda) to the smaller endpoint value and
start the first cell at the specified initial stationary
occupation. The interpolation and constant response are
**assumed**, not enclosures of the true evolution. Integrate
this scalar driven equation analytically, and verify the
elementary update against an independently integrated
moderate-stiffness test. The test agrees in absolute value
to (5.35×10^{-16}).

On the fixed comoving band up to 1 GeV, even the least
integrated response in a cell is about (2.7868×10^{10}).
At the q32 final state, the maximum absolute lag in the
sampled occupations is (8.014×10^{-25}), while the
absolute energy-weighted lag is (1.982×10^{-19}) of
the retained dark energy. The heavy stationary proxy
itself is (1.237×10^{-6}) of that energy, so the ratio
of the two is around (1.6×10^{-13}).
A second radial order checks the integral. These
numbers justify the approximation only for this
particular driven interpolation of the resolved
decay-only source; they do not bound short-lived
forcing, additional production, or an unknown
initial heavy population.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_heavy_relaxation_lag_reference.py EOS.dat \
 "$D"/warm-early-q{16,32}-time32-checkpoint-2026-11-26.json.gz \
 "$D"/warm-light-occupation-next-q{16,32}-tape-2026-11-26.json.gz
~~~
