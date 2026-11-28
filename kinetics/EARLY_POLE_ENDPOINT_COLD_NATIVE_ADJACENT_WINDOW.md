# Direct extension of the native D2 cold band to q=1.8e-5

On the *same* held fine and third D2 x=.05 endpoint states, add the separate band [1.5e-5,1.8e-5] with no extrapolation in momentum. The two independent meshes agree on the scalar-pole optics, exact state hashes, panel topology, forward/inverse excess closure and total ledgers. The added band contributes 2.28e-7 and 6.59e-7 of the expanded forward counts on the fine and third charts, respectively. It contributes up to 1.284e-5 of the forward centered excess.

On the expanded retained band [1e-6,1.8e-5], the third chart still has lower mean q, RMS q and mean speed at common mass for all three weights. With each chart's native pole mass instead, forward and inverse mean-speed intervals still straddle zero; the third-minus-fine excess-weighted speed remains negative, in [-3.91795e-4,-1.31354e-5]. This is a bounded numerical comparison on a specified momentum band, not a statement about speed or abundance outside it.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28
python -W error -m unittest -q test_early_pole_endpoint_cold_native_adjacent_window
python -W error early_pole_endpoint_cold_native_adjacent_window.py --root results \
 --archive "results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json" \
 --predecessor "results/early-pole-cold-native-placement-eight-$d.json" \
 --accepted "results/early-pole-cold-native-moment-sign-refinement-$d.json"
```
