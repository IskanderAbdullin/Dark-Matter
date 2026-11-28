# Closed FD collision step with the tabulated visible equation of state

Replace the ideal visible heat capacity by the thermodynamically consistent interpolation of the [Saikawa–Shirai authors' central table](https://member.ipmu.jp/satoshi.shirai/DM2020/). During the step, determine visible temperature from the unused energy under the fixed total-energy constraint; resolve the collision rate at the actual state corresponding to the averaged entropy forces. This is the vacuum timelike charged-lepton FD test only. No expansion, elastic heat, medium corrections, or abundance calculation is claimed.

At three finite displacements, the error-halving ratios were 4.0016–4.1102 against an independent high-accuracy differential integration. The maximum scaled discrete entropy identity error was 1.92e-13. The original ideal-reservoir results reproduce byte for byte with the generalized solver. As a separate equilibrium check across 13 temperatures, the derivative of the interpolated visible energy reproduces its tabulated-EOS heat capacity to 1.73e-8 in relative units; two collision quadratures agree on both mode rates to 2.11e-9.

The ideal proxy exaggerates the finite-capacity change in the fast rate over the sampled range. At T=0.68825 GeV the tabulated bath raises it by 5.60% above an infinite reservoir, versus 35.77% with the ideal 10.75-degree proxy. At T=0.05 GeV the tabulated rise is 27.67%, versus 40.09% with the proxy. These are local rate changes, not the effect on final relic density.

## Reproduction

Save the authors' [central numeric table](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat) locally as `eos2020.dat`; it is not redistributed. Normalize the held interval to compute the same input hash:

```sh
python - <<'PY'
from pathlib import Path
r=[s.split() for s in Path('eos2020.dat').read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
out=[' '.join(x) for x in r if 0.008<=float(x[0])<=1.2]
assert len(out)==102 and all(len(x.split())==3 for x in out)
Path('eos-range.dat').write_text('\n'.join(out)+'\n')
PY
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error fd_portal_tabulated_bath_discrete_entropy.py eos-range.dat
python -W error fd_portal_tabulated_capacity_modes.py eos-range.dat
```
