# Keep the far-away resonance out of sub-pole thermal energy bins

The [global-width vacuum envelope](INDEPENDENT-WARM-LATE-COLD-VISIBLE-DOSE-2026-11-26.md)
allows every thermally occupied energy to use the largest
vacuum cross section, including a 10 GeV pole far beyond
most pair energies. Keep the same three charged-lepton
source, recoil threshold and entropy clock, but partition
the total incoming lab energy. For any energy bin ending
at Q_hi, s≤Q_hi². The factor s/[(s−V²)²+(V Gamma_V)²]
in the conservative cross-section upper increases until
s=sqrt(V⁴+(V Gamma_V)²), so bound it only on the bin's
allowed invariant-mass range. Integrate each bin's
Boltzmann polynomial tail using the difference of its
two incomplete gamma functions, without subtracting
underflowed ordinary probabilities. The fixed-bin integral
upper is monotone in T along the adiabatic history:
cooling raises the total recoil threshold and decreases
the Boltzmann factor at each held total energy. The binned
upper can therefore keep the earlier one-sided time rule.

In a separate conditional column, also apply the maximum
cold daughter fraction over open unpolarized decays.
For a cut p_c below the rest-frame recoil P, the maximum
evaluated at gamma=E_s,rest/E_cut simplifies to
p_c²/[2P(P+sqrt(P²−p_c²))]. This increases with p_c,
so the hotter endpoint is sufficient on each time bin.
No angular ceiling is assumed in the kinematics-only column.

| Starting T (MeV) | Prior log10 bound / held yield | Energy binned | Also open decay-angle cap |
| --- | ---: | ---: | ---: |
| 7.39405 | −118.89 | −126.59 | −130.52 |
| 10.442 | −50.38 | −58.38 | −62.01 |
| 12 | −33.69 | −41.79 | −45.31 |
| 25 | +6.04 | −2.61 | −5.48 |
| 50 | +14.09 | +5.08 | +2.82 |

All bins again end at .1 MeV and use the same held q64
cold yield for comparison. Starting at 25 MeV, the energy
bound alone is below 0.25% of that yield; with the angular
assumption it is below 3.4e−6 of the held yield. The 50 MeV
upper remains too loose to decide the contribution. It is
not a lower or a measured source. No stimulated decay,
hadronic production, prior heavy population, or scatter-in
from other momenta enters this ceiling.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_late_cold_visible_offpole_dose_reference.py EOS.dat \
 results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
~~~

[Energy-panel and angular ceilings](results/warm-late-cold-visible-offpole-dose-2026-11-26.json).
