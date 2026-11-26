# Matched q32/q64 support on a short coupled warm continuation

The [expanding prompt-only support test](INDEPENDENT-WARM-PROMPT-REENTRY-EXPANDING-2026-11-26.md)
left out both below-pole elastic scattering and signed scalar pair
conversion. Repeat the support comparison with those two operators
enabled, starting from the same actual q32 checkpoint. Widen the
represented grid with five additional initially empty shells up to
twice the q32 momentum edge; do not copy known collisionless
cohorts into those bins. Use the same stored inherited exterior
gravity in both runs. New positive pair births beyond each
expanded domain are recorded separately as collisionless receipts
by the original integrator, not returned as scattering inputs.

Run matched two- and four-block continuations to (b=1.001)
and visible temperature 10.380162 MeV. The archived collision
energy/number ledgers and independent retained-plus-new-receipt
number checks close to (10^{-14}). At four blocks, the new
explicit q32-to-q64 shell populations contain about
(1.144	imes10^{-11}) of the initial scalar number and
(2.213	imes10^{-11}) of the initial light number. They carry
(3.03	imes10^{-11}) and (1.30	imes10^{-10}) of the
initial total dark energy, respectively. These are distinct from
new collisionless births at the **outermost** cut.

The four-block shared-shell scalar/light number (L^1)
difference is (4.38	imes10^{-8}/6.60	imes10^{-9}) normalized
to each initial species count. The corresponding full energy
changes relative to initial total energy are (6.14	imes10^{-9})
and (5.87	imes10^{-9}). Sampled light occupations in the
shared shells differ by at most (5.34	imes10^{-13}) absolute.
The time-discretization change is larger than the tested support
difference in the aggregate moments: refining two to four
blocks moves the scalar/light count (L^1) by about
(5.21	imes10^{-5}/3.08	imes10^{-5}). The support comparison
also rebuilds the event quadratures; its difference is not a
proof of the exact physical effect of the added tail alone.

These are short selected-operator runs with the original
closed prompt geometry at each chosen outer cutoff. No older
unknown exterior phase distribution, heavy dynamics, complete
visible source, or future abundance is reconstructed.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_matched_explicit_support_coupled_reference.py EOS.dat run2 --execute
python -W error warm_matched_explicit_support_coupled_reference.py EOS.dat run4 --steps 4 --execute
python -W error warm_matched_explicit_support_coupled_reference.py EOS.dat run2 --refined-work run4
~~~

The [radiation table](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat)
digest is retained in the archived plans. [Matched metrics](results/warm-matched-explicit-support-coupled-2026-11-26.json.gz);
exact final checkpoints: [q32, two blocks](results/warm-matched-explicit-support-q32-2-blocks-2026-11-26.json.gz),
[q64, two blocks](results/warm-matched-explicit-support-q64-2-blocks-2026-11-26.json.gz),
[q32, four blocks](results/warm-matched-explicit-support-q32-4-blocks-2026-11-26.json.gz),
[q64, four blocks](results/warm-matched-explicit-support-q64-4-blocks-2026-11-26.json.gz).
