# The light fermion's innermost shell is sensitive to far-tail scalar support

A [sharp shell regularity test](INDEPENDENT-WARM-SHELL-LIPSCHITZ-CAP-2026-11-26.md)
locates the smallest slope allowance for certifying (f_{m light}<1/2)
in the innermost light shell. Inspect the (Hleftrightarrow L+phi)
contribution directly at fixed light momentum. This is a frozen
diagnostic, not an evolution or proof of that cap.

In the light rest frame the required scalar energy is
(E_*=(M^2-m^2-m_phi^2)/(2m)), about 0.240 GeV. With physical
light momentum (p), the allowed lab scalar energies have center
(E_L E_*/m) and half-width (p k_*/m), where
(k_*=(M/m)P_*) and (P_*) is the scalar momentum in the heavy
rest frame. Evaluating an angular average instead of dividing
a small numerical energy interval by (p) gives the exact finite
rest limit. For equal per-charge, per-spin degeneracies of (H)
and (L), the collision term is
[
 C_L(p)=Gamma_H(M/m)^3(m/E_L)
        leftlangle
           (1-f_L)f_H(1+f_phi)-f_L(1-f_H)f_phi
        ightangle .
]
The average is uniform in allowed scalar energy for fixed (p);
(E_H=E_L+E_phi). No scalar state is filled outside the archived
support. The rate is normalized by the proper heavy decay width.
As independent checks, thermal detailed balance cancels pointwise
within (1.2	imes10^{-15}), and directly integrating the cold
kernel with initially empty daughters agrees with the isotropic
parent-coordinate calculation of daughter number and energy
to (9	imes10^{-14}).

At the first matched saved temperature, 10.391 MeV, the
innermost shell reaches only (4.72	imes10^{-5}) GeV in
light momentum. Nevertheless, its scalar recoil window is
0.239434–0.240567 GeV, well beyond the q16 scalar edge
(0.151125 GeV). The inverse term and stimulated part vanish
in q16; a frozen heavy population still gives spontaneous gain.
In q32 the window is present with scalar occupation up to
(6.11	imes10^{-15}). Its inverse loss per Hubble time,
(2.3893	imes10^{-3}) at rest, slightly exceeds spontaneous
gain (2.3682	imes10^{-3}). Thus the net q32 contribution
is negative at every probed point in this shell at all six
saved times, whereas q16 is positive. These are local rates
in occupation per Hubble time, not integrated yields.

The sign flips to positive when cropping only the scalar
support of the **same** q32 reconstruction to q16's cutoff,
both when holding its former stationary heavy shape fixed and
when recomputing it with the crop. This isolates the support
effect without substituting the q16 light spectrum. The q32
frozen local relaxation is about 41 Hubble rates at rest in
the first frame, with equilibrium (f_L=5.76	imes10^{-5}).
For its scalar window, a sufficient uniform heavy occupation
bound to make this operator nonpositive at (f_L=1/2) is
(ar f_phi/(1+2ar f_phi)simeq5.76	imes10^{-15}) at
rest. This follows from
((1-f_L)f_H+(f_H-f_L)f_phi). The frozen parent is lower than
that across the sample, but this does not bound the true parent.

Only the selected decay and inverse decay with a frozen
heavy stationary response were probed. No dynamical heavy
profile, exterior scattering or elastic gain is inferred.
Tiny occupation is not a safe criterion for dropping a
kinematically selected high-energy tail from a local rate.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error warm_near_rest_light_inverse_decay_reference.py EOS.dat
~~~

Use the same [tabulated radiation EOS](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat)
as the warm restart; its SHA-256 is in the
[diagnostics](results/warm-near-rest-light-inverse-decay-2026-11-26.json.gz).
