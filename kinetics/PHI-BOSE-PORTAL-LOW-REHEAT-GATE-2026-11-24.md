# Low-reheating dilute lepton-portal injection and heavy-depletion gate (2026-11-24 UTC)

## Why test the opposite initial condition
The [new variational entropy gate](PHI-BOSE-ENTROPY-ENERGY-GATE-2026-11-24.md) disfavors an initially fully thermal, subsequently isolated 9-state benchmark. But a colder/empty dark sector can evade that entropy bound. Test an explicitly different branch: instantaneous reheat into SM at T_R=5..100 MeV, initially empty dark sector; only injection via charged-lepton annihilation through the same mV=10 GeV, gD=.1, epsilon=1e-5 vector with GammaV=.005305164626 GeV. Use fixed g*S=g*rho=10.75 and RD as controlled approximation at 5–20 MeV and illustrative beyond that. No thermal dark population assumption. This is forward injection only, **not** a physical coupled Boltzmann solution, and omits hadrons, direct scalar loop production, backreaction, QCD/reheat dynamics, finite-T masses/width and e± entropy evolution.

## Physical MB reaction density
For each charged lepton l mass m_l, each distinct Dirac chi mass m, standard spin-averaged s-channel cross section:
    sigma_l(s)= C*s/(12*pi*D(s)) * (1+2*m_l²/s)*(1+2*m²/s) * beta_chi/beta_l,
    C=(gD*eps*e)^2, D=(s-mV²)^2+mV² GammaV².
Multiplying by the exact incoming lambda(s,m_l²,m_l²)/sqrt(s) and degeneracies g_l=g_antil=2 gives the MB forward event density
    gamma_l(T)=C*T^8/(48*pi^5) * integral_{2*max(m_l,m)/T}^inf du
       u^6*K1(u)*(1+2*m_l²/s)*(1+2*m²/s)*beta_l*beta_chi/D(s)
with s=T²u². This counts one chi and one antichi per event. Shift the Bessel integral by the threshold and use kve to keep the exponential explicit. Independent analytic contact/massless limit integral u^6 K1(u)=384 is reproduced at T=.001GeV to 6.79e-7 (remaining finite-mV effect). Vary cutoff threshold+100→150 and verify finite electron mass rate <= massless e rate for both m at T=.005,.03,.1 GeV (12 local assertions passed). NWA pole at u=mV/T is highly Boltzmann suppressed at these low temperatures; no high-temperature inference from this quadrature.

At a prescribed SM T and vacuum kinematics, thermal FD f_l <= exp(-E/T) and >= exp(-E/T)/2; thus MB is an upper bound on the forward fermionic SM source (Pauli final factors <=1), and in the **empty final** limit FD forward is >= MB/4. This does not bound actual *net* production after dark backreaction; other initial SM processes are not included. The MB integral uses physical masses for e,mu,tau, not three fictitious massless lepton species.

Y_i per single charge = integral gamma_l/(s_SM H) dln(T_R/T), sum three leptons, upper cutoff ln(T_R/T)=5 (low-T tail exponentially suppressed in this branch). For each injected chi1 or antichi1, the one-pass decay chi1->chi2+phi yields one stable heavy and one scalar. Count both charges.

| T_R MeV | MB injection-only Omega_heavy h² | phi from one decay per chi1 only | phi if additionally all injected heavies convert one-to-one into phi |
|---:|---:|---:|---:|
|5|0.00327246|5.81e-20|1.40716e-8|
|8|0.0680109|1.46e-13|2.92447e-7|
|10|0.202833|1.98e-11|8.72203e-7|
|15|1.05826|1.39e-8|4.56444e-6|
|20|3.01785|3.77e-7|1.33539e-5|
|25|6.97312|2.81e-6|3.27902e-5|
|30|14.5121|1.10e-5|7.33730e-5|
|50|120.998|1.99e-4|7.18834e-4|
|100|1482.0694|0.00304077|0.00941367|

The MB injection-only heavy abundance crosses 0.120 at T_R=8.94495 MeV (conditional; **not** observational reheat bound or actual heavy relic). For T_R=20 MeV even 100% one-to-one conversion leaves scalar abundance ~9000 times below 0.120. This is **not** a bound for arbitrary number-changing dark reactions: in the full minimal network, phi phi -> chi1 antichi1 -> chi2 antichi2+2phi -> 4phi can multiply scalar count if energies permit. Such reactions, direct SM/hadron production, and thermalization must be separately solved or bounded. The last column assumes only one pass of chi1 decays and no further number multiplication, not a claim that physical conversion is efficient.

## A necessary depletion criterion for any one-pass heavy-fed scalar branch
Let N_H^in=N_chi1^in+N_chi2^in including both charges and N_H^f be surviving heavy count. With at most one phi per chi1 decay plus one phi per heavy removed (pair->two phi), N_phi^f <= 2 N_H^in-N_H^f. If the desired scalar rest-mass density fraction F of total dark density is at least F, then N_H^f/N_phi^f <= q=(1-F)/F*(mphi/m2), so N_H^f/N_H^in <=2q/(1+q).
For F=.5/.9/.99, survival fraction must be at most 8.599963e-6/9.555551e-7/8.686868e-8. In the decay-only branch without depletion, each scalar accompanies a stable heavy and Omega_phi/Omega_H <= mphi/m2 =4.3e-6. Thus small initial entropy is not itself a scalar production mechanism; an efficient, physically derived multiplication/depletion history is necessary. No universal model exclusion is implied.

[Source](phi_bose_portal_low_reheat_gate.py), [machine benchmark](results/phi-bose-portal-low-reheat-gate-2026-11-24.json). Local scipy 1.17.0, no Agent PC, no new CI claim. Next decisive calculations: quantify true FD/BE portal freeze-in + number multiplication through inverse phi-pair production and kinematic threshold; test viable scalar f(q) and structure suppression. Keep CLASS-ready=false.
