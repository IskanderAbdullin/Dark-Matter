# Direct pole: metric feedback, biased slow trigger, and a conditional new candidate

Date: 2026-11-25. Base revision: `14999f922ab393a94099c76e24839ff45c478be2`, branch `research/early-bootstrap-20261125`. This is an independent pilot and a prospective low-energy model modification, not a cosmological detection, fit to all data, UV completion, or universal exclusion. Existing research code and candidate parameters were not modified.

## Decision and negative results

The original fixed homogeneous background does not justify treating eta as smooth. In an independent Newtonian-gauge calculation with the signed varying-mass collisionless eta distribution, self-consistent linear metric feedback, and standard free-streaming radiation fraction .409, a seed early potential psi=1e-5 gives |delta n_eta/n_eta|=.1 at z=3.39075 for k=.03 Mpc^-1. The fine direct-evolved potential run agrees with the algebraic metric constraint at relative 1.25e-8. Rest-frame number contrast is -.09906 at that threshold (not a pure gauge effect); eta background fraction is 9.64e-5. At initial amplitude 1e-8, the threshold is still reached at z=2.60365. A no-late-metric control, a radiation-fluid control, seed half/double, a sudden photon-baryon decoupling calculation (threshold z=3.37107), and independent grids were performed. The tested unstable band includes k=.003--.04 Mpc^-1 (eta or chi diagnostic threshold); .05 and .06 do not cross the 10% thresholds. Four even potential shape changes do not fix the instability. Simply increasing v by 2--3 times and retuning the initial seed shifts rather than eliminates growth; larger initial seeds also destroy the early contact proxy. The result is not an observational exclusion: the nonlinear fate of a subdominant component is not solved.

## New candidate: a weak explicit bias and slower field

Use the earlier renormalized vev-curvature potential `U0(x)=X[(x^2-1)^2/4+CW*(x^4 log(x^2)-1.5*x^4+2*x^2-.5)]`, with x=chi/v, add `epsilon*X*(x-1)^2`, and use v=10 v_ref. Keep X=2.036548001081445e-10 eV^4, v_ref=9.049809269161514e24 eV, m_eta(x=1)=.012664880357500833 eV, T_eta,0=.0039033442679353545/36 eV, CW=.2, seed=1e-8; retune epsilon=.0032765520306085277 for first x^4=.5 at z~3.8. This is an explicitly biased EFT, not a symmetry-preserving modification. The linear spurion and even quadratic shift are prescribed renormalized inputs. The positive vacuum is unchanged; a negative metastable minimum remains. Larger inertia reduces the rate of changes of the eta mass and the scalar curvature in H units.

- Max extra fraction z<=100 = .0187981; max H/H_LCDM-1 = .00953366; H0/H0_ref=1.00307119. Background continuity residual ~1.18e-5 on the stated sampled integral. No standard cosmological parameters refitted in these results.
- Halo A4 exposure (fixed illustrative halo since z=2) remains .99337 of reference. Retuned sigma0=.22236922 cm^2/g for the same initial Rayleigh e-fold diagnostic, not a demonstrated halo solution. Late contact proxy .00837964, early contact proxy z=100--1e5 3.60e-19 for seed 1e-8. All-scale finite two-loop ratio ~.7540; the matching scales are Q/M=.5,1,2; this does not control high-scale counterterms.
- An independent NR conservative quantum elastic solve on both reconstructed 220-bin production PSDs gives tail changes at q>2 q_peak = -.007786538 and -.007790657. Four ordered time intervals at tighter tolerance on mesh 240 give -.007786546. The independent old-background solve reproduces archived -.008102175 to 1.04e-9. Number and energy drift <=6.67e-16 and <=2.23e-16; detailed balance residual <=6.34e-16 relative to dilute scale; greatest accepted negative number fraction <=1.1e-16. Only homogeneous z<=100 AA elastic scattering, not chemistry or in-medium scattering. The spectra are transcribed from the public archive, with geometric face reconstruction and weighted number checks. The prior kernel is independently reproduced, not merely reused.
- Exact eta angular/radial evolution starts at z=100 (max earlier |m_eta|/T_eta=1.96e-5). Adiabatic-like initial psi=1e-5. On k=.0003,.001,.003,.006,.01,.02 Mpc^-1 the pressureless-matter/free-radiation pilot gives peak eta number contrast .00005395,.00027131,.00228528,.00527970,.00662622,.00938104; none hits .1 chi or eta. Acoustic approximation with instantaneous z_dec=1090 yields max eta .00826098 at k=.03, approximately .00791 at .05 and .00722 at .06; fine .03 run agrees within 2.3e-7 in eta contrast and density constraint <5e-14. Matter perturbations are <.185 in these runs. These are finite modes and simplified standard matter (the 86 eV relic represented as pressureless matter); no full recombination, polarization, nonlinear evolution, or full momentum-dependent dark matter transfer.
- Independent frozen collisionless root scan 480x100, then 720x160: post-half max integrated positive root .22090/.22040, versus 13.514 original; max additional eta feedback contribution .03959/.04282. No fast subhorizon roots meeting w=k/(aH)>=10, gamma/H>=5, gamma/H>=3|d ln |x|/d ln a| on either biased grid. Frozen roots do not establish actual growth or stability; direct mode integrations are primary.
- Compensated eta-number and absolute chi-x trial modes at k=.003, initial amplitude 1e-5, max eta contrast 2.17e-6 and 6.80e-6. The latter has relative density constraint residual 5.46e-6 because the physical metric is tiny; these trial modes are not an inflationary spectrum. For larger initial seed, homogeneous results demonstrate sign selection but not a domain-wall network calculation.

## Early bound and matching caveats

Finite-field alternative (constant A quartic and small S quartic): kappa=g0/v falls by ~10, so the earlier constructive tree bound lambda_A>=9(kappa chi_cap)^4/(4 lambda_S M^4) gives about 1.68e-21 for chi_cap=.1 reduced Mpl, lambda_S=3.41708e-6. Peak chi~.0466 Mpl is inside the cap. The global quartic obstruction remains in the earlier restricted positive-portal test. At fixed x, vacuum scattering amplitude for positive equal masses and subthreshold S is bounded by lambda_A + g0^2*x^2*(1/detuning+2/mS^2). Reuse the previous conditional phase-space rearrangement envelope, with F=2.1 f_peak, N<=2 N_ref r^3, m>=86 eV and H>=.5 H0 sqrt(Or) r^2. The bound integrates amplitude(r)^2 min(a^2 r,b^2/r) d ln r (a=2N_ref/m, b=F p_star^2/(4 pi^2), p_star=(6 pi^2*2N_ref/F)^(1/3)). It is a positive tree-vacuum bound, not an in-medium pole calculation. For seed 1e-8, L1/N_ref <=1.332e-4 and 1.325e-4 for r>=101; positive mode gain <=7.69e6 through the previous finite production start. At seed 1e-7 the bound is ~1.59e-4, at 1e-6 it is 1.46e-2, and at 1e-5 it is uninformative (~120). Failure of this conservative bound is not a failure of the model. Constant quartic and square stabilizer are separate matching alternatives; the square-equivalent constant quartic is not a viable early interaction. Caps, absence of primordial S/coherent fields, vacuum detuning, production, and high-scale matching remain unresolved.

All 12 homogeneous seeds from -0.1 to +0.1 eventually select positive x in this fixed background. For |seed|=.01, early z=100--1e5 contact proxy ~1.18 (not a lower bound). The bias does not automatically solve wall cosmology. The negative local vacuum x=-.99337, barrier x=-.03291, vacuum difference 2.66e-12 eV^4, unbiased-wall tension 1.02e21 eV^3. The dimensional pressure-balance estimate H/H0=1.81 (z~1.02 in fixed LCDM), one-wall-per-Hubble-volume fraction ~.022 at balance. Actual network formation and evolution require spatial initial conditions and simulation. For small X changes at fixed thermal matching, CW and m_eta^2 scale with X. X/X_ref=.8,.85,.9 with factor10 and tuned bias give max extra fractions .01695,.01748,.01797 and finite threshold ratios .9419,.8866,.8375. These are homogeneous screening rows only, not selected replacements; main candidate has more loop margin. The larger prospective factor15 scan was interrupted; no result claimed.

## Observational diagnostic: DESI DR2 and growth

Pinned public [DESI DR2 covariance and means](https://github.com/CobayaSampler/bao_data/tree/bb0c1c9009dc76d1391300e169e8df38fd1096db/desi_bao_dr2) and [published analysis](https://doi.org/10.1103/tr6y-kpc6), all 13 measurements. Profile a common BAO ruler, and separately profile Omega_m while holding H0_ref=67.4, Or, physical trigger parameters and bias fixed. No retuning activation during the profile. With Omega_m=.315, chi2 LCDM/unbiased/biased =14.30805/16.29699/16.28280 and rulers 148.63535/147.68338/147.72116 Mpc. The prior reference is reproduced to ~2.4e-10 in chi2; prior unbiased result to ~9.4e-5. With Omega_m and ruler free, chi2 LCDM/biased=10.28401/10.26311, Omega_m=.2971366/.2931377, rulers 150.66630/150.22466. These are limited BAO-only fits, not evidence in favor of the new model. Fixed-redshift drag estimates are not a calibrated ruler; fixed-parameter biased chi2=26.05845 vs 31.56755 LCDM under that approximation. Some distance quadratures flagged roundoff at very tight tolerances; no precision cosmological likelihood claimed. At fixed early psi and k=.03, present linear matter amplitude biased=-.08811976, LCDM=-.08987456, massless-eta control=-.08955027: ~1.95% and ~1.60% suppression respectively. Standard species approximated; no CMB likelihood, exact recombination or observational claim.

## Reproducibility, next plan, and sources

Run `python -m pytest -q test_metric_research.py` in the extracted research directory with NumPy, SciPy, Numba and pytest. The new suite completed 12 passed in 1.68s. Main numerical source files and scientific JSON snapshots are included in the accompanying standalone reproducibility package. No original files are replaced. Earlier output revisions predate added metadata; original slow-trigger JSON `max_x` was signed for negative trajectories; the final code uses `max_abs_x`. The compact perturbation mode stores sampled states and exact event endpoints, not all dense polynomials (earlier large parallel trial exhausted memory); a sliced collision run reported a tool timeout after its complete output had been written and was subsequently checked. Raw execution logs/cache are omitted. No credentials or personal data are included.

Next: implement a full Einstein-Boltzmann treatment (recombination, eta hierarchy, and actual nonthermal A phase-space transfer); jointly fit CMB+BAO and lensing/structure; rederive in-medium and number-changing early bounds and portal-era production; test scalar initial conditions, domain-wall avoidance, high-scale matching, and nonlinear halo evolution. The biased model is a promising conditional research candidate, not established dark matter physics. Do not discard the archived unmodified benchmark.

References: [base checkpoint](https://github.com/IskanderAbdullin/Dark-Matter/blob/14999f922ab393a94099c76e24839ff45c478be2/RESEARCH-CHECKPOINT-2026-11-25.md); [eta perturbations](https://github.com/IskanderAbdullin/Dark-Matter/blob/14999f922ab393a94099c76e24839ff45c478be2/kinetics/DIRECT-POLE-ETA-PERTURBATIONS-2026-11-25.md); [direct-pole kinetic benchmark](https://github.com/IskanderAbdullin/Dark-Matter/blob/14999f922ab393a94099c76e24839ff45c478be2/kinetics/DIRECT-POLE-KINETICS-2026-11-25.md); [matching](https://github.com/IskanderAbdullin/Dark-Matter/blob/14999f922ab393a94099c76e24839ff45c478be2/kinetics/DIRECT-POLE-MATCHING-2026-11-25.md); [initial conditions and geometry](https://github.com/IskanderAbdullin/Dark-Matter/blob/14999f922ab393a94099c76e24839ff45c478be2/kinetics/DIRECT-POLE-INITIAL-AND-GEOMETRY-2026-11-25.md); [Ma & Bertschinger](https://arxiv.org/abs/astro-ph/9506072); [Ichiki & Keum](https://doi.org/10.1088/1475-7516/2008/06/005); [Martin two-loop effective potential](https://arxiv.org/abs/hep-ph/0111209). Public production PSD inputs: [240](https://github.com/IskanderAbdullin/Dark-Matter/blob/14999f922ab393a94099c76e24839ff45c478be2/kinetics/results/direct-pole-evolution-240.json), [320](https://github.com/IskanderAbdullin/Dark-Matter/blob/14999f922ab393a94099c76e24839ff45c478be2/kinetics/results/direct-pole-evolution-320.json).

---

## Complete scientific reproducibility appendix

All embedded files are plain text, with SHA-256 checksums below. Extract each block into a fresh local directory; the prior repository files remain unchanged. The JSON snapshots contain numerical model outputs only, not execution logs.

```text
b318156eea12ef118ea6b3fadd0d3aeb184be852e85bdbd4112b48eb0f9f1291  REPORT.md
84f3e6f42cf3e4672a3611c5c51eeaecabe3c5acec00e1d30640e2e1e79fb058  acoustic.json
6e32d1f82567d8fbf53066afbfa1cb70adcc352e45720b1cc8cd23600e7f4ccb  band-.003.json
efe99f93c6ed17d69a1d479ef40b51b1d743c5f47b6c458136a7b39e66b5afae  band-.006.json
d0c918bce834e1818cfe98c9a8927d8503f226399903bc4e0bcfd6bf15dbecbb  band-.01.json
fa7bd5219aaa506e8b2a8e0b1263e706d053f27736d718d3192dcaa33eb5976c  band-.015.json
32c09c7e3b6212f5ae1c282083df58e088cd1e9754351a0f90df1f0b4e6db2db  band-.02.json
925e39bf43ed4a428fd8715dba86fd03fde91b059c5f823b16a06e1f90202aeb  band-.03.json
5a2c939e05423409cdefc726d4c8995a2a35ee183493874f1a4e64cc16cdf89f  band-.04.json
f220df17f4b8de18f26aac14b1d426070a2e29610ce752b8facd0bdb4cb6783b  band-.05.json
01a09c9a6783bd464e3cd3b072cfccc80d3795ee61f5ac96bb4df123e3981e09  band-.06.json
aa88c76640a65301aa4feea489e6f9ff5dd74ca729b601c1587977b4cfb5ba5f  bias-1.json
cd409c41c987bcbbf383b5928ff2e8b2ab1191df4070b022af338093134f8f2e  bias-2.json
d33eb7280b706e5280235b491334d958c1e19d7c0583f726c7914cd9e514e758  bias-3.json
bd542ce6c7877c6765c66a8e479d91fb0f8e11bed71f036b0342f565e1f66f86  bias-bao.json
ef51023745f286c3f7a376d946f8fccf7c033d3b01a9bf5a0bc3d2c3c3567c50  bias-early-bound.json
95fab31dd5242a1eb50b8a3cef8ef5f1ab256481dc397b62e51baac67dcef226  bias-late32.json
ad3763c794d7853920eb6953b67f6859daf9a1aec0f5679dbd0e0da11e6bc202  bias-late35.json
476d4ab675269118e52ae95a5db0b082891738556750642b670ff07715bbe5c4  bias-late38.json
11735f8bfb8028d8e3ec82304f224d8052414aff9c6c0c9d673c9cad148f3630  bias-local-modes.json
16eb1f22bb000896c476c365fdc5de8674e6042b92d600a7662582fea2039fb5  bias-seed-wall.json
997bbee3d2e92d6fdc3aa7264c243b945513581e50f8424bb5f8b76e14febcc4  bias-window-0.json
38fb2a4df58c5c36e810f630d6e17400efc44ca0bc2e363b1f54735eca497bdb  bias10-acoustic-fine-k.03.json
b25ddab00aaa1f51140e31088c66cd7a8aeb3ae2f21c3277e5c7fb11329ee3c6  bias10-acoustic-fine-k.05.json
9a84304a4a2739fe4e805b7cfc171da898e8f6b49e16e0d39f2d1455152a0f32  bias10-acoustic-fine-k.06.json
5752431918cf091dae0c0b08280d50e70a21ec22dd944e11b51bc333882bd748  bias10-acoustic-k03.json
b6b34f393ef0b0d0539d674b1f49e44eee6fe7b05f14c2038cd73d29440f72fa  bias10-chi_iso-k003.json
2e6a6dd0dc6f0c1f982526ac801744d52bbbb3e98d4876c39b75e26cf91c5a59  bias10-early-k.0003.json
ba441e69e43646237e61b44ea99582e33016e954c1a2461340c5071cb9bee1ca  bias10-early-k.001.json
67c576999a56c369ee9d2c65f10b3867c4569cf5583437e087a5879063c0bac1  bias10-early-k.003.json
0f85fb13894582826f8d043b077c5db9c4cc90e1fadb1a8e6b31b6ba445bc42d  bias10-early-k.006.json
efa1f9681506b3e6f7f5abca49fa489c353950f2af567f8a9869482dfd0e207e  bias10-early-k.01.json
15a9e521e695c0351a88ad6f5259f1ae5dba51b400ed0f8ccfb20d6929877b67  bias10-early-k.02.json
0bd42a31c57b7a7161258cddbe1ebc752d88eefa31c37a2167610ad6e142b714  bias10-eta_iso-k003.json
7c7a380550539a9502d8253fff78c5c6a4f9c48049ddb6a1f1f5ee139086becf  bias10-k.001.json
7ba9fe1bdf58c593af4c688dcc9340c2a7eb8d3c2dcb8861cdb7c145f3a62ac3  bias10-k.002.json
ae9883c17bd51dfbb9fd819901441151f338f660a198f802b4cfb9b4d2ee6a09  bias10-k.003.json
d000c54943ad263cd5f982e5d69fc2ea0c35b7cedbb7717bd7b2f645466ac6be  bias10-k.006.json
07c53f344295819ee5e563434de5c1f249231ab0fcdf6f221f846ecc8f973d81  bias10-k.01.json
498ea19d3c254acc6643ccf5e810095820ad713c73ba0a4905919b56513e681c  control-no-late-metric.json
9add24639604f45109e9922e25d964279a586d8d9be8d815bda5afb96b69aa1d  direct_pole_bias_background_checks.py
16f7f54775a39241b7f17f95c14fa1b14e3b2e4d206018527c0c1c306d3ce462  direct_pole_bias_bao.py
e40c0e7968919404cdf83ff6bc4e531f563df3b64f3bbb495182f4bfb5030ebb  direct_pole_bias_early_bound.py
0d73f6a9db0e6f5e63bed0cff2b471f2677900613e82386bba008dd19e99a9e6  direct_pole_bias_elastic.py
361f78c8dfde4b321f89ad7e02b02d122697c835c1bac79e2e5a0a11a604b2ef  direct_pole_bias_local_modes.py
2beea9f18fe384333ede52a038ffc318b2a14b642a2ee28cd3a30b65de3f66ae  direct_pole_bias_parameter_window.py
34c74a2bdebee8b34fbb4268f470948ea63a30b3aa01ee3aaee31dbbc3ddb65b  direct_pole_metric_acoustic.py
6396d68172dbbfdca0748620ee338f30ece28d22170088f0079f44b61439d7b7  direct_pole_metric_pilot.py
f4d73c83de8ca1f62f8738834400ce569a0b3b060e114d1918f0edee100b9773  direct_pole_slow_trigger.py
1ebfe5662c7ad363d7e131d134a5300577b98ad6a992f459cc23e7ad8c209128  direct_pole_trigger_shape.py
4489d28caf2ef81a29f37ce06e3f547558bf5cae56571a62f86e16d79ebfdc67  elastic-base240.json
895588bb3e1490ceb78b9dcccc2ca80c8abe377bf8799b8699cf87b0173f8a98  elastic-bias240-sliced.json
2225e3b1345d8a5fb42122ee5e8051b00c7a781d3673a1290b246c1c86e30ce4  elastic-bias240.json
ccbc5e0ec6b111fb81f5f54556a62fdfbd9dbaf7e970b9ed183318fe0367e9c7  elastic-bias320.json
e1ae55b9ea1457342f63b7bb244d9a776d8581f02183326bf3a652098e8f330e  evolved-phi-fine.json
e72fd4a779a03bd4b4d6c13bbb95df656f288fe101d516f4e994ccf0780c1fcd  evolved-phi.json
70271fcfd81ec107f3d32c72bca7b8c84f53281b0df2a0a61c1f697875620265  free-rad-early.json
775e158f5ec26e2d34a92d0f510575afa4231aa0df12f3904eee62fc294d7ee1  free-rad-fine.json
88b764abb0d9d7af0d5f50b859c2e7b4af58ea686433d0c29e9586b32b9aa9db  free-rad.json
b1cf4bf13dec4ecbde12df2590a93a42e9c6b40544c204f0e39dad117d0102f0  growth-k03-full.json
bd4649687c785e0b3116106cf9c61e77ea2677d1cd01cd8044d5f71042d54f07  growth-k03-lcdm.json
dee6ba9b88489346d04cfdb1a0db4e06df0fd9731f813bb23e478196f2917932  growth-k03-massless_eta.json
a36788bc0404b0141c1533f89528a4df7bf0c5277f2092f4244dfeaf75407761  k02.json
871440446b99274aa3409dd72bb452d5a38dc0057821379c59852764c4080fcf  k04.json
4cb1c0a61e74834e1053ea36bf1630891d9374b266e6a105d238a0c836f77660  pilot-1.json
b8bdb169fcd5a4d09a09388eb07e30dd13b95d4c6069ec7417203daebf347e46  pilot-2.json
53ff3ce4bf70057919d99db138ad4a416f64e5601f1cf70bf8ae65a0691d5583  psd_reference.py
07b23feb9390d907c4a0144512a430d6480b54a5282bfc08407b83353d238455  reproduce.sh
a1f530e13f46bbed2c482046c8fb5ea82f21202b991cdb29d3dd64ddd0cf286f  requirements.txt
11fe4ea2eeea2479a5a6cff9e6fe86f601309e79aef3f6474135039408bfd364  seed-double.json
d82edbff47cc411d50af771979ecb9793b0e23c7dfdf2cfddad7745fe05dd527  seed-half.json
190a8b5b85afc4375280c6495f163db4af35ad0cab7690e168626d20f242877e  shape--0.1-k03.json
d37bcb4132b847f31bc4004f134ee86286134931146a13e3a04d42e9d91c2559  shape--0.14-k03.json
edb0e77020e184eaf0ff4b3c97db0ea806bc453e88ae1bd7086fd63582fd7253  shape-0.1-k03.json
a7eb998414a3e4ef38eadd2a9b152d6f8b21699a810e57f63da765e645827ba2  shape-0.2-k03.json
13faadc716b937d27f927d0fa754721d700e114531340f5ee93e87e42a24c4ff  shape-screen.json
f6bfcc3c3d2b5dcc269ff5f8b69be170768611e42739d03f829d0b2c7ddbba32  slow-1.json
3ac18e9dd04d3978c08348f20c8efb2ea6edebe66bfa59f81775fcf2c52349bb  slow-2.json
24a30841d56306248402b771beafadd68820a6c73c926cf6b97e395496e382b1  slow-3.json
4d23d72d17bc2dbe334303dcbe31a7e320efce8ceab1791fb2325d555dd25881  slow-4.json
ce6a544a5ab004c2e0364aed3f17e605e92e3a970d4054db44a96d15545c3b9b  slow-5.json
bd71eeee58ec15c8c39d806bb52237c0fcde1a6b396477c22516ebf70c45fe7d  slow2-k015.json
002f294147e12da7c5c144b25904cde44a6811f0360629718914d8cc2a2d3ba3  slow25-k012.json
211ef6ea0a9317e8d9dd4dcb6ea7056fd7cb12d96b572eaadad2bade8856c151  slow275-k011.json
16e95c2a81672fc7bbf41169ffcd6c8aefea536e89846222808521ff728225fb  slow3-k01.json
7883fe1e8b18671744dcf0a5d45b3eac12b962bc97cf3104ee780cc29aa06af2  small-amplitude-fine.json
48e0442a7b7f3e4bcea5974205a7cd78130926001b51b98d05d9332b8524925c  small-amplitude.json
9726ed97d9f3cbda32fa2d1d0b15c03c9104acc312a92201eeff53a97e5019bd  test_metric_research.py
6542946e0dc832e92d51e2d7a873653b6a128c568ddcc5368f35e56995d56768  threshold-point1-new.json
7e267c714e6f7eacbd29165b09a37f1672b1b91931e89f69c5ad5b34261a65c4  threshold-point1.json
3550536c691ca76488011e608019e81674c8530f029564227d991f28da8e25ef  tiny-amplitude-fine.json
d16ef0adcf9dde8f81300ac3e58ba027d615e90e16a9f135ff58e4b90ed265ff  tiny-amplitude.json
```

To extract the blocks from a local copy of this Markdown document (Python 3, in the desired fresh research directory):

```python
from pathlib import Path
import re
p = Path("DIRECT-POLE-METRIC-BIASED-TRIGGER-2026-11-25.md")
for name, body in re.findall(r"^<!-- FILE: ([A-Za-z0-9_.-]+) -->\n```[^\n]*\n(.*?)^```\n<!-- END FILE -->", p.read_text(), re.M|re.S):
    Path(name).write_text(body)
```

### acoustic.json

<!-- FILE: acoustic.json -->
```json
{
  "k": 0.03,
  "order": 5,
  "angles": 96,
  "rtol": 1e-07,
  "zdec": 1090.0,
  "zmatch": 12.0,
  "baryon_fraction": 0.15555555555555556,
  "nu_fraction": 0.409,
  "amplitude": 1e-05,
  "z_at_stop": 3.3710650055612494,
  "reached_eta_threshold": true,
  "reached_chi_threshold": false,
  "max_eta": 0.10000000000000023,
  "max_chi_rel": 0.053317906841659,
  "eta_rest_at_stop": -0.0988114845501542,
  "max_matter": 0.025772276547787517,
  "max_density_constraint_rel": 1.2900548674324235e-08,
  "match_eta_number_per_psi": -0.9228254582370714,
  "match_matter_per_psi": -876.4707550929982,
  "scope": "Sudden decoupling, ideal tight coupling before z_dec, collisionless photons thereafter; no visibility function, polarization or full recombination. This is a sensitivity diagnostic."
}
```
<!-- END FILE -->

### band-.003.json

<!-- FILE: band-.003.json -->
```json
{
  "k": 0.003,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000001480981622,
  "match_delta_m_per_psi": -28.50163219086404,
  "match_eta_delta_n_per_psi": -2.7514352116858034,
  "z_at_stop": 2.1256213476622774,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": true,
  "max_eta_number_fraction": 0.034658000582110875,
  "max_matter_fraction": 0.0011103734230368148,
  "max_phi": 1.641248758650519e-05,
  "max_delta_x": 0.0982075756627036,
  "max_delta_x_over_x_floor": 0.10130809381929246,
  "at_stop": {
    "k_over_aH": 12.868154488558813,
    "eta_number": -0.034658000582110875,
    "eta_number_rest_frame": -0.0323784649982151,
    "eta_background_fraction": 0.00014104285757827434,
    "eta_density_over_total": 9.079901104491782e-06,
    "delta_x_over_x": 0.09999999999999896
  },
  "max_constraint_derivative_rel": 9.172904201389575e-06,
  "rms_constraint_derivative_rel": 1.1500381479765627e-06,
  "max_density_constraint_rel": 9.399488755815045e-06,
  "early_nfev": 23771,
  "late_nfev": 10649,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### band-.006.json

<!-- FILE: band-.006.json -->
```json
{
  "k": 0.006,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.0000059239264862,
  "match_delta_m_per_psi": -96.30005700481057,
  "match_eta_delta_n_per_psi": -2.4556046560449896,
  "z_at_stop": 1.1906976452466325,
  "reached_eta_threshold": true,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.1,
  "max_matter_fraction": 0.005456201412441699,
  "max_phi": 1.2320158870722718e-05,
  "max_delta_x": 0.0694290352851755,
  "max_delta_x_over_x_floor": 0.07053023299931377,
  "at_stop": {
    "k_over_aH": 28.993668565101405,
    "eta_number": -0.1,
    "eta_number_rest_frame": -0.09850729555803513,
    "eta_background_fraction": 0.0001323303222854483,
    "eta_density_over_total": -2.1206697465437323e-05,
    "delta_x_over_x": -0.06040578227933636
  },
  "max_constraint_derivative_rel": 1.1135152514885423e-05,
  "rms_constraint_derivative_rel": 1.0440272290101446e-06,
  "max_density_constraint_rel": 4.184060742104312e-06,
  "early_nfev": 15236,
  "late_nfev": 18209,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### band-.01.json

<!-- FILE: band-.01.json -->
```json
{
  "k": 0.01,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.00001645535135,
  "match_delta_m_per_psi": -225.9309068529801,
  "match_eta_delta_n_per_psi": -2.1207389608963827,
  "z_at_stop": 2.604190164243615,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": true,
  "max_eta_number_fraction": 0.08659851664454692,
  "max_matter_fraction": 0.007984838286436905,
  "max_phi": 8.780553706491333e-06,
  "max_delta_x": 0.08810795189485644,
  "max_delta_x_over_x_floor": 0.10000000000000281,
  "at_stop": {
    "k_over_aH": 40.40873611259203,
    "eta_number": -0.08659851664454692,
    "eta_number_rest_frame": -0.08439507295052949,
    "eta_background_fraction": 0.00012529341923424723,
    "eta_density_over_total": 1.4427219530673517e-06,
    "delta_x_over_x": 0.10000000000000281
  },
  "max_constraint_derivative_rel": 5.368424759140705e-06,
  "rms_constraint_derivative_rel": 5.559536742044799e-07,
  "max_density_constraint_rel": 1.094327330628226e-06,
  "early_nfev": 11621,
  "late_nfev": 8663,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### band-.015.json

<!-- FILE: band-.015.json -->
```json
{
  "k": 0.015,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000037024540537,
  "match_delta_m_per_psi": -421.9577432873175,
  "match_eta_delta_n_per_psi": -1.7615186712637083,
  "z_at_stop": 3.0992874049248846,
  "reached_eta_threshold": true,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.10000000000000006,
  "max_matter_fraction": 0.013189874348015503,
  "max_phi": 6.804824754767071e-06,
  "max_delta_x": 0.08521610582199565,
  "max_delta_x_over_x_floor": 0.09548537416232805,
  "at_stop": {
    "k_over_aH": 57.251076341223346,
    "eta_number": -0.10000000000000006,
    "eta_number_rest_frame": -0.09808142109454332,
    "eta_background_fraction": 0.0001449136151640442,
    "eta_density_over_total": -2.4759315824060854e-05,
    "delta_x_over_x": -0.07160870360325228
  },
  "max_constraint_derivative_rel": 5.446411755899075e-06,
  "rms_constraint_derivative_rel": 3.544780129534272e-07,
  "max_density_constraint_rel": 5.278879900552543e-07,
  "early_nfev": 9791,
  "late_nfev": 6998,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### band-.02.json

<!-- FILE: band-.02.json -->
```json
{
  "k": 0.02,
  "radial_order": 5,
  "angles": 64,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.0000658214053995,
  "match_delta_m_per_psi": -634.9334516627928,
  "match_eta_delta_n_per_psi": -1.4554274788053378,
  "z_at_stop": 3.3025060747277535,
  "reached_eta_threshold": true,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.09999999999999983,
  "max_matter_fraction": 0.01894900396481864,
  "max_phi": 5.438868385610656e-06,
  "max_delta_x": 0.05848664133644327,
  "max_delta_x_over_x_floor": 0.06814170864282199,
  "at_stop": {
    "k_over_aH": 74.66932221912955,
    "eta_number": -0.09999999999999983,
    "eta_number_rest_frame": -0.09860085322077695,
    "eta_background_fraction": 0.0001454770191683912,
    "eta_density_over_total": -2.1046219887068692e-05,
    "delta_x_over_x": -0.04505906856936527
  },
  "max_constraint_derivative_rel": 7.164023253915697e-06,
  "rms_constraint_derivative_rel": 6.078532122380861e-07,
  "max_density_constraint_rel": 2.712464819033737e-07,
  "early_nfev": 8918,
  "late_nfev": 6443,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### band-.03.json

<!-- FILE: band-.03.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4088035686345,
  "match_eta_delta_n_per_psi": -1.123519554436656,
  "z_at_stop": 3.390754825534185,
  "reached_eta_threshold": true,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.0999999999999999,
  "max_matter_fraction": 0.03133711969214545,
  "max_phi": 3.81645902950761e-06,
  "max_delta_x": 0.020241700435602535,
  "max_delta_x_over_x_floor": 0.0319516155407187,
  "at_stop": {
    "k_over_aH": 110.96188931506147,
    "eta_number": -0.0999999999999999,
    "eta_number_rest_frame": -0.09905977652824009,
    "eta_background_fraction": 9.636805137606811e-05,
    "eta_density_over_total": -6.75207191384444e-06,
    "delta_x_over_x": 0.0319516155407187
  },
  "max_constraint_derivative_rel": 3.3709966818684456e-06,
  "rms_constraint_derivative_rel": 2.428638382288629e-07,
  "max_density_constraint_rel": 5.1942830824648834e-08,
  "early_nfev": 8234,
  "late_nfev": 6155,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### band-.04.json

<!-- FILE: band-.04.json -->
```json
{
  "k": 0.04,
  "radial_order": 5,
  "angles": 96,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.0002632856215983,
  "match_delta_m_per_psi": -1491.0348602663244,
  "match_eta_delta_n_per_psi": -0.8563670555362346,
  "z_at_stop": 3.1713720143100494,
  "reached_eta_threshold": true,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.1000000000000003,
  "max_matter_fraction": 0.045895146676407456,
  "max_phi": 2.906900908651542e-06,
  "max_delta_x": 0.010026352624137962,
  "max_delta_x_over_x_floor": 0.016347829048476897,
  "at_stop": {
    "k_over_aH": 151.45967034715196,
    "eta_number": -0.1000000000000003,
    "eta_number_rest_frame": -0.09953102149856655,
    "eta_background_fraction": 0.00010084605022982947,
    "eta_density_over_total": -8.525479381426938e-06,
    "delta_x_over_x": 0.016347829048476897
  },
  "max_constraint_derivative_rel": 1.745173760959804e-05,
  "rms_constraint_derivative_rel": 1.0658459959319337e-06,
  "max_density_constraint_rel": 2.999051901714456e-08,
  "early_nfev": 8204,
  "late_nfev": 6713,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### band-.05.json

<!-- FILE: band-.05.json -->
```json
{
  "k": 0.05,
  "radial_order": 5,
  "angles": 128,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000411383783747,
  "match_delta_m_per_psi": -1886.0533676333198,
  "match_eta_delta_n_per_psi": -0.7090928112157937,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.03778992857540149,
  "max_matter_fraction": 0.18964891317984528,
  "max_phi": 2.3416828411348733e-06,
  "max_delta_x": 0.00011554115016058865,
  "max_delta_x_over_x_floor": 0.0001488987153455544,
  "at_stop": {
    "k_over_aH": 221.68048512465805,
    "eta_number": -0.03778992857540149,
    "eta_number_rest_frame": -0.03778965257301596,
    "eta_background_fraction": 5.193108415566764e-05,
    "eta_density_over_total": -1.9614503037887207e-06,
    "delta_x_over_x": 8.586107561498444e-06
  },
  "max_constraint_derivative_rel": 5.755776622080978e-08,
  "rms_constraint_derivative_rel": 7.172486633876381e-09,
  "max_density_constraint_rel": 6.469103654323335e-10,
  "early_nfev": 8459,
  "late_nfev": 47885,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### band-.06.json

<!-- FILE: band-.06.json -->
```json
{
  "k": 0.06,
  "radial_order": 5,
  "angles": 144,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000592392648596,
  "match_delta_m_per_psi": -2253.973619123714,
  "match_eta_delta_n_per_psi": -0.5810447323133531,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.01440174153007267,
  "max_matter_fraction": 0.22664000850803007,
  "max_phi": 1.9435429964942427e-06,
  "max_delta_x": 1.0064660006315838e-05,
  "max_delta_x_over_x_floor": 1.3123952125933007e-05,
  "at_stop": {
    "k_over_aH": 266.01658214958974,
    "eta_number": -0.01440174153007267,
    "eta_number_rest_frame": -0.014401268106317677,
    "eta_background_fraction": 5.193108415566764e-05,
    "eta_density_over_total": -7.475239951248087e-07,
    "delta_x_over_x": 3.831780673610502e-06
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 2.5639301247586045e-10,
  "early_nfev": 8834,
  "late_nfev": 42329,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias-1.json

<!-- FILE: bias-1.json -->
```json
[
  {
    "vev_factor": 2.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 1.2888956122800165e-08,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.01754745221126468,
    "positive_local_root_integral": 5.440288704994039,
    "minimum_curvature_after_half": -4097.399666977045,
    "late_contact_proxy": 0.009175370615111393,
    "early_contact_proxy_100_to_1e5": 1.172777357476351e-24,
    "early_representative_pair_proxy": 1.163169420559652e-24,
    "early_x_z1e5": 9.999593183592548e-09,
    "max_abs_x": 1.2949179371999262,
    "today_x": 0.9669979477266565,
    "activation_integral": 4.96074981033187e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 3.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 2.6780974881820943e-06,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.01675713797917659,
    "positive_local_root_integral": 2.9444928039056775,
    "minimum_curvature_after_half": -1731.1820722762955,
    "late_contact_proxy": 0.00954913354203596,
    "early_contact_proxy_100_to_1e5": 1.1752947357454208e-24,
    "early_representative_pair_proxy": 1.1656848763840047e-24,
    "early_x_z1e5": 9.999819195208545e-09,
    "max_abs_x": 1.287975231829415,
    "today_x": 1.0069364384181283,
    "activation_integral": 4.960633056726472e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  }
]
```
<!-- END FILE -->

### bias-2.json

<!-- FILE: bias-2.json -->
```json
[
  {
    "vev_factor": 4.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 3.931350584814036e-05,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.016068714586544824,
    "positive_local_root_integral": 1.7570284230317448,
    "minimum_curvature_after_half": -914.4040428887374,
    "late_contact_proxy": 0.009842965827726988,
    "early_contact_proxy_100_to_1e5": 1.0359466543968231e-23,
    "early_representative_pair_proxy": 1.0349892763313231e-23,
    "early_x_z1e5": 9.99989832589195e-09,
    "max_abs_x": 1.2810990254609314,
    "today_x": 1.0241694287383043,
    "activation_integral": 4.979737041103193e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 6.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 0.0006047914848040583,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.014906171617716306,
    "positive_local_root_integral": 0.713641049605639,
    "minimum_curvature_after_half": -311.4192534926019,
    "late_contact_proxy": 0.010413843944601622,
    "early_contact_proxy_100_to_1e5": 2.2840092111049025e-20,
    "early_representative_pair_proxy": 2.2840082054538315e-20,
    "early_x_z1e5": 9.99995500920055e-09,
    "max_abs_x": 1.2679420489927624,
    "today_x": 1.0290003606865066,
    "activation_integral": 4.97903744183086e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  }
]
```
<!-- END FILE -->

### bias-3.json

<!-- FILE: bias-3.json -->
```json
[
  {
    "vev_factor": 8.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 0.002585274451356568,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.014039079216465127,
    "positive_local_root_integral": 0.3119074657968981,
    "minimum_curvature_after_half": -100.686102109237,
    "late_contact_proxy": 0.010911767294022934,
    "early_contact_proxy_100_to_1e5": 8.042707430732143e-19,
    "early_representative_pair_proxy": 8.042707165654377e-19,
    "early_x_z1e5": 9.999975082277226e-09,
    "max_abs_x": 1.2558819244438566,
    "today_x": 0.9709928551315407,
    "activation_integral": 4.969259543257345e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 10.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 0.006704850663097823,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.013461782233096992,
    "positive_local_root_integral": 0.018969767898436678,
    "minimum_curvature_after_half": -7.175357291342785,
    "late_contact_proxy": 0.011224906880376506,
    "early_contact_proxy_100_to_1e5": 6.199695712482126e-18,
    "early_representative_pair_proxy": 6.199695572686891e-18,
    "early_x_z1e5": 9.999984571332809e-09,
    "max_abs_x": 1.2451445964985957,
    "today_x": 0.9729327743493299,
    "activation_integral": 5.010446903061306e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  }
]
```
<!-- END FILE -->

### bias-bao.json

<!-- FILE: bias-bao.json -->
```json
{
  "source": "https://github.com/CobayaSampler/bao_data/tree/bb0c1c9009dc76d1391300e169e8df38fd1096db/desi_bao_dr2",
  "source_cov_sha256": "252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509",
  "source_mean_sha256": "9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585",
  "results": {
    "lcdm": {
      "fixed_omega_m": {
        "model": "lcdm",
        "omega_m": 0.315,
        "distance_Mpc": [
          1184.5194606096752,
          1984.8265220326516,
          3343.3704213374767,
          2602.3796799863026,
          2965.4864523582373,
          3233.7453438009443,
          2583.2051784064874,
          4128.7297720933875,
          2068.0077199099997,
          4451.313482962551,
          1893.7967155810231,
          1266.8227730977144,
          5760.928209727279
        ],
        "ruler_Mpc": 148.6353547413654,
        "ruler_sigma_local": 0.41171988521313746,
        "chi2": 14.308052790055788,
        "predicted": [
          7.969298170484481,
          13.35366357140514,
          22.493776310186398,
          17.508483661337518,
          19.95142042428845,
          21.756232556028532,
          17.37948002278073,
          27.777575390980356,
          13.91329622422915,
          29.947878085318983,
          12.741226465778245,
          8.5230245206604,
          38.75880149612887
        ],
        "approx_fixed_drag_ruler": 146.94434306452612,
        "approx_fixed_drag_chi2": 31.567545693000763,
        "D_M_z1090": 13864.3590556941,
        "r_s_z1090": 144.28738638432625,
        "theta_s_z1090": 0.010407072249406823,
        "H0_ratio": 1.0
      },
      "profiled_omega_m": {
        "model": "lcdm",
        "omega_m": 0.29713655087939334,
        "distance_Mpc": [
          1190.5888999878223,
          1997.286054882685,
          3385.366159878759,
          2623.7157711705245,
          3013.3024631866874,
          3266.298257750106,
          2633.127674646798,
          4180.292608792536,
          2115.526458088131,
          4510.462988462775,
          1939.3174424048311,
          1301.1426231397425,
          5853.758882782577
        ],
        "ruler_Mpc": 150.6663044670646,
        "ruler_sigma_local": 0.4173391724574264,
        "chi2": 10.284009526928266,
        "predicted": [
          7.9021577133597445,
          13.256355241123527,
          22.46929843971048,
          17.414084592113056,
          19.99984318886238,
          21.679022853209446,
          17.476553128190616,
          27.745371624923216,
          14.041138564930964,
          29.93677321825302,
          12.871606888245955,
          8.635923126555282,
          38.852475365931596
        ],
        "approx_fixed_drag_ruler": 149.07832625770337,
        "approx_fixed_drag_chi2": 25.072165821561136,
        "D_M_z1090": 14190.239044533242,
        "r_s_z1090": 146.36212072659995,
        "theta_s_z1090": 0.010314281547144593,
        "H0_ratio": 1.0
      }
    },
    "unbiased": {
      "fixed_omega_m": {
        "model": "unbiased",
        "omega_m": 0.315,
        "distance_Mpc": [
          1179.2279130095797,
          1975.4229245204424,
          3322.5549632565812,
          2588.8732628456946,
          2944.615870768731,
          3215.557597703105,
          2563.1549497533597,
          4103.221665975197,
          2050.293530490785,
          4423.011676363565,
          1877.2448836862472,
          1255.1747859862658,
          5720.7848782824
        ],
        "ruler_Mpc": 147.6833811685567,
        "ruler_sigma_local": 0.40908604024050305,
        "chi2": 16.296986957851953,
        "predicted": [
          7.984838264663521,
          13.376067834374787,
          22.497825665735686,
          17.5298888904156,
          19.93870838728921,
          21.773320547374695,
          17.355743953532137,
          27.783909289644665,
          13.8830348700556,
          29.949285026968692,
          12.711280503143923,
          8.499092965332956,
          38.73682220041434
        ],
        "approx_fixed_drag_ruler": 146.45277228861642,
        "approx_fixed_drag_chi2": 25.498914921281393,
        "D_M_z1090": 13798.975995413157,
        "r_s_z1090": 143.80007493093626,
        "theta_s_z1090": 0.01042106856180749,
        "H0_ratio": 1.0032365535031618
      }
    },
    "biased": {
      "fixed_omega_m": {
        "model": "biased",
        "omega_m": 0.315,
        "distance_Mpc": [
          1179.4541858296275,
          1975.808582284189,
          3323.5428637126993,
          2589.428469487214,
          2945.656334731121,
          3216.300343510529,
          2564.165872942781,
          4104.260619983086,
          2050.628547057564,
          4424.166623100028,
          1878.1394399472572,
          1255.4921351068772,
          5722.364807922196
        ],
        "ruler_Mpc": 147.72116413428972,
        "ruler_sigma_local": 0.409190677574323,
        "chi2": 16.282802403590136,
        "predicted": [
          7.984327721364382,
          13.375257322559612,
          22.49875895028384,
          17.52916370963085,
          19.94065205208711,
          21.772779563168537,
          17.35814829222277,
          27.783836148569762,
          13.881751874046884,
          29.949443257014448,
          12.714085019258892,
          8.499067431972982,
          38.7376097491361
        ],
        "approx_fixed_drag_ruler": 146.45277228815053,
        "approx_fixed_drag_chi2": 26.058452992850683,
        "D_M_z1090": 13800.711614134258,
        "r_s_z1090": 143.80007493051696,
        "theta_s_z1090": 0.010419757977062677,
        "H0_ratio": 1.0030711878986183
      },
      "profiled_omega_m": {
        "model": "biased",
        "omega_m": 0.2931376600285811,
        "distance_Mpc": [
          1187.1371455348362,
          1991.4574032224423,
          3375.6783836742047,
          2616.0450635474385,
          3004.3520767134873,
          3256.734317349317,
          2625.5347196758084,
          4168.042219943653,
          2109.6303723871224,
          4497.255460029974,
          1933.4945388750975,
          1297.056115428739,
          5836.723539436663
        ],
        "ruler_Mpc": 150.22466136158005,
        "ruler_sigma_local": 0.41611580671484544,
        "chi2": 10.263107573909608,
        "predicted": [
          7.902411859511415,
          13.256527824210876,
          22.470866987339633,
          17.414218410190355,
          19.99906040381896,
          21.679092419523514,
          17.47738817234764,
          27.745392681641487,
          14.043169432143985,
          29.93686535398739,
          12.87068661929824,
          8.634109098151452,
          38.85329803066146
        ],
        "approx_fixed_drag_ruler": 149.05768533420255,
        "approx_fixed_drag_chi2": 18.251679762899048,
        "D_M_z1090": 14201.426605203584,
        "r_s_z1090": 146.332284619324,
        "theta_s_z1090": 0.010304055267638113,
        "H0_ratio": 1.0028723321657982
      }
    }
  },
  "scope": "DESI DR2 only, fixed H0=67.4, fixed Or, physical trigger parameters held fixed during Omega_m profile; fixed-redshift sound horizon diagnostic not a fitted ruler."
}
```
<!-- END FILE -->

### bias-early-bound.json

<!-- FILE: bias-early-bound.json -->
```json
{
  "matching": {
    "g2_eV2": 5.1650631107508515e-11,
    "kappa": 7.941424108611536e-32,
    "lambda_S": 3.4170831596658935e-06,
    "field_cap_over_Mpl": 0.1,
    "lambda_A": 1.683149091809942e-21,
    "lambda_threshold_bound": 1.6136637954716584e-08,
    "original_finite_loop_ratio": 0.753975227799072,
    "constant_A4_finite_loop_ratio": 0.6736618798910903,
    "scope": "Finite-field tree bound and finite scale-local thresholds only. Global quartic witness unchanged under v rescaling up to the slight coupling retune."
  },
  "rows": [
    {
      "seed": 0.0,
      "mesh": 240,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012999474738298878,
      "constant_only_rge101_L1_bound": 0.00012999477012871268,
      "combined_amplitude_rge101_L1_bound": 0.0001329320868385461,
      "positive_mode_gain_to_prior_production_start": 7672293.84503358,
      "occupation_cap": 21525656357.32168,
      "rearrangement_transition_r": 1511583871.2165565,
      "lambda_at_z100": 4.622107727166666e-18,
      "peak_abs_x_history": 1.6921326880609287e-05,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 1e-08,
      "mesh": 240,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012999474738298878,
      "constant_only_rge101_L1_bound": 0.00012999477012871268,
      "combined_amplitude_rge101_L1_bound": 0.00013319018237583486,
      "positive_mode_gain_to_prior_production_start": 7687012.01535689,
      "occupation_cap": 21525656357.32168,
      "rearrangement_transition_r": 1511583871.2165565,
      "lambda_at_z100": 4.627035149698633e-18,
      "peak_abs_x_history": 1.6930347297256456e-05,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 1e-07,
      "mesh": 240,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012999474738298878,
      "constant_only_rge101_L1_bound": 0.00012999477012871268,
      "combined_amplitude_rge101_L1_bound": 0.00015914079037237968,
      "positive_mode_gain_to_prior_production_start": 9213924.690200044,
      "occupation_cap": 21525656357.32168,
      "rearrangement_transition_r": 1511583871.2165565,
      "lambda_at_z100": 4.6715001230095655e-18,
      "peak_abs_x_history": 1.7011531047278363e-05,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 1e-06,
      "mesh": 240,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012999474738298878,
      "constant_only_rge101_L1_bound": 0.00012999477012871268,
      "combined_amplitude_rge101_L1_bound": 0.014574810224478998,
      "positive_mode_gain_to_prior_production_start": 859973508.527087,
      "occupation_cap": 21525656357.32168,
      "rearrangement_transition_r": 1511583871.2165565,
      "lambda_at_z100": 5.127848739287627e-18,
      "peak_abs_x_history": 1.7823368568430626e-05,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 3e-06,
      "mesh": 240,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012999474738298878,
      "constant_only_rge101_L1_bound": 0.00012999477012871268,
      "combined_amplitude_rge101_L1_bound": 0.9903837546175369,
      "positive_mode_gain_to_prior_production_start": 58452065419.96164,
      "occupation_cap": 21525656357.32168,
      "rearrangement_transition_r": 1511583871.2165565,
      "lambda_at_z100": 6.218111138985368e-18,
      "peak_abs_x_history": 1.9627452093810485e-05,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 1e-05,
      "mesh": 240,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012999474738298878,
      "constant_only_rge101_L1_bound": 0.00012999477012871268,
      "combined_amplitude_rge101_L1_bound": 119.73260784591447,
      "positive_mode_gain_to_prior_production_start": 7066619095434.186,
      "occupation_cap": 21525656357.32168,
      "rearrangement_transition_r": 1511583871.2165565,
      "lambda_at_z100": 1.0861223917530709e-17,
      "peak_abs_x_history": 2.5941745980524912e-05,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 0.0001,
      "mesh": 240,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012999474738298878,
      "constant_only_rge101_L1_bound": 0.00012999477012871268,
      "combined_amplitude_rge101_L1_bound": 1194856.8148779587,
      "positive_mode_gain_to_prior_production_start": 7.05204741667078e+16,
      "occupation_cap": 21525656357.32168,
      "rearrangement_transition_r": 1511583871.2165565,
      "lambda_at_z100": 1.8518450713925466e-16,
      "peak_abs_x_history": 0.00010712573723895841,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 0.0,
      "mesh": 320,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012931924717677148,
      "constant_only_rge101_L1_bound": 0.0001293192146424482,
      "combined_amplitude_rge101_L1_bound": 0.00013224507637227422,
      "positive_mode_gain_to_prior_production_start": 7603108.977363672,
      "occupation_cap": 21441710365.93471,
      "rearrangement_transition_r": 1509616347.5374174,
      "lambda_at_z100": 4.622107727166666e-18,
      "peak_abs_x_history": 1.6921326880609287e-05,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 1e-08,
      "mesh": 320,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012931924717677148,
      "constant_only_rge101_L1_bound": 0.0001293192146424482,
      "combined_amplitude_rge101_L1_bound": 0.00013250184195051225,
      "positive_mode_gain_to_prior_production_start": 7617694.426666811,
      "occupation_cap": 21441710365.93471,
      "rearrangement_transition_r": 1509616347.5374174,
      "lambda_at_z100": 4.627035149698633e-18,
      "peak_abs_x_history": 1.6930347297256456e-05,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 1e-07,
      "mesh": 320,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012931924717677148,
      "constant_only_rge101_L1_bound": 0.0001293192146424482,
      "combined_amplitude_rge101_L1_bound": 0.00015831769321736266,
      "positive_mode_gain_to_prior_production_start": 9130838.174838435,
      "occupation_cap": 21441710365.93471,
      "rearrangement_transition_r": 1509616347.5374174,
      "lambda_at_z100": 4.6715001230095655e-18,
      "peak_abs_x_history": 1.7011531047278363e-05,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 1e-06,
      "mesh": 320,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012931924717677148,
      "constant_only_rge101_L1_bound": 0.0001293192146424482,
      "combined_amplitude_rge101_L1_bound": 0.014499073100245646,
      "positive_mode_gain_to_prior_production_start": 852218702.1301724,
      "occupation_cap": 21441710365.93471,
      "rearrangement_transition_r": 1509616347.5374174,
      "lambda_at_z100": 5.127848739287627e-18,
      "peak_abs_x_history": 1.7823368568430626e-05,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 3e-06,
      "mesh": 320,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012931924717677148,
      "constant_only_rge101_L1_bound": 0.0001293192146424482,
      "combined_amplitude_rge101_L1_bound": 0.98523694734756,
      "positive_mode_gain_to_prior_production_start": 57924974240.5972,
      "occupation_cap": 21441710365.93471,
      "rearrangement_transition_r": 1509616347.5374174,
      "lambda_at_z100": 6.218111138985368e-18,
      "peak_abs_x_history": 1.9627452093810485e-05,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 1e-05,
      "mesh": 320,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012931924717677148,
      "constant_only_rge101_L1_bound": 0.0001293192146424482,
      "combined_amplitude_rge101_L1_bound": 119.1103826972792,
      "positive_mode_gain_to_prior_production_start": 7002895896495.527,
      "occupation_cap": 21441710365.93471,
      "rearrangement_transition_r": 1509616347.5374174,
      "lambda_at_z100": 1.0861223917530709e-17,
      "peak_abs_x_history": 2.5941745980524912e-05,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    },
    {
      "seed": 0.0001,
      "mesh": 320,
      "n": 12000,
      "constant_only_all_r_L1_bound": 0.00012931924717677148,
      "constant_only_rge101_L1_bound": 0.0001293192146424482,
      "combined_amplitude_rge101_L1_bound": 1188647.3952133728,
      "positive_mode_gain_to_prior_production_start": 6.988455617765458e+16,
      "occupation_cap": 21441710365.93471,
      "rearrangement_transition_r": 1509616347.5374174,
      "lambda_at_z100": 1.8518450713925466e-16,
      "peak_abs_x_history": 0.00010712573723895841,
      "scope": "Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes."
    }
  ]
}
```
<!-- END FILE -->

### bias-late32.json

<!-- FILE: bias-late32.json -->
```json
[
  {
    "vev_factor": 6.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 3.746105861179259e-05,
    "target_z": 3.2,
    "z_half": 3.1999613995056224,
    "max_extra_zle100": 0.03126137991741585,
    "positive_local_root_integral": 1.5992929670396392,
    "minimum_curvature_after_half": -794.2185759937413,
    "late_contact_proxy": 0.005269279396626546,
    "early_contact_proxy_100_to_1e5": 1.5550264370917033e-24,
    "early_representative_pair_proxy": 1.5453204968044202e-24,
    "early_x_z1e5": 9.99995481078295e-09,
    "max_abs_x": 1.2802569934505472,
    "today_x": 0.9824510982130213,
    "activation_integral": 4.9120513346435206e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 8.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 0.0003191830151647568,
    "target_z": 3.2,
    "z_half": 3.1999613995056224,
    "max_extra_zle100": 0.02978498670495141,
    "positive_local_root_integral": 0.8192432927900579,
    "minimum_curvature_after_half": -372.41801644742395,
    "late_contact_proxy": 0.005526848343011549,
    "early_contact_proxy_100_to_1e5": 1.9477132050309212e-22,
    "early_representative_pair_proxy": 1.947615488099219e-22,
    "early_x_z1e5": 9.999974636472981e-09,
    "max_abs_x": 1.270887636329997,
    "today_x": 1.0127607787575512,
    "activation_integral": 4.881109862988473e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 10.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 0.0011805813334504938,
    "target_z": 3.2,
    "z_half": 3.1999613995056224,
    "max_extra_zle100": 0.028496910847557542,
    "positive_local_root_integral": 0.4891806093783537,
    "minimum_curvature_after_half": -176.530201538893,
    "late_contact_proxy": 0.0056599941524043745,
    "early_contact_proxy_100_to_1e5": 6.0934210661552874e-21,
    "early_representative_pair_proxy": 6.0934112803104886e-21,
    "early_x_z1e5": 9.999983875793544e-09,
    "max_abs_x": 1.2620289077961653,
    "today_x": 1.0367178351482274,
    "activation_integral": 4.936701178081746e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  }
]
```
<!-- END FILE -->

### bias-late35.json

<!-- FILE: bias-late35.json -->
```json
[
  {
    "vev_factor": 6.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 9.338521785529019e-05,
    "target_z": 3.5,
    "z_half": 3.499550093326829,
    "max_extra_zle100": 0.025411033900905326,
    "positive_local_root_integral": 1.3307325682058249,
    "minimum_curvature_after_half": -623.7834176850814,
    "late_contact_proxy": 0.006449762851932295,
    "early_contact_proxy_100_to_1e5": 1.475395781767655e-23,
    "early_representative_pair_proxy": 1.4744239693941882e-23,
    "early_x_z1e5": 9.999954830341811e-09,
    "max_abs_x": 1.2771335633922256,
    "today_x": 0.9618497113358746,
    "activation_integral": 4.906034647322819e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 8.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 0.000626132965162596,
    "target_z": 3.5,
    "z_half": 3.499550093326829,
    "max_extra_zle100": 0.024114291145544493,
    "positive_local_root_integral": 0.6631220066018523,
    "minimum_curvature_after_half": -275.3583805493867,
    "late_contact_proxy": 0.006694745961431656,
    "early_contact_proxy_100_to_1e5": 2.8110298814441683e-21,
    "early_representative_pair_proxy": 2.811020154509639e-21,
    "early_x_z1e5": 9.999974696858703e-09,
    "max_abs_x": 1.26700151214899,
    "today_x": 0.973829765753506,
    "activation_integral": 4.931515349684038e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 10.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 0.002042548439903863,
    "target_z": 3.5,
    "z_half": 3.499550093326829,
    "max_extra_zle100": 0.02300786170835057,
    "positive_local_root_integral": 0.35096298325031816,
    "minimum_curvature_after_half": -115.23087828483374,
    "late_contact_proxy": 0.006910881020958008,
    "early_contact_proxy_100_to_1e5": 5.423635732198151e-20,
    "early_representative_pair_proxy": 5.423634655062988e-20,
    "early_x_z1e5": 9.999983984320477e-09,
    "max_abs_x": 1.257550229633473,
    "today_x": 0.9617037354876904,
    "activation_integral": 4.950334192124894e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  }
]
```
<!-- END FILE -->

### bias-late38.json

<!-- FILE: bias-late38.json -->
```json
[
  {
    "vev_factor": 6.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 0.000200410332337724,
    "target_z": 3.8,
    "z_half": 3.799893428172381,
    "max_extra_zle100": 0.020856197138952634,
    "positive_local_root_integral": 1.0764253689192667,
    "minimum_curvature_after_half": -489.4534341366171,
    "late_contact_proxy": 0.007673599920230731,
    "early_contact_proxy_100_to_1e5": 2.8013950006769155e-22,
    "early_representative_pair_proxy": 2.8012989808685603e-22,
    "early_x_z1e5": 9.99995486777268e-09,
    "max_abs_x": 1.2739416984017884,
    "today_x": 1.0319701559249967,
    "activation_integral": 4.968276942560356e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 8.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 0.0011085805186184896,
    "target_z": 3.8,
    "z_half": 3.799893428172381,
    "max_extra_zle100": 0.01972597550437843,
    "positive_local_root_integral": 0.532154497398876,
    "minimum_curvature_after_half": -199.6035961303608,
    "late_contact_proxy": 0.008004359232810176,
    "early_contact_proxy_100_to_1e5": 2.7256640187081325e-20,
    "early_representative_pair_proxy": 2.725663002834991e-20,
    "early_x_z1e5": 9.999974791769783e-09,
    "max_abs_x": 1.2630676063623405,
    "today_x": 1.0145081233129327,
    "activation_integral": 4.973955786403827e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 10.0,
    "tuned": true,
    "seed": 1e-08,
    "bias": 0.0032765520306085277,
    "target_z": 3.8,
    "z_half": 3.799893428172381,
    "max_extra_zle100": 0.018798096416401823,
    "positive_local_root_integral": 0.22157938904376978,
    "minimum_curvature_after_half": -67.91903215600301,
    "late_contact_proxy": 0.008379641243694214,
    "early_contact_proxy_100_to_1e5": 3.602592912309753e-19,
    "early_representative_pair_proxy": 3.6025927397922043e-19,
    "early_x_z1e5": 9.999984139689059e-09,
    "max_abs_x": 1.2530861690045436,
    "today_x": 1.0238711134981502,
    "activation_integral": 4.925513462838452e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  }
]
```
<!-- END FILE -->

### bias-local-modes.json

<!-- FILE: bias-local-modes.json -->
```json
[
  {
    "vev_factor": 1.0,
    "bias": 0.0,
    "order": 20,
    "ns": 480,
    "nk": 100,
    "z_half_sample": 4.295402197021065,
    "max_gamma_over_H": 126.29540630600347,
    "k_peak_integral": 3.053855508833412e-05,
    "max_frozen_positive_integral_all": 43.70953312872516,
    "max_frozen_positive_integral_post_half": 13.514017105699212,
    "max_unforced_post_half": 13.509723125661043,
    "max_eta_increment_post_half": 0.8735052249301437,
    "max_fast_subhorizon_integral": 3.5116594521059716,
    "min_curvature_after_half": -16329.417091098203,
    "k": [
      1e-05,
      1.1497569953977357e-05,
      1.3219411484660286e-05,
      1.5199110829529332e-05,
      1.747528400007683e-05,
      2.0092330025650458e-05,
      2.310129700083158e-05,
      2.6560877829466838e-05,
      3.053855508833412e-05,
      3.511191734215127e-05,
      4.037017258596558e-05,
      4.641588833612782e-05,
      5.3366992312063123e-05,
      6.135907273413175e-05,
      7.054802310718646e-05,
      8.111308307896872e-05,
      9.326033468832199e-05,
      0.00010722672220103231,
      0.0001232846739442066,
      0.00014174741629268049,
      0.00016297508346206434,
      0.0001873817422860383,
      0.00021544346900318845,
      0.0002477076355991711,
      0.0002848035868435802,
      0.00032745491628777284,
      0.00037649358067924675,
      0.00043287612810830614,
      0.0004977023564332113,
      0.0005722367659350221,
      0.0006579332246575682,
      0.000756463327554629,
      0.0008697490026177834,
      0.001,
      0.0011497569953977356,
      0.0013219411484660286,
      0.0015199110829529332,
      0.0017475284000076847,
      0.002009233002565048,
      0.0023101297000831605,
      0.0026560877829466868,
      0.0030538555088334154,
      0.003511191734215131,
      0.004037017258596553,
      0.004641588833612782,
      0.005336699231206312,
      0.006135907273413175,
      0.007054802310718645,
      0.008111308307896872,
      0.0093260334688322,
      0.010722672220103232,
      0.012328467394420659,
      0.014174741629268047,
      0.01629750834620645,
      0.01873817422860385,
      0.021544346900318843,
      0.024770763559917114,
      0.02848035868435802,
      0.03274549162877728,
      0.037649358067924674,
      0.043287612810830614,
      0.049770235643321135,
      0.0572236765935022,
      0.06579332246575682,
      0.07564633275546291,
      0.08697490026177834,
      0.1,
      0.11497569953977356,
      0.13219411484660287,
      0.1519911082952933,
      0.1747528400007683,
      0.2009233002565046,
      0.2310129700083158,
      0.26560877829466895,
      0.30538555088334185,
      0.35111917342151344,
      0.4037017258596558,
      0.4641588833612782,
      0.5336699231206312,
      0.6135907273413176,
      0.7054802310718645,
      0.8111308307896873,
      0.9326033468832199,
      1.072267222010323,
      1.232846739442066,
      1.4174741629268048,
      1.6297508346206435,
      1.873817422860383,
      2.1544346900318865,
      2.4770763559917137,
      2.848035868435805,
      3.2745491628777317,
      3.7649358067924714,
      4.328761281083062,
      4.9770235643321135,
      5.72236765935022,
      6.5793322465756825,
      7.56463327554629,
      8.697490026177835,
      10.0
    ],
    "all_integral": [
      43.70946261603528,
      43.709473448228145,
      43.709484885597696,
      43.70949646105158,
      43.70950803451699,
      43.709518547364034,
      43.70952722566546,
      43.709532697066614,
      43.70953312872516,
      43.7095256506652,
      43.70950685644337,
      43.70947126652224,
      43.70941225397588,
      43.70932023693493,
      43.709182725506786,
      43.70898245766663,
      43.70869657444575,
      43.708294316861995,
      43.70773452257034,
      43.706962184603285,
      43.705903875957794,
      43.70446166027074,
      43.70250507756423,
      43.69985992299256,
      43.69629429550284,
      43.69149872122391,
      43.68505981980816,
      43.67642427608113,
      43.66484819183116,
      43.649325721738414,
      43.62926449492294,
      43.603144943607326,
      43.56813699001671,
      43.521068265085866,
      43.45772024451682,
      43.37668726462908,
      43.26697578379951,
      43.12405605596608,
      42.93552907934605,
      42.6882981130096,
      42.36471162302247,
      41.94208877835605,
      41.39674552113604,
      40.68997969103535,
      39.787447765906585,
      38.64045887039888,
      37.198276204174604,
      35.40682172533424,
      33.21285004888203,
      30.5706093675422,
      27.460704337402625,
      23.894676208287173,
      19.96519201747451,
      16.004695976758892,
      13.060111885264186,
      11.025840381644407,
      8.817502001781897,
      6.348106904718059,
      3.7145020527392227,
      0.9241768471425619,
      0.19194152861585237,
      0.040939227673089815,
      0.007131124016737793,
      0.0009457524566113228,
      9.27925219674948e-05,
      6.111416346206631e-06,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "post_half_integral": [
      13.509837454849361,
      13.50985433123308,
      13.509873792015682,
      13.509896012303376,
      13.509921701728805,
      13.509950925804302,
      13.509984397670836,
      13.510022711984004,
      13.510066637586508,
      13.510116745441778,
      13.510174179870825,
      13.510239479697411,
      13.510313977135056,
      13.510398616102638,
      13.510494827075954,
      13.510603758433598,
      13.510726901516646,
      13.510865707074267,
      13.511021624583162,
      13.511196031257226,
      13.511390127850479,
      13.51160478829324,
      13.511840447065163,
      13.51209641492312,
      13.512371024426532,
      13.51266070034715,
      13.512959267918593,
      13.513256996307735,
      13.513538754155922,
      13.513782644346291,
      13.513956994789888,
      13.514017105699212,
      13.513899688165603,
      13.513517682030884,
      13.512750505006426,
      13.511433207038891,
      13.50934106499818,
      13.506166902666411,
      13.50149553888845,
      13.494764894634258,
      13.485215699005174,
      13.471823603349756,
      13.453208548626039,
      13.427512687840432,
      13.392235281151702,
      13.34401009831894,
      13.278298241104933,
      13.188960402516843,
      13.06764498666929,
      12.902884512190804,
      12.678703317213756,
      12.372344485437832,
      11.950248969389044,
      11.360080407174447,
      10.513162900334153,
      9.28769073219224,
      7.624506987292627,
      5.614810714541194,
      3.380809004160514,
      0.8735052249301437,
      0.1862289225837175,
      0.04053024579902959,
      0.007109977793108864,
      0.0009440467413458453,
      9.26638154686606e-05,
      6.111416346206631e-06,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "post_half_unforced": [
      13.509723125661043,
      13.50972287990382,
      13.509722555027228,
      13.509722125559477,
      13.50972155782837,
      13.509720807321205,
      13.509719815194837,
      13.509718503662032,
      13.509716769892629,
      13.509714477951121,
      13.50971144813866,
      13.509707442903698,
      13.509702148216721,
      13.509695148948506,
      13.509685896321491,
      13.50967366488203,
      13.509657495619507,
      13.509636120772182,
      13.509607864422772,
      13.509570511088356,
      13.509521131997827,
      13.50945585543031,
      13.509369563097716,
      13.509255488747836,
      13.50910468748499,
      13.508905334144343,
      13.508641795611652,
      13.508293404183007,
      13.507832835491689,
      13.507223963301454,
      13.506419022058493,
      13.505354853134737,
      13.503947937653813,
      13.502087821555383,
      13.499628408820296,
      13.496376425183783,
      13.492076121492943,
      13.486388971120904,
      13.478866688267749,
      13.468915308298644,
      13.455747260541532,
      13.43831722663485,
      13.415235955059257,
      13.384653839374952,
      13.344102536936274,
      13.290277463771522,
      13.218735298106616,
      13.12346606741021,
      12.996273694510531,
      12.825850410988048,
      12.596331097485924,
      12.284886300168594,
      11.857305152239677,
      11.258443619548206,
      10.38282365012238,
      8.991253638704164,
      7.123912017521213,
      4.894213776392046,
      2.8018651672726813,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "fast_integral": [
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      3.5116594521059716,
      2.2961175765956257,
      2.292162730750568,
      2.486530782937942,
      2.446279041162603,
      2.381270986335595,
      2.2585425329363367,
      2.2422759121731612,
      2.220595268745975,
      2.1916162223588973,
      2.1527290943827273,
      1.0339769709690336,
      1.000226889849598,
      0.9538693674276065,
      0.8890950956316671,
      0.7959938366003367,
      0.6550299392711593,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "scope": "Positive frozen roots can overestimate or underestimate true growth across time-dependent phases. No metric feedback or actual time evolution in this diagnostic."
  },
  {
    "vev_factor": 10.0,
    "bias": 0.0032765520306085277,
    "order": 20,
    "ns": 480,
    "nk": 100,
    "z_half_sample": 3.762891565693618,
    "max_gamma_over_H": 10.573713633210257,
    "k_peak_integral": 1.3219411484660286e-05,
    "max_frozen_positive_integral_all": 3.7867550482103236,
    "max_frozen_positive_integral_post_half": 0.22089732439567208,
    "max_unforced_post_half": 0.22021520046748339,
    "max_eta_increment_post_half": 0.039590306782110486,
    "max_fast_subhorizon_integral": 0.0,
    "min_curvature_after_half": -66.73887402673338,
    "k": [
      1e-05,
      1.1497569953977357e-05,
      1.3219411484660286e-05,
      1.5199110829529332e-05,
      1.747528400007683e-05,
      2.0092330025650458e-05,
      2.310129700083158e-05,
      2.6560877829466838e-05,
      3.053855508833412e-05,
      3.511191734215127e-05,
      4.037017258596558e-05,
      4.641588833612782e-05,
      5.3366992312063123e-05,
      6.135907273413175e-05,
      7.054802310718646e-05,
      8.111308307896872e-05,
      9.326033468832199e-05,
      0.00010722672220103231,
      0.0001232846739442066,
      0.00014174741629268049,
      0.00016297508346206434,
      0.0001873817422860383,
      0.00021544346900318845,
      0.0002477076355991711,
      0.0002848035868435802,
      0.00032745491628777284,
      0.00037649358067924675,
      0.00043287612810830614,
      0.0004977023564332113,
      0.0005722367659350221,
      0.0006579332246575682,
      0.000756463327554629,
      0.0008697490026177834,
      0.001,
      0.0011497569953977356,
      0.0013219411484660286,
      0.0015199110829529332,
      0.0017475284000076847,
      0.002009233002565048,
      0.0023101297000831605,
      0.0026560877829466868,
      0.0030538555088334154,
      0.003511191734215131,
      0.004037017258596553,
      0.004641588833612782,
      0.005336699231206312,
      0.006135907273413175,
      0.007054802310718645,
      0.008111308307896872,
      0.0093260334688322,
      0.010722672220103232,
      0.012328467394420659,
      0.014174741629268047,
      0.01629750834620645,
      0.01873817422860385,
      0.021544346900318843,
      0.024770763559917114,
      0.02848035868435802,
      0.03274549162877728,
      0.037649358067924674,
      0.043287612810830614,
      0.049770235643321135,
      0.0572236765935022,
      0.06579332246575682,
      0.07564633275546291,
      0.08697490026177834,
      0.1,
      0.11497569953977356,
      0.13219411484660287,
      0.1519911082952933,
      0.1747528400007683,
      0.2009233002565046,
      0.2310129700083158,
      0.26560877829466895,
      0.30538555088334185,
      0.35111917342151344,
      0.4037017258596558,
      0.4641588833612782,
      0.5336699231206312,
      0.6135907273413176,
      0.7054802310718645,
      0.8111308307896873,
      0.9326033468832199,
      1.072267222010323,
      1.232846739442066,
      1.4174741629268048,
      1.6297508346206435,
      1.873817422860383,
      2.1544346900318865,
      2.4770763559917137,
      2.848035868435805,
      3.2745491628777317,
      3.7649358067924714,
      4.328761281083062,
      4.9770235643321135,
      5.72236765935022,
      6.5793322465756825,
      7.56463327554629,
      8.697490026177835,
      10.0
    ],
    "all_integral": [
      3.786738491540557,
      3.7867492955518838,
      3.7867550482103236,
      3.7867531980783875,
      3.7867396191336096,
      3.786708983346067,
      3.786653870112926,
      3.786564174085232,
      3.786426196071929,
      3.786221535853034,
      3.78592509869812,
      3.785503548567681,
      3.7849122273298925,
      3.7840912697857187,
      3.782961091003733,
      3.7814154222028824,
      3.779313893185774,
      3.7764700736262373,
      3.772638184952449,
      3.7674944306738873,
      3.760613799722795,
      3.7514417710658696,
      3.739258047585584,
      3.723134024485066,
      3.7018836204747005,
      3.674008565351224,
      3.6376450062091203,
      3.5905162294481388,
      3.5299046590323595,
      3.452656319405248,
      3.3552397094569786,
      3.2338846539515362,
      3.0848240547417234,
      2.9046329457384483,
      2.6906358430492334,
      2.4415325400224064,
      2.158595107821498,
      1.8453295634243392,
      1.50661192305111,
      1.1554636300214638,
      0.8031770165291426,
      0.4493770998469868,
      0.14477891472787027,
      0.028404792188616912,
      0.005750022492546037,
      0.0011307496825666372,
      0.0001995636727217173,
      2.6363476941258702e-05,
      1.805031895402859e-06,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "post_half_integral": [
      0.2203205584588948,
      0.220335181922823,
      0.22035162357317664,
      0.22037033532184405,
      0.22039135103340987,
      0.2204149936399722,
      0.22044148818534826,
      0.2204710387652281,
      0.2205038076253868,
      0.22053988517937995,
      0.22057924778570637,
      0.22062169881990598,
      0.22066688384267727,
      0.22071390312114042,
      0.22076151801694194,
      0.22080737544834234,
      0.2208487646832582,
      0.22088080122152526,
      0.22089732439567208,
      0.22088924869642512,
      0.22084363731162676,
      0.22074277471824014,
      0.22056225953374256,
      0.22026838060053558,
      0.2198148642724704,
      0.21913832221112475,
      0.21815201743000082,
      0.21673727368265608,
      0.21473162273405313,
      0.21191147934424012,
      0.20796830505628908,
      0.20247198743293426,
      0.19481621038026722,
      0.1841391462472596,
      0.16924949006556617,
      0.14886402683726896,
      0.1227670149240158,
      0.09203386660880614,
      0.05927522802084154,
      0.035025357253995894,
      0.022280179168033272,
      0.014264737281962835,
      0.008570095344377913,
      0.004558746493638444,
      0.0020140083619980854,
      0.0006843574937833653,
      0.0001634761128240026,
      2.4710939253627046e-05,
      1.805031895402859e-06,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "post_half_unforced": [
      0.22021520046748339,
      0.22021389893879015,
      0.2202121783842725,
      0.22020990389466366,
      0.22020689712216165,
      0.22020292229152355,
      0.2201976677043727,
      0.22019072128342193,
      0.22018153823365505,
      0.2201693982753202,
      0.22015334907875445,
      0.2201321314365996,
      0.22010408025501096,
      0.22006699351089312,
      0.22001795874285768,
      0.21995312319531685,
      0.21986738910893122,
      0.21975400941568235,
      0.21960405064375144,
      0.21940567828921503,
      0.2191432039803948,
      0.2187958115095876,
      0.21833584724319238,
      0.2177265147789587,
      0.2169187461490737,
      0.21584691895464406,
      0.21442292666660426,
      0.21252784335768404,
      0.20999996658619005,
      0.20661718905070708,
      0.2020700245159524,
      0.19591815870181695,
      0.18751518120306895,
      0.1758633509191527,
      0.15928014223489026,
      0.13432569239929923,
      0.09893944424482805,
      0.059566626409171854,
      0.019684921238731057,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "fast_integral": [
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "scope": "Positive frozen roots can overestimate or underestimate true growth across time-dependent phases. No metric feedback or actual time evolution in this diagnostic."
  },
  {
    "vev_factor": 10.0,
    "bias": 0.0032765520306085277,
    "order": 28,
    "ns": 720,
    "nk": 160,
    "z_half_sample": 3.788361689165205,
    "max_gamma_over_H": 10.583541099319634,
    "k_peak_integral": 1.2977993346475677e-05,
    "max_frozen_positive_integral_all": 3.792377182620672,
    "max_frozen_positive_integral_post_half": 0.22039545593462798,
    "max_unforced_post_half": 0.2196594347744514,
    "max_eta_increment_post_half": 0.04282463124673804,
    "max_fast_subhorizon_integral": 0.0,
    "min_curvature_after_half": -66.99411290889468,
    "k": [
      1e-05,
      1.0907766915109564e-05,
      1.1897937907435857e-05,
      1.2977993346475677e-05,
      1.4156092644919915e-05,
      1.5441135899948328e-05,
      1.684283113011666e-05,
      1.8371767615786392e-05,
      2.0039495897155573e-05,
      2.185861503424674e-05,
      2.384286778806726e-05,
      2.6007244442001166e-05,
      2.8368096047762683e-05,
      3.094325795144362e-05,
      3.3752184532845694e-05,
      3.68160961760047e-05,
      4.01581395812115e-05,
      4.380356262962907e-05,
      4.777990512153967e-05,
      5.211720682918044e-05,
      5.684823443592555e-05,
      6.20087290762581e-05,
      6.76376763465999e-05,
      7.377760082683312e-05,
      8.047488733750888e-05,
      8.778013135972489e-05,
      9.574852126495768e-05,
      0.00010444025524245698,
      0.00011392099607392682,
      0.00012426236719115054,
      0.00013554249376408254,
      0.0001478465929071304,
      0.0001612676174624066,
      0.0001759069582234984,
      0.00019187520990478324,
      0.00020929300664290976,
      0.00022829193334233326,
      0.0002490155197497901,
      0.0002716203247675567,
      0.000296277119197087,
      0.00032317217584819515,
      0.00035250867676009135,
      0.00038450824816527683,
      0.0004194126347923944,
      0.00045748552615674007,
      0.0004990145486353981,
      0.0005443134383763516,
      0.0005937244114571097,
      0.0006476207491984753,
      0.000706409618164559,
      0.000770535146153055,
      0.0008404817774137396,
      0.0009167779324426061,
      0.001,
      0.0010907766915109554,
      0.0011897937907435858,
      0.0012977993346475663,
      0.0014156092644919916,
      0.0015441135899948314,
      0.001684283113011666,
      0.0018371767615786374,
      0.0020039495897155573,
      0.002185861503424672,
      0.0023842867788067264,
      0.0026007244442001138,
      0.002836809604776268,
      0.0030943257951443585,
      0.0033752184532845697,
      0.0036816096176004666,
      0.00401581395812115,
      0.0043803562629629025,
      0.004777990512153967,
      0.0052117206829180385,
      0.005684823443592555,
      0.006200872907625803,
      0.006763767634659991,
      0.007377760082683305,
      0.008047488733750887,
      0.00877801313597248,
      0.009574852126495767,
      0.010444025524245688,
      0.011392099607392682,
      0.012426236719115042,
      0.013554249376408254,
      0.014784659290713025,
      0.01612676174624066,
      0.017590695822349824,
      0.019187520990478324,
      0.020929300664290953,
      0.022829193334233326,
      0.024901551974978984,
      0.02716203247675567,
      0.029627711919708668,
      0.032317217584819515,
      0.035250867676009094,
      0.03845082481652768,
      0.04194126347923944,
      0.04574855261567406,
      0.04990145486353981,
      0.05443134383763522,
      0.05937244114571097,
      0.06476207491984753,
      0.0706409618164559,
      0.07705351461530549,
      0.08404817774137396,
      0.0916777932442606,
      0.1,
      0.10907766915109542,
      0.11897937907435857,
      0.1297799334647565,
      0.14156092644919915,
      0.15441135899948297,
      0.1684283113011666,
      0.18371767615786355,
      0.20039495897155574,
      0.21858615034246698,
      0.23842867788067262,
      0.2600724444200111,
      0.2836809604776268,
      0.30943257951443554,
      0.33752184532845697,
      0.3681609617600463,
      0.40158139581211505,
      0.4380356262962898,
      0.4777990512153967,
      0.5211720682918034,
      0.5684823443592555,
      0.6200872907625796,
      0.6763767634659991,
      0.7377760082683312,
      0.8047488733750887,
      0.8778013135972489,
      0.9574852126495768,
      1.04440255242457,
      1.1392099607392683,
      1.2426236719115056,
      1.3554249376408254,
      1.478465929071304,
      1.612676174624066,
      1.7590695822349842,
      1.9187520990478324,
      2.0929300664290977,
      2.2829193334233326,
      2.490155197497901,
      2.716203247675567,
      2.96277119197087,
      3.2317217584819518,
      3.5250867676009134,
      3.845082481652768,
      4.194126347923945,
      4.574855261567401,
      4.990145486353981,
      5.443134383763517,
      5.937244114571097,
      6.476207491984746,
      7.06409618164559,
      7.705351461530541,
      8.404817774137397,
      9.167779324426052,
      10.0
    ],
    "all_integral": [
      3.792364956304476,
      3.792370630801789,
      3.7923747932121814,
      3.792377182620672,
      3.792376702919007,
      3.7923726755919778,
      3.7923641783901094,
      3.7923500116059783,
      3.7923287115453324,
      3.7922984960293884,
      3.7922571986487337,
      3.7922022586930866,
      3.7921304876413413,
      3.7920378409278435,
      3.791919771480863,
      3.791770642323041,
      3.791583787550402,
      3.7913507495162753,
      3.791061782877962,
      3.7907050840509093,
      3.7902662684142143,
      3.789728299235946,
      3.789070758272944,
      3.7882691388797207,
      3.7872939415465448,
      3.7861103144460353,
      3.784676559336595,
      3.7829430234374533,
      3.7808506563317636,
      3.778329459865357,
      3.775295941309091,
      3.771651831774164,
      3.7672806084464665,
      3.7620448574194794,
      3.7557828872117835,
      3.748304300704957,
      3.739386732144386,
      3.728770144525208,
      3.716152318336215,
      3.701183177073652,
      3.6834595561644785,
      3.662519687485154,
      3.6378384266915025,
      3.6088229953400814,
      3.574811121489698,
      3.5350704683422896,
      3.4888010530476037,
      3.4351427365488725,
      3.373186256063313,
      3.301991216054946,
      3.220609613847438,
      3.128118266827665,
      3.023661090114551,
      2.906509460459674,
      2.7761572157505134,
      2.6324484446827645,
      2.4756040280126452,
      2.306016255922601,
      2.1241419591906374,
      1.930425895563258,
      1.7255606124535208,
      1.5120389859398187,
      1.2945119844708126,
      1.0764634358408958,
      0.8567758760455462,
      0.6353710438082675,
      0.41690837882513937,
      0.21730421816399773,
      0.08357398672771654,
      0.030181520505565784,
      0.011162747358384864,
      0.0041193880854242665,
      0.0014969898770842355,
      0.0005266690906652832,
      0.00017358934165928193,
      5.0998187912816237e-05,
      1.2310290760481929e-05,
      2.1991705082734542e-06,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "post_half_integral": [
      0.21973004172656485,
      0.21973628034766668,
      0.21974316106100605,
      0.21975095493624328,
      0.21975941959876857,
      0.2197687021938044,
      0.2197789556335751,
      0.2197902824551855,
      0.21980279367058628,
      0.2198166082560527,
      0.21983185200614677,
      0.21984872439599407,
      0.2198673569856009,
      0.21988769454957163,
      0.2199099790823409,
      0.21993431767667027,
      0.21996092836482345,
      0.2199895972464885,
      0.22002044393386863,
      0.22005339653306671,
      0.22008829730384208,
      0.22012488121356086,
      0.22016288545171425,
      0.220201632844848,
      0.22024024316355878,
      0.22027786660121818,
      0.22031326158370618,
      0.22034486887193377,
      0.22037075343158283,
      0.2203886718423146,
      0.22039545593462798,
      0.22038769716636517,
      0.2203611825483376,
      0.220310621088901,
      0.22023005234899767,
      0.22011160357774323,
      0.21994653542613402,
      0.2197239138061675,
      0.21943083633926236,
      0.21905173635374148,
      0.2185679173887143,
      0.21795690608297838,
      0.2171917247365562,
      0.21623964625841446,
      0.21506136126653516,
      0.2136093371944339,
      0.21182549730749034,
      0.2096393929650034,
      0.20696475405345485,
      0.20369602025768313,
      0.1997032917057353,
      0.19482789460724886,
      0.18887889439444036,
      0.18163790036303065,
      0.1728857379397043,
      0.16245513432587927,
      0.15023035446277822,
      0.1360502087815217,
      0.11976215758177026,
      0.10121495573153294,
      0.08056003940279574,
      0.059812427746271005,
      0.04282463124673804,
      0.031456785037984936,
      0.02380346195292344,
      0.018135605163225373,
      0.013651636614779526,
      0.01000270057767965,
      0.00703667745626655,
      0.004685334849083834,
      0.0029067355235327694,
      0.0016506652843367347,
      0.0008411697462237672,
      0.00037609104997543856,
      0.00014437874935751014,
      4.6358647385647686e-05,
      1.1957649346679131e-05,
      2.1991705082734542e-06,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "post_half_unforced": [
      0.2196594347744514,
      0.21965816598915983,
      0.21965665636680598,
      0.21965486018574817,
      0.21965272304156053,
      0.2196501801969723,
      0.2196471546177517,
      0.21964355463456042,
      0.2196392711592675,
      0.2196341743704206,
      0.21962810976607217,
      0.21962089346239586,
      0.2196123065928283,
      0.21960208863400865,
      0.21958992945054406,
      0.21957545980934323,
      0.21955824006436078,
      0.2195377466521117,
      0.219513355964739,
      0.2194843250775618,
      0.2194497686977431,
      0.21940863156459425,
      0.21935965536287771,
      0.21930133899861468,
      0.21923189081909505,
      0.2191491710165815,
      0.21905062201237785,
      0.218933184036472,
      0.21879319234140213,
      0.21862625143138736,
      0.21842708021460774,
      0.21818931988091625,
      0.21790529320755916,
      0.21756569927589792,
      0.2171592201139185,
      0.21667200342826126,
      0.21608696409584777,
      0.21538280739309673,
      0.21453259829960517,
      0.21390469227416553,
      0.21320545109924038,
      0.2123707636666315,
      0.21137373699007683,
      0.21018184478358878,
      0.20875561108438043,
      0.20704691977884987,
      0.20499680926730685,
      0.20253253905117946,
      0.19956359229283757,
      0.19597605912626284,
      0.19162442719713843,
      0.18631893879962977,
      0.17980466971414805,
      0.17172310559651713,
      0.16152900922360652,
      0.14824868752273376,
      0.13173134735015748,
      0.11661877169532102,
      0.09495641508484445,
      0.0764170440345351,
      0.0489522931472351,
      0.01917875639116245,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "fast_integral": [
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "scope": "Positive frozen roots can overestimate or underestimate true growth across time-dependent phases. No metric feedback or actual time evolution in this diagnostic."
  }
]
```
<!-- END FILE -->

### bias-seed-wall.json

<!-- FILE: bias-seed-wall.json -->
```json
{
  "wall": {
    "vev_over_Mpl": 0.03716554114645386,
    "mu_eV": 1.5769136014685922e-31,
    "negative_local_minimum": -0.9933679185620146,
    "barrier": -0.0329061503382041,
    "delta_vacuum_eV4": 2.660336448884598e-12,
    "unbiased_wall_tension_eV3": 1.0206089636780571e+21,
    "pressure_balance_H_over_H0": 1.813017798989102,
    "pressure_balance_z_lcdm": 1.0210969002746813,
    "one_wall_per_hubble_density_fraction_at_balance": 0.022012177493029715,
    "one_wall_per_hubble_density_fraction_today": 0.039908469589370194,
    "bias_tadpole_eV3": 1.4746952758692566e-38,
    "bias_mass_squared_eV2": 1.6295318851574984e-64,
    "scope": "Pressure balance H~delta V/tension; one wall per Hubble volume dimensional estimate, not a network or a prediction that walls form."
  },
  "seeds": [
    {
      "seed": 1e-08,
      "z_half": 3.7994556129568933,
      "today_x": 1.0238711134981502,
      "sign_at_half": 1.0,
      "max_extra_zle100": 0.01879809729182447,
      "h0_ratio": 1.0030711878986183,
      "max_h_ratio_zle100": 1.0095336723804147,
      "continuity_rel": 1.1757367557725832e-05,
      "late_contact_proxy": 0.008379599570069503,
      "early_contact_proxy_100_1e5": 3.60259128666909e-19,
      "early_pair_proxy_100_1e5": 3.6025911141521085e-19,
      "high_contact_extrapolation_1e5_1e10": 1.2197949313713516e-14,
      "high_pair_NR_extrapolation_1e5_1e10": 1.4763499294158233e-22,
      "max_abs_x_zge100": 1.6859638115038405e-05,
      "activation_integral": 4.9255379635292475e+32,
      "scope": "High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds."
    },
    {
      "seed": 0.0,
      "z_half": 3.7994556129568933,
      "today_x": 1.0238711038596455,
      "sign_at_half": 1.0,
      "max_extra_zle100": 0.018798098088106476,
      "h0_ratio": 1.0030711880399255,
      "max_h_ratio_zle100": 1.009533672790052,
      "continuity_rel": 1.1757369114436376e-05,
      "late_contact_proxy": 0.008379599266257632,
      "early_contact_proxy_100_1e5": 3.591246596832972e-19,
      "early_pair_proxy_100_1e5": 3.5912465213832203e-19,
      "high_contact_extrapolation_1e5_1e10": 5.002020593636506e-55,
      "high_pair_NR_extrapolation_1e5_1e10": 4.907852711109451e-55,
      "max_abs_x_zge100": 1.685061627570351e-05,
      "activation_integral": 4.9255379546842845e+32,
      "scope": "High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds."
    },
    {
      "seed": -1e-05,
      "z_half": 3.7994556129568933,
      "today_x": 1.0238614654120584,
      "sign_at_half": 1.0,
      "max_extra_zle100": 0.0187988940691728,
      "h0_ratio": 1.003071329332345,
      "max_h_ratio_zle100": 1.0095340822726349,
      "continuity_rel": 1.1757932721778788e-05,
      "late_contact_proxy": 0.008379295561003594,
      "early_contact_proxy_100_1e5": 1.1823515963803602e-12,
      "early_pair_proxy_100_1e5": 1.1726719968288265e-12,
      "high_contact_extrapolation_1e5_1e10": 0.012197971229435625,
      "high_pair_NR_extrapolation_1e5_1e10": 1.4763525800680305e-10,
      "max_abs_x_zge100": 1e-05,
      "activation_integral": 4.925529113965726e+32,
      "scope": "High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds."
    },
    {
      "seed": 1e-05,
      "z_half": 3.7994556129568933,
      "today_x": 1.0238807367302445,
      "sign_at_half": 1.0,
      "max_extra_zle100": 0.01879730189703339,
      "h0_ratio": 1.0030710467758068,
      "max_h_ratio_zle100": 1.0095332631999319,
      "continuity_rel": 1.175679311757716e-05,
      "late_contact_proxy": 0.008379902980667677,
      "early_contact_proxy_100_1e5": 1.1823567386174323e-12,
      "early_pair_proxy_100_1e5": 1.1726771737804827e-12,
      "high_contact_extrapolation_1e5_1e10": 0.012197927431773224,
      "high_pair_NR_extrapolation_1e5_1e10": 1.4763472791253423e-10,
      "max_abs_x_zge100": 2.5872458056881117e-05,
      "activation_integral": 4.925546799483408e+32,
      "scope": "High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds."
    },
    {
      "seed": -0.0001,
      "z_half": 3.7976099752474353,
      "today_x": 1.023774466650897,
      "sign_at_half": 1.0,
      "max_extra_zle100": 0.018806048619504283,
      "h0_ratio": 1.0030726022674932,
      "max_h_ratio_zle100": 1.0095377628646454,
      "continuity_rel": 1.1762558470855839e-05,
      "late_contact_proxy": 0.008376562558978696,
      "early_contact_proxy_100_1e5": 1.182372929722336e-08,
      "early_pair_proxy_100_1e5": 1.1726931741446198e-08,
      "high_contact_extrapolation_1e5_1e10": 121.98167826446084,
      "high_pair_NR_extrapolation_1e5_1e10": 1.4763763747236018e-06,
      "max_abs_x_zge100": 0.0001,
      "activation_integral": 4.925449729560418e+32,
      "scope": "High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds."
    },
    {
      "seed": 0.0001,
      "z_half": 3.7994556129568933,
      "today_x": 1.0239671811877966,
      "sign_at_half": 1.0,
      "max_extra_zle100": 0.01879012672722356,
      "h0_ratio": 1.0030697766716918,
      "max_h_ratio_zle100": 1.0095295720498367,
      "continuity_rel": 1.175116271459439e-05,
      "late_contact_proxy": 0.008382636818291368,
      "early_contact_proxy_100_1e5": 1.1823311908904872e-08,
      "early_pair_proxy_100_1e5": 1.1726517828638281e-08,
      "high_contact_extrapolation_1e5_1e10": 121.97729846823567,
      "high_pair_NR_extrapolation_1e5_1e10": 1.4763233648968857e-06,
      "max_abs_x_zge100": 0.00010706925286472294,
      "activation_integral": 4.9256265859624186e+32,
      "scope": "High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds."
    },
    {
      "seed": -0.001,
      "z_half": 3.7920773194835977,
      "today_x": 1.0228790934925431,
      "sign_at_half": 1.0,
      "max_extra_zle100": 0.018877840676886658,
      "h0_ratio": 1.0030854703777157,
      "max_h_ratio_zle100": 1.009574697851037,
      "continuity_rel": 1.176345875140194e-05,
      "late_contact_proxy": 0.00834924741843122,
      "early_contact_proxy_100_1e5": 0.0001182559566317634,
      "early_pair_proxy_100_1e5": 0.00011728782852939883,
      "high_contact_extrapolation_1e5_1e10": 1220009.0150012895,
      "high_pair_NR_extrapolation_1e5_1e10": 0.014766090386073391,
      "max_abs_x_zge100": 0.001,
      "activation_integral": 4.924673644467944e+32,
      "scope": "High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds."
    },
    {
      "seed": 0.001,
      "z_half": 3.8049967872723274,
      "today_x": 1.0248064715283458,
      "sign_at_half": 1.0,
      "max_extra_zle100": 0.01871855029485763,
      "h0_ratio": 1.0030572007905052,
      "max_h_ratio_zle100": 1.0094927529278248,
      "continuity_rel": 1.1767190678660774e-05,
      "late_contact_proxy": 0.008410016143485893,
      "early_contact_proxy_100_1e5": 0.00011821350729001314,
      "early_pair_proxy_100_1e5": 0.00011724572676669243,
      "high_contact_extrapolation_1e5_1e10": 1219571.0051805372,
      "high_pair_NR_extrapolation_1e5_1e10": 0.014760789037867567,
      "max_abs_x_zge100": 0.001,
      "activation_integral": 4.926442344618336e+32,
      "scope": "High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds."
    },
    {
      "seed": -0.01,
      "z_half": 3.735278809167183,
      "today_x": 1.011433512191885,
      "sign_at_half": 1.0,
      "max_extra_zle100": 0.01962405263524883,
      "h0_ratio": 1.0032310147107726,
      "max_h_ratio_zle100": 1.0099588430428972,
      "continuity_rel": 1.125281979309523e-05,
      "late_contact_proxy": 0.008072230754739722,
      "early_contact_proxy_100_1e5": 1.1839461326512513,
      "early_pair_proxy_100_1e5": 1.1742535004869255,
      "high_contact_extrapolation_1e5_1e10": 12214394551.265474,
      "high_pair_NR_extrapolation_1e5_1e10": 147.83403379606233,
      "max_abs_x_zge100": 0.01,
      "activation_integral": 4.918906309251056e+32,
      "scope": "High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds."
    },
    {
      "seed": 0.01,
      "z_half": 3.8626315476720503,
      "today_x": 1.0304663358698094,
      "sign_at_half": 1.0,
      "max_extra_zle100": 0.01801060805453761,
      "h0_ratio": 1.0029424104229065,
      "max_h_ratio_zle100": 1.0091288022486145,
      "continuity_rel": 1.1743821881670427e-05,
      "late_contact_proxy": 0.008688460416201423,
      "early_contact_proxy_100_1e5": 1.179814068942845,
      "early_pair_proxy_100_1e5": 1.1701552648737124,
      "high_contact_extrapolation_1e5_1e10": 12171765294.826372,
      "high_pair_NR_extrapolation_1e5_1e10": 147.31808067937988,
      "max_abs_x_zge100": 0.01,
      "activation_integral": 4.936133828298445e+32,
      "scope": "High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds."
    },
    {
      "seed": -0.1,
      "z_half": 2.4108281361496013,
      "today_x": 0.9380295318252013,
      "sign_at_half": 1.0,
      "max_extra_zle100": 0.053875546920452676,
      "h0_ratio": 1.0099731355061274,
      "max_h_ratio_zle100": 1.0280775281815435,
      "continuity_rel": 1.1856570671803312e-05,
      "late_contact_proxy": 0.004072383847984544,
      "early_contact_proxy_100_1e5": 12199.36122077751,
      "early_pair_proxy_100_1e5": 12099.48850328899,
      "high_contact_extrapolation_1e5_1e10": 125856914276151.3,
      "high_pair_NR_extrapolation_1e5_1e10": 1523279.3766831195,
      "max_abs_x_zge100": 0.1,
      "activation_integral": 4.773791155412614e+32,
      "scope": "High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds."
    },
    {
      "seed": 0.1,
      "z_half": 4.524956667319193,
      "today_x": 0.976230292285621,
      "sign_at_half": 1.0,
      "max_extra_zle100": 0.011805635019509524,
      "h0_ratio": 1.0018748955673422,
      "max_h_ratio_zle100": 1.0059556019435492,
      "continuity_rel": 1.1214733342187533e-05,
      "late_contact_proxy": 0.013257902556246065,
      "early_contact_proxy_100_1e5": 11612.060938609866,
      "early_pair_proxy_100_1e5": 11516.996282363258,
      "high_contact_extrapolation_1e5_1e10": 119797924726451.05,
      "high_pair_NR_extrapolation_1e5_1e10": 1449945.8305869158,
      "max_abs_x_zge100": 0.1,
      "activation_integral": 5.0152339916649346e+32,
      "scope": "High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds."
    }
  ]
}
```
<!-- END FILE -->

### bias-window-0.json

<!-- FILE: bias-window-0.json -->
```json
[
  {
    "xfactor": 0.8,
    "vev_factor": 10.0,
    "target_z": 3.8,
    "bias": 0.0014710145900123108,
    "tuned": true,
    "max_extra": 0.016948282200962075,
    "max_loop_over_X": 0.9418552096276805,
    "cosmic_contact_proxy": 0.008471637291072712,
    "sigma0": 0.2220796617718809,
    "positive_local_root_integral": 0.19812870995154072,
    "min_curvature_after_half": -61.76456344295539,
    "max_abs_x": 1.2694304149391313,
    "vev_over_Mpl": 0.033241870569614515,
    "passes": true,
    "X": 1.629238400865156e-10,
    "m_eta": 0.01132781336250948,
    "CW": 0.16000000000000003
  },
  {
    "xfactor": 0.85,
    "vev_factor": 10.0,
    "target_z": 3.8,
    "bias": 0.0017996516500521193,
    "tuned": true,
    "max_extra": 0.017481755544424587,
    "max_loop_over_X": 0.8866390797053887,
    "cosmic_contact_proxy": 0.0084458903297833,
    "sigma0": 0.22217342753424824,
    "positive_local_root_integral": 0.20513289016375794,
    "min_curvature_after_half": -63.568175186671425,
    "max_abs_x": 1.265264627775706,
    "vev_over_Mpl": 0.03426493588790794,
    "passes": true,
    "X": 1.7310658009192282e-10,
    "m_eta": 0.011676442750227436,
    "CW": 0.17
  },
  {
    "xfactor": 0.9,
    "vev_factor": 10.0,
    "target_z": 3.8,
    "bias": 0.0022004752604101205,
    "tuned": true,
    "max_extra": 0.017965620674700818,
    "max_loop_over_X": 0.837531955769679,
    "cosmic_contact_proxy": 0.008421996638363058,
    "sigma0": 0.22225335033179147,
    "positive_local_root_integral": 0.21141567412195603,
    "min_curvature_after_half": -65.18187353874987,
    "max_abs_x": 1.2611589720853937,
    "vev_over_Mpl": 0.03525832814864992,
    "passes": true,
    "X": 1.8328932009733007e-10,
    "m_eta": 0.01201496046696906,
    "CW": 0.18000000000000002
  }
]
```
<!-- END FILE -->

### bias10-acoustic-fine-k.03.json

<!-- FILE: bias10-acoustic-fine-k.03.json -->
```json
{
  "k": 0.03,
  "control": "full",
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "order": 6,
  "angles": 200,
  "rtol": 5e-08,
  "zdec": 1090.0,
  "zmatch": 100.0,
  "baryon_fraction": 0.15555555555555556,
  "nu_fraction": 0.409,
  "amplitude": 1e-05,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "reached_chi_threshold": false,
  "max_eta": 0.008260753886014145,
  "max_chi_rel": 4.3612041556582705e-07,
  "eta_rest_at_stop": -0.008259392481749379,
  "max_matter": 0.08811974814243732,
  "matter_at_stop": -0.08811974814243732,
  "phi_at_stop": 2.3381982706399044e-06,
  "max_density_constraint_rel": 4.648387198638723e-14,
  "match_eta_number_per_psi": -0.8856597459372366,
  "match_matter_per_psi": -115.32914490707827,
  "scope": "Sudden decoupling, ideal tight coupling before z_dec, collisionless photons thereafter; no visibility function, polarization or full recombination. This is a sensitivity diagnostic."
}
```
<!-- END FILE -->

### bias10-acoustic-fine-k.05.json

<!-- FILE: bias10-acoustic-fine-k.05.json -->
```json
{
  "k": 0.05,
  "control": "full",
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "order": 5,
  "angles": 224,
  "rtol": 5e-08,
  "zdec": 1090.0,
  "zmatch": 100.0,
  "baryon_fraction": 0.15555555555555556,
  "nu_fraction": 0.409,
  "amplitude": 1e-05,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "reached_chi_threshold": false,
  "max_eta": 0.007913350189585748,
  "max_chi_rel": 9.107057017254492e-08,
  "eta_rest_at_stop": -0.007912853043273816,
  "max_matter": 0.1587328338600584,
  "matter_at_stop": -0.1587328338600584,
  "phi_at_stop": 1.5163415531801434e-06,
  "max_density_constraint_rel": 1.4468119365440887e-13,
  "match_eta_number_per_psi": -0.6139270795629026,
  "match_matter_per_psi": -207.03837814163523,
  "scope": "Sudden decoupling, ideal tight coupling before z_dec, collisionless photons thereafter; no visibility function, polarization or full recombination. This is a sensitivity diagnostic."
}
```
<!-- END FILE -->

### bias10-acoustic-fine-k.06.json

<!-- FILE: bias10-acoustic-fine-k.06.json -->
```json
{
  "k": 0.06,
  "control": "full",
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "order": 5,
  "angles": 256,
  "rtol": 5e-08,
  "zdec": 1090.0,
  "zmatch": 100.0,
  "baryon_fraction": 0.15555555555555556,
  "nu_fraction": 0.409,
  "amplitude": 1e-05,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "reached_chi_threshold": false,
  "max_eta": 0.0072155987124624945,
  "max_chi_rel": 3.588353537825827e-08,
  "eta_rest_at_stop": -0.007215276829058508,
  "max_matter": 0.184375450731183,
  "matter_at_stop": -0.184375450731183,
  "phi_at_stop": 1.223134499195067e-06,
  "max_density_constraint_rel": 9.96712653935988e-14,
  "match_eta_number_per_psi": -0.5074497833011785,
  "match_matter_per_psi": -240.3852835884915,
  "scope": "Sudden decoupling, ideal tight coupling before z_dec, collisionless photons thereafter; no visibility function, polarization or full recombination. This is a sensitivity diagnostic."
}
```
<!-- END FILE -->

### bias10-acoustic-k03.json

<!-- FILE: bias10-acoustic-k03.json -->
```json
{
  "k": 0.03,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "order": 5,
  "angles": 160,
  "rtol": 1e-07,
  "zdec": 1090.0,
  "zmatch": 100.0,
  "baryon_fraction": 0.15555555555555556,
  "nu_fraction": 0.409,
  "amplitude": 1e-05,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "reached_chi_threshold": false,
  "max_eta": 0.0082609768128561,
  "max_chi_rel": 4.361213124599871e-07,
  "eta_rest_at_stop": -0.008259614516103598,
  "max_matter": 0.08811976162003239,
  "max_density_constraint_rel": 3.9751448223307225e-13,
  "match_eta_number_per_psi": -0.8856598814067413,
  "match_matter_per_psi": -115.3291625476508,
  "scope": "Sudden decoupling, ideal tight coupling before z_dec, collisionless photons thereafter; no visibility function, polarization or full recombination. This is a sensitivity diagnostic."
}
```
<!-- END FILE -->

### bias10-chi_iso-k003.json

<!-- FILE: bias10-chi_iso-k003.json -->
```json
{
  "k": 0.003,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 100.0,
  "zstart": 1000000.0,
  "initial_mode": "chi_iso",
  "initial_amplitude": 1e-05,
  "psi_initial": 0.0,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 1.957993422274414e-05,
  "initial_radiation_delta_per_psi": 3.9453256476279136e-22,
  "match_delta_m_per_psi": -8.634423022243559e-08,
  "match_eta_delta_n_per_psi": -5.429692262537219e-08,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 6.798966027239725e-06,
  "max_matter_fraction": 3.279499843931256e-08,
  "max_phi": 3.732749354760408e-10,
  "max_delta_x": 1.041487074636297e-06,
  "max_delta_x_over_x_floor": 5.207435373181485e-06,
  "at_stop": {
    "k_over_aH": 13.303021872731664,
    "eta_number": 6.798966027239725e-06,
    "eta_number_rest_frame": 6.746778448639705e-06,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": 3.7072530551780756e-10,
    "delta_x_over_x": 1.0541946397575973e-08
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 5.455375885486006e-06,
  "early_nfev": 8381,
  "late_nfev": 23087,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like or compensated trial mode; no primordial spectrum imposed. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias10-early-k.0003.json

<!-- FILE: bias10-early-k.0003.json -->
```json
{
  "k": 0.0003,
  "radial_order": 5,
  "angles": 32,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 100.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 1.957993422274414e-05,
  "initial_radiation_delta_per_psi": -2.000000014809816,
  "match_delta_m_per_psi": -1.995548496957794,
  "match_eta_delta_n_per_psi": -1.9952537255538811,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 5.394687697376263e-05,
  "max_matter_fraction": 5.549575219796173e-05,
  "max_phi": 1.0117107897143378e-05,
  "max_delta_x": 5.952101257426843e-05,
  "max_delta_x_over_x_floor": 6.435310229237661e-05,
  "at_stop": {
    "k_over_aH": 1.330302187273166,
    "eta_number": -5.394687697376263e-05,
    "eta_number_rest_frame": -2.778826552157792e-05,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": -1.3151590714613722e-09,
    "delta_x_over_x": 2.9822083609166428e-05
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 1.7346924971569176e-08,
  "early_nfev": 146525,
  "late_nfev": 23087,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias10-early-k.001.json

<!-- FILE: bias10-early-k.001.json -->
```json
{
  "k": 0.001,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 100.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 1.957993422274414e-05,
  "initial_radiation_delta_per_psi": -2.0000001645535135,
  "match_delta_m_per_psi": -2.340684879392428,
  "match_eta_delta_n_per_psi": -2.3057837807168857,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.0002713063452470629,
  "max_matter_fraction": 0.00034623172884016353,
  "max_phi": 1.0036639942449093e-05,
  "max_delta_x": 6.220520189824563e-05,
  "max_delta_x_over_x_floor": 6.363010116268745e-05,
  "at_stop": {
    "k_over_aH": 4.434340624243888,
    "eta_number": -0.0002713063452470629,
    "eta_number_rest_frame": -0.00024472287218637077,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": -1.3182137600203186e-08,
    "delta_x_over_x": 2.9231469100564872e-05
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 1.7980576513709315e-09,
  "early_nfev": 53333,
  "late_nfev": 23087,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias10-early-k.003.json

<!-- FILE: bias10-early-k.003.json -->
```json
{
  "k": 0.003,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 100.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 1.957993422274414e-05,
  "initial_radiation_delta_per_psi": -2.000001480981622,
  "match_delta_m_per_psi": -5.260624124037303,
  "match_eta_delta_n_per_psi": -3.2830151029935846,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.002285282206497307,
  "max_matter_fraction": 0.0027158248705016572,
  "max_phi": 9.445452968781343e-06,
  "max_delta_x": 8.137102122391728e-05,
  "max_delta_x_over_x_floor": 8.521204113152598e-05,
  "at_stop": {
    "k_over_aH": 13.303021872731664,
    "eta_number": -0.002285282206497307,
    "eta_number_rest_frame": -0.0022529292438560364,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": -1.2387862334203072e-07,
    "delta_x_over_x": 1.004636331654952e-05
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 3.471297961842272e-10,
  "early_nfev": 22541,
  "late_nfev": 23087,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias10-early-k.006.json

<!-- FILE: bias10-early-k.006.json -->
```json
{
  "k": 0.006,
  "radial_order": 5,
  "angles": 64,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 100.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 1.957993422274414e-05,
  "initial_radiation_delta_per_psi": -2.0000059239264862,
  "match_delta_m_per_psi": -13.97638130514229,
  "match_eta_delta_n_per_psi": -2.3572888199025686,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.005279698736399854,
  "max_matter_fraction": 0.009551739920432234,
  "max_phi": 8.30003244675006e-06,
  "max_delta_x": 2.4156931873484713e-05,
  "max_delta_x_over_x_floor": 2.913439564447508e-05,
  "at_stop": {
    "k_over_aH": 26.606043745463328,
    "eta_number": -0.005279698736399854,
    "eta_number_rest_frame": -0.00526265071422933,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": -2.8773472321914063e-07,
    "delta_x_over_x": -5.592723569105626e-06
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 2.4543886932637162e-11,
  "early_nfev": 14006,
  "late_nfev": 23087,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias10-early-k.01.json

<!-- FILE: bias10-early-k.01.json -->
```json
{
  "k": 0.01,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 100.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 1.957993422274414e-05,
  "initial_radiation_delta_per_psi": -2.00001645535135,
  "match_delta_m_per_psi": -30.751782718845558,
  "match_eta_delta_n_per_psi": -2.030609882780203,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.00662621775512236,
  "max_matter_fraction": 0.02260929517027632,
  "max_phi": 7.094705967938269e-06,
  "max_delta_x": 5.595680154810113e-06,
  "max_delta_x_over_x_floor": 6.82642096800752e-06,
  "at_stop": {
    "k_over_aH": 44.34340624243888,
    "eta_number": -0.00662621775512236,
    "eta_number_rest_frame": -0.006617530775556206,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": -3.6084310236070593e-07,
    "delta_x_over_x": -2.4522930964991046e-06
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 3.5761341552041693e-12,
  "early_nfev": 10391,
  "late_nfev": 23087,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias10-early-k.02.json

<!-- FILE: bias10-early-k.02.json -->
```json
{
  "k": 0.02,
  "radial_order": 5,
  "angles": 128,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 100.0,
  "zstart": 1000000.0,
  "initial_mode": "adiabatic",
  "initial_amplitude": 1e-05,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 1.957993422274414e-05,
  "initial_radiation_delta_per_psi": -2.0000658214053995,
  "match_delta_m_per_psi": -83.97815475169192,
  "match_eta_delta_n_per_psi": -1.4350364451763165,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.00938103552023841,
  "max_matter_fraction": 0.06378898756028636,
  "max_phi": 5.010268970557125e-06,
  "max_delta_x": 1.4777436359205798e-06,
  "max_delta_x_over_x_floor": 1.7037731091757795e-06,
  "at_stop": {
    "k_over_aH": 88.68681248487776,
    "eta_number": -0.00938103552023841,
    "eta_number_rest_frame": -0.00937771190258827,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": -5.106688198790347e-07,
    "delta_x_over_x": -7.011025667678715e-07
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 8.9201209497146e-13,
  "early_nfev": 7616,
  "late_nfev": 23087,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like or compensated trial mode; no primordial spectrum imposed. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias10-eta_iso-k003.json

<!-- FILE: bias10-eta_iso-k003.json -->
```json
{
  "k": 0.003,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 100.0,
  "zstart": 1000000.0,
  "initial_mode": "eta_iso",
  "initial_amplitude": 1e-05,
  "psi_initial": 0.0,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 1.957993422274414e-05,
  "initial_radiation_delta_per_psi": -0.017922450992313113,
  "match_delta_m_per_psi": 0.0007032949290023173,
  "match_eta_delta_n_per_psi": -0.18678481354940205,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 2.169582305722289e-06,
  "max_matter_fraction": 4.030260661415028e-07,
  "max_phi": 1.3994072392941698e-09,
  "max_delta_x": 5.053337433094679e-09,
  "max_delta_x_over_x_floor": 1.1152875223231449e-08,
  "at_stop": {
    "k_over_aH": 13.303021872731664,
    "eta_number": 2.1993139688988262e-07,
    "eta_number_rest_frame": 2.1822148346164186e-07,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": 1.1864122740547977e-11,
    "delta_x_over_x": -1.9554492830890212e-09
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 1.6098630072749567e-10,
  "early_nfev": 15323,
  "late_nfev": 23087,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like or compensated trial mode; no primordial spectrum imposed. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias10-k.001.json

<!-- FILE: bias10-k.001.json -->
```json
{
  "k": 0.001,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 0.07333376759425343,
  "initial_radiation_delta_per_psi": -2.0000001645535135,
  "match_delta_m_per_psi": -5.1290800954625215,
  "match_eta_delta_n_per_psi": -3.2250664650472,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.00027125325788898117,
  "max_matter_fraction": 0.00034623150563074627,
  "max_phi": 9.857598828498858e-06,
  "max_delta_x": 6.24487487518958e-05,
  "max_delta_x_over_x_floor": 6.362122461239318e-05,
  "at_stop": {
    "k_over_aH": 4.434340624243888,
    "eta_number": -0.00027125325788898117,
    "eta_number_rest_frame": -0.0002446711600622243,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": -1.3179334701590117e-08,
    "delta_x_over_x": 2.9229871728200537e-05
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 1.7895149598071457e-09,
  "early_nfev": 54563,
  "late_nfev": 12827,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias10-k.002.json

<!-- FILE: bias10-k.002.json -->
```json
{
  "k": 0.002,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 0.07333376759425343,
  "initial_radiation_delta_per_psi": -2.000000658214054,
  "match_delta_m_per_psi": -14.194073958487161,
  "match_eta_delta_n_per_psi": -2.7570592647467835,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.0009853405389636885,
  "max_matter_fraction": 0.0012664736333545222,
  "max_phi": 9.543382823684435e-06,
  "max_delta_x": 7.233551719827297e-05,
  "max_delta_x_over_x_floor": 7.44915131159623e-05,
  "at_stop": {
    "k_over_aH": 8.868681248487777,
    "eta_number": -0.0009853405389636885,
    "eta_number_rest_frame": -0.0009562349866645981,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": -5.232667445047593e-08,
    "delta_x_over_x": 2.431630596506926e-05
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 5.93055919076295e-10,
  "early_nfev": 31868,
  "late_nfev": 12827,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias10-k.003.json

<!-- FILE: bias10-k.003.json -->
```json
{
  "k": 0.003,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 0.07333376759425343,
  "initial_radiation_delta_per_psi": -2.000001480981622,
  "match_delta_m_per_psi": -28.5015200162751,
  "match_eta_delta_n_per_psi": -2.751439045143794,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.0022852583320930116,
  "max_matter_fraction": 0.0027158248338918123,
  "max_phi": 9.200213041493877e-06,
  "max_delta_x": 8.13355655682631e-05,
  "max_delta_x_over_x_floor": 8.538892697577029e-05,
  "at_stop": {
    "k_over_aH": 13.303021872731664,
    "eta_number": -0.0022852583320930116,
    "eta_number_rest_frame": -0.0022529054262068445,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": -1.238773234366561e-07,
    "delta_x_over_x": 1.0046373090508922e-05
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 3.6223880505648794e-10,
  "early_nfev": 23771,
  "late_nfev": 12827,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias10-k.006.json

<!-- FILE: bias10-k.006.json -->
```json
{
  "k": 0.006,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 0.07333376759425343,
  "initial_radiation_delta_per_psi": -2.0000059239264862,
  "match_delta_m_per_psi": -96.29962436825966,
  "match_eta_delta_n_per_psi": -2.455603533688307,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.005279686382285053,
  "max_matter_fraction": 0.009551739937570873,
  "max_phi": 8.165270071880977e-06,
  "max_delta_x": 2.416668672538127e-05,
  "max_delta_x_over_x_floor": 2.9135296099428418e-05,
  "at_stop": {
    "k_over_aH": 26.606043745463328,
    "eta_number": -0.005279686382285053,
    "eta_number_rest_frame": -0.00526263835996935,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": -2.8773405094437326e-07,
    "delta_x_over_x": -5.592724216319457e-06
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 2.525702858613002e-11,
  "early_nfev": 15236,
  "late_nfev": 12827,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### bias10-k.01.json

<!-- FILE: bias10-k.01.json -->
```json
{
  "k": 0.01,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1e-08,
  "shape": 0.0,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 0.07333376759425343,
  "initial_radiation_delta_per_psi": -2.00001645535135,
  "match_delta_m_per_psi": -225.92986333943946,
  "match_eta_delta_n_per_psi": -2.1207304895082606,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.00662620935872804,
  "max_matter_fraction": 0.022609295212194087,
  "max_phi": 6.971323133314149e-06,
  "max_delta_x": 5.594083745589983e-06,
  "max_delta_x_over_x_floor": 6.831064703216804e-06,
  "at_stop": {
    "k_over_aH": 44.34340624243888,
    "eta_number": -0.00662620935872804,
    "eta_number_rest_frame": -0.006617522379172894,
    "eta_background_fraction": 5.443758442287916e-05,
    "eta_density_over_total": -3.608426454272608e-07,
    "delta_x_over_x": -2.4522928394308397e-06
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 3.362260556269527e-12,
  "early_nfev": 11621,
  "late_nfev": 12827,
  "background_today_x": 1.0238711134981502,
  "background_max_x": 1.2530893329619288,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### control-no-late-metric.json

<!-- FILE: control-no-late-metric.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 64,
  "rtol": 1e-07,
  "metric": false,
  "seed": 1e-08,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.0,
  "eta_number_threshold": 1.0,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.0001276663934515,
  "match_delta_m_per_psi": -1000.0911987884457,
  "match_eta_delta_n_per_psi": -1.0498416727301105,
  "z_at_stop": 2.124689152331134,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 1.0000000000000009,
  "max_matter_fraction": 0.020058869430154672,
  "max_phi": 0.0,
  "max_delta_x": 0.19707855245068948,
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "early_nfev": 8234,
  "late_nfev": 10742,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### direct_pole_bias_background_checks.py

<!-- FILE: direct_pole_bias_background_checks.py -->
```python
"""Background, seed dependence, and explicit-bias wall diagnostics.
Wall tension is calculated in the exactly degenerate unbiased potential as
an estimate; no network evolution or inflation initial spectrum assumed.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid,quad,cumulative_trapezoid
from scipy.optimize import brentq
from direct_pole_metric_pilot import background,potential,X,VE,M0,T0,MPL,H0,OM,OR,C,HBAR,MPC

BIAS=.0032765520306085277; FACTOR=10.; VREF=np.sqrt(2)*44.1e5; V0=17.834903110469433
BASE_RATIO=52.43035455726708; BASE_PROXY=.008746239106334854

def run_seed(seed):
    sol,loc=background(seed=seed,vev_factor=FACTOR,bias=BIAS)
    t=np.r_[np.linspace(-np.log(1e10+1),-np.log(1e5+1),6000)[:-1],np.linspace(-np.log(1e5+1),-np.log(101),5000)[:-1],np.linspace(-np.log(101),0,12000)]
    x,dx=sol.sol(t); H2,dh,G,er,ep,rm,rr=np.array([loc(a,b,c) for a,b,c in zip(t,x,dx)]).T; H=np.sqrt(H2); U=np.array([potential(a,bias=BIAS)[0] for a in x]);K=.5*(VE*FACTOR)**2*H2*dx**2; extra=U+K+er; p=K-U+ep
    r=np.exp(-t); half=np.flatnonzero(abs(x)**4>=.5); late=r<=101.000001;halo=r<=3.; mid=(r<=1e5+1.01)&(r>=101-1e-7); high=r>=1e5+1-1e-7
    A4=trapezoid(x[halo]**4/H[halo],t[halo]); weight=r**4*x**4/H; pair=(10/(10+(V0*r/VREF)**2))**2
    # High-z extrapolation of the representative NR pair proxy is NOT a
    # physical collision rate once the pair becomes relativistic.
    fac=BASE_PROXY/(BASE_RATIO*A4); bstd=H0*np.sqrt(OM*r**3+OR*r**4+1-OM-OR)
    R=extra[late]/r[late]**4; ref=R[0]+cumulative_trapezoid((extra[late]-3*p[late])/r[late]**4,t[late],initial=0.)
    return dict(seed=seed,z_half=float(r[half[0]]-1) if half.size else None,today_x=float(x[-1]),sign_at_half=float(np.sign(x[half[0]])) if half.size else None,
        max_extra_zle100=float(max(extra[late]/(3*MPL*MPL*H2[late]))),h0_ratio=float(H[-1]/H0),max_h_ratio_zle100=float(max(H[late]/bstd[late])),continuity_rel=float(max(abs(R-ref))/max(abs(R))),
        late_contact_proxy=float(fac*trapezoid(weight[late],t[late])),early_contact_proxy_100_1e5=float(fac*trapezoid(weight[mid],t[mid])),early_pair_proxy_100_1e5=float(fac*trapezoid((weight*pair)[mid],t[mid])),
        high_contact_extrapolation_1e5_1e10=float(fac*trapezoid(weight[high],t[high])),high_pair_NR_extrapolation_1e5_1e10=float(fac*trapezoid((weight*pair)[high],t[high])),
        max_abs_x_zge100=float(max(abs(x[~late]))),activation_integral=float(A4),scope='High redshift proxies extrapolate frozen NR scaling; do not interpret as collision bounds.')

def wall():
    ve=VE*FACTOR; mu=np.sqrt(X)/ve; grad=lambda x:potential(x,bias=BIAS)[1]; neg=brentq(grad,-1.1,-.7); barrier=brentq(grad,-.2,.2);delta=potential(neg,bias=BIAS)[0]-potential(1,bias=BIAS)[0]
    tension=ve*quad(lambda x:np.sqrt(max(0.,2*potential(x)[0])),-1,1,epsabs=1e-15)[0]
    h=delta/tension; H0r=h/H0
    zann=None if H0r<1 else brentq(lambda z:np.sqrt(OM*(1+z)**3+OR*(1+z)**4+1-OM-OR)-H0r,0,1e4)
    return dict(vev_over_Mpl=ve/MPL,mu_eV=mu,negative_local_minimum=neg,barrier=barrier,delta_vacuum_eV4=delta,unbiased_wall_tension_eV3=tension,pressure_balance_H_over_H0=H0r,pressure_balance_z_lcdm=zann,one_wall_per_hubble_density_fraction_at_balance=tension/(3*MPL*MPL*h),one_wall_per_hubble_density_fraction_today=tension/(3*MPL*MPL*H0),bias_tadpole_eV3=2*BIAS*X/ve,bias_mass_squared_eV2=2*BIAS*X/ve**2,scope='Pressure balance H~delta V/tension; one wall per Hubble volume dimensional estimate, not a network or a prediction that walls form.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();r={'wall':wall(),'seeds':[]};a.output.write_text(json.dumps(r,indent=2)+'\n')
    for seed in (1e-8,0.,-1e-5,1e-5,-1e-4,1e-4,-1e-3,1e-3,-.01,.01,-.1,.1):
        d=run_seed(seed);r['seeds'].append(d);a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print(seed,d['z_half'],d['today_x'],d['early_contact_proxy_100_1e5'],d['max_extra_zle100'],flush=True)
    print(r['wall'])
```
<!-- END FILE -->

### direct_pole_bias_bao.py

<!-- FILE: direct_pole_bias_bao.py -->
```python
"""DESI DR2 compressed BAO, fixed H0 and Or, common ruler profiled.
Public input pinned at CobayaSampler/bao_data bb0c1c9009dc76d1391300e169e8df38fd1096db.
Covariance block entries transcribed exactly from published data file.
No CMB, likelihood of recombination, SN or H0 parameter refit.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize_scalar
import direct_pole_metric_pilot as m

BIAS=.0032765520306085277; OM0=.315
Z=np.array([.295,.51,.51,.706,.706,.934,.934,1.321,1.321,1.484,1.484,2.33,2.33]);KIND='V M H M H M H M H M H H M'.split()
DATA=np.array([7.94167639,13.58758434,21.86294686,17.35069094,19.45534918,21.57563956,17.64149464,27.60085612,14.17602155,30.51190063,12.81699964,8.631545674846294,38.988973961958784])
COV=np.diag([5.78998687e-03,2.83473742e-02,1.83928040e-01,3.23752442e-02,1.11469198e-01,2.61732816e-02,4.04183878e-02,1.05336516e-01,5.04233092e-02,5.83020277e-01,2.68336193e-01,1.02136194e-02,2.82685779e-01])
for i,c in zip((1,3,5,7,9,11),(-3.26062007e-02,-2.37445646e-02,-1.12938006e-02,-2.90308418e-02,-1.95215562e-01,-2.31395216e-02)):COV[i,i+1]=COV[i+1,i]=c
CHO=cho_factor(COV); C0=m.C*m.HBAR/m.MPC

def profile(d):
    w=cho_solve(CHO,d);normal=d@w;scale=DATA@w/normal;res=d*scale-DATA
    return dict(ruler_Mpc=float(1/scale),ruler_sigma_local=float(1/(scale**2*np.sqrt(normal))),chi2=float(res@cho_solve(CHO,res)),predicted=(d*scale).tolist())

def evaluate(om,model,full=False):
    m.OM=om
    try:
        if model=='lcdm':
            def hh(t):return m.H0*np.sqrt(om*np.exp(-3*t)+m.OR*np.exp(-4*t)+1-om-m.OR)
        else:
            sol,loc=m.background(vev_factor=10. if model=='biased' else 1.,bias=BIAS if model=='biased' else 0.)
            def hh(t):x,dx=sol.sol(t);return np.sqrt(loc(t,x,dx)[0])
        d={}
        for z in np.unique(Z):
            t=-np.log1p(z);dh=C0/hh(t);dm=quad(lambda s:C0*np.exp(-s)/hh(s),t,0,epsabs=2e-5,epsrel=2e-8,limit=150)[0];d[z]={'H':dh,'M':dm,'V':np.cbrt(z*dm*dm*dh)}
        ds=np.array([d[z][k] for z,k in zip(Z,KIND)]);r=dict(model=model,omega_m=om,distance_Mpc=ds.tolist(),**profile(ds))
        if full:
            def rs(z):return quad(lambda s:C0*np.exp(-s)/(hh(s)*np.sqrt(3*(1+3*.02237/(4*2.4728e-5)*np.exp(s)))),-np.log(1e14),-np.log1p(z),epsabs=2e-6,epsrel=2e-8,limit=200)[0]
            drag=rs(1060);star=rs(1090);dm=quad(lambda s:C0*np.exp(-s)/hh(s),-np.log1p(1090),0,epsabs=1e-5,epsrel=2e-8,limit=200)[0];res=ds/drag-DATA
            r.update(approx_fixed_drag_ruler=drag,approx_fixed_drag_chi2=float(res@cho_solve(CHO,res)),D_M_z1090=dm,r_s_z1090=star,theta_s_z1090=star/dm,H0_ratio=hh(0)/m.H0)
        return r
    finally:m.OM=OM0

def run(output):
    result={'source':'https://github.com/CobayaSampler/bao_data/tree/bb0c1c9009dc76d1391300e169e8df38fd1096db/desi_bao_dr2','source_cov_sha256':'252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509','source_mean_sha256':'9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585','results':{},'scope':'DESI DR2 only, fixed H0=67.4, fixed Or, physical trigger parameters held fixed during Omega_m profile; fixed-redshift sound horizon diagnostic not a fitted ruler.'}
    for model in ('lcdm','unbiased','biased'):
        r=evaluate(OM0,model,True);result['results'][model]={'fixed_omega_m':r};output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(model,'fixed',r['chi2'],r['ruler_Mpc'],r['approx_fixed_drag_chi2'],r['theta_s_z1090'],flush=True)
    for model in ('lcdm','biased'):
        opt=minimize_scalar(lambda om:evaluate(om,model)['chi2'],bounds=(.23,.4),method='bounded',options={'xatol':1e-5});assert opt.success
        r=evaluate(float(opt.x),model,True);result['results'][model]['profiled_omega_m']=r;output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(model,'profile',r['chi2'],r['omega_m'],r['ruler_Mpc'],r['approx_fixed_drag_chi2'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.output)
```
<!-- END FILE -->

### direct_pole_bias_early_bound.py

<!-- FILE: direct_pole_bias_early_bound.py -->
```python
"""Conditional early elastic bound and finite-field matching rescaling.
Uses the earlier rearrangement bound with a time-dependent *upper* vacuum
matrix element (all s,t,u exchange denominators bounded). Positive mass,
occupation and total number caps are assumptions. Medium shifts, coherent
fields, portal production, number-changing reactions not covered.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from direct_pole_metric_pilot import background,VE,X,MPL,H0,OR
from psd_reference import get_psd

M=86.; HBARC=1.973269804e-5; CGS=HBARC**3/1.7826619216278976e-33; N0=.12*1.8783416e-29*CGS/M
DET=.0032008304597957726; MS2=4*M*M-DET; BASE_KAPPA=7.92823162096704e-31; LS=3.4170831596658935e-6
BIAS=.0032765520306085277; FACTOR=10.; SIG_RATIO=.22236922325533096/.22089528043464016

def matching():
    g2=(BASE_KAPPA*VE)**2*np.sqrt(SIG_RATIO);kappa=np.sqrt(g2)/(VE*FACTOR); cap=.1*MPL
    la=9*(kappa*cap)**4/(4*LS*M**4)
    return dict(g2_eV2=g2,kappa=kappa,lambda_S=LS,field_cap_over_Mpl=.1,lambda_A=la,lambda_threshold_bound=g2*(1/DET+2/MS2),original_finite_loop_ratio=.7514722644010279*np.sqrt(SIG_RATIO),constant_A4_finite_loop_ratio=.671425531844289*np.sqrt(SIG_RATIO),scope='Finite-field tree bound and finite scale-local thresholds only. Global quartic witness unchanged under v rescaling up to the slight coupling retune.')

def integrate(seed,mesh=240,n=12000):
    mat=matching(); F={240:21525656357.32168,320:21441710365.93471}[mesh]; npref=2*N0; a=npref/M; pstar=(6*np.pi*np.pi*npref/F)**(1/3);b=F*pstar*pstar/(4*np.pi*np.pi);h=.5*H0*np.sqrt(OR); pref=4*(1+2*F)/(32*np.pi*N0*h)
    sol,_=background(seed=seed,vev_factor=FACTOR,bias=BIAS);t=np.linspace(-np.log(1e14),-np.log(101),n);r=np.exp(-t);x=sol.sol(t)[0];envelope=np.minimum(a*a*r,b*b/r);amp=mat['lambda_A']+mat['lambda_threshold_bound']*x*x
    # Beyond the starting history only a constant |x|<=|seed| envelope
    # is assumed. Its analytic tail uses the relativistic branch.
    tail_amp=mat['lambda_A']+mat['lambda_threshold_bound']*seed*seed; tail=tail_amp**2*b*b/r[0]; upper=pref*(trapezoid(amp*amp*envelope,t)+tail)
    contact=pref*mat['lambda_A']**2*(trapezoid(envelope,t)+b*b/r[0]); constant_all=pref*mat['lambda_A']**2*2*a*b
    # Positive mode gain bound, I <= F p^2/(sqrt(p^2+M^2)+M).
    p=pstar*r; I=F*p*p/(np.hypot(p,M)+M); gain=amp*amp*(1+2*F)*I*I/(64*np.pi**3*M*h*r*r)
    # Large-r tail using I<=F*p; for this conservative envelope I^2/r^2
    # tends to a constant, so integrate only to prior finite production rstart.
    rstart={240:11159824101379422000,320:11159719732901884000}[mesh]; gain_tail=tail_amp**2*(1+2*F)*(F*pstar)**2/(64*np.pi**3*M*h)*np.log(rstart/r[0])
    return dict(seed=seed,mesh=mesh,n=n,constant_only_all_r_L1_bound=constant_all,constant_only_rge101_L1_bound=contact,combined_amplitude_rge101_L1_bound=upper,positive_mode_gain_to_prior_production_start=float(trapezoid(gain,t)+gain_tail),occupation_cap=F,rearrangement_transition_r=b/a,lambda_at_z100=float(amp[-1]),peak_abs_x_history=float(max(abs(x))),scope='Vacuum tree amplitude, positive mass >=86 eV, number <=2 N0 r^3 and f<=cap, H>=0.5 H0 sqrt(Or) r^2; frozen seed bound beyond history. No in-medium pole or primordial S/eta/chi, coherence or number-changing processes.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();r={'matching':matching(),'rows':[]}
    for mesh in (240,320):
        for seed in (0.,1e-8,1e-7,1e-6,3e-6,1e-5,1e-4):
            d=integrate(seed,mesh);r['rows'].append(d);a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n'); print(mesh,seed,d['combined_amplitude_rge101_L1_bound'],d['positive_mode_gain_to_prior_production_start'],flush=True)
```
<!-- END FILE -->

### direct_pole_bias_elastic.py

<!-- FILE: direct_pole_bias_elastic.py -->
```python
"""Independent NR conservative finite-volume Bose operator on archived PSDs.
Same nodal graph, interpolated reverse occupation, and normalization as the
published direct-pole benchmark, but analytic NR energy and separate code.
Only the post-production homogeneous collision integral at z<=100 is tested.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,cumulative_trapezoid,trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre
from numba import njit
from direct_pole_metric_pilot import background
from psd_reference import get_psd

BIAS=.0032765520306085277
VREF=np.sqrt(2)*44.1e5

def geometry(e):
    vol=np.diff(e**3)/(6*np.pi**2); K=.6*np.diff(e**5)/np.diff(e**3); p=np.sqrt(K); dK=np.diff(e**2)
    ii=[];jj=[];kk=[];ll=[];ww=[];rr=[];ss=[];los=[];his=[]
    for i in range(len(K)):
        for j in range(i,len(K)):
            s=K[i]+K[j]; lo=max(K[0],s-K[-1]); hi=min(s/2,K[-1])
            if hi<=lo:continue
            L=np.maximum(K[:-1],lo);H=np.minimum(K[1:],hi);width=np.maximum(H-L,0);t=np.clip(((L+H)/2-K[:-1])/np.diff(K),0,1)
            hats=np.r_[width*(1-t),0]+np.r_[0,width*t]; k=np.flatnonzero(hats>0);k4=s-K[k];ok=(k4>=K[0])&(k4<=K[-1]);k=k[ok];k4=k4[ok]
            if not len(k):continue
            l=np.clip(np.searchsorted(K,k4)-1,0,len(K)-2);w=np.clip((k4-K[l])/(K[l+1]-K[l]),0,1);p4=np.sqrt(k4)
            low=np.maximum(abs(p[i]-p[j]),abs(p[k]-p4));high=np.minimum(p[i]+p[j],p[k]+p4);rate=2*(1 if i==j else 2)*dK[i]*dK[j]*hats[k]*np.maximum(high-low,0)
            ok=rate>0;n=ok.sum();ii.extend([i]*n);jj.extend([j]*n);kk.extend(k[ok]);ll.extend(l[ok]);ww.extend(w[ok]);rr.extend(rate[ok]);ss.extend([2*s]*n);los.extend(low[ok]);his.extend(high[ok])
    idx=np.array([ii,jj,kk,ll,np.array(ll)+1],np.int64).T; w=np.asarray(ww); coef=np.stack([np.ones_like(w)]*3+[1-w,w],axis=1)
    return vol,K,idx,coef,np.asarray(rr),np.asarray(ss),np.asarray(los),np.asarray(his)

def pole(s,l,h,b):
    d=b+s;root=np.sqrt(d);width=h-l;dh=b+np.maximum(s-h*h,0);dl=b+np.maximum(s-l*l,0)
    lg=np.log1p(width/(root+l))+np.log1p(width*(root+h)/dh)
    return b*b/(2*d)*((d+h*l)/(dh*dl)+lg/(2*root*width))

@njit(cache=True)
def flux_jac(f,idx,coef,rate,vol,jac):
    n=len(f);ff=np.maximum(f,1e-90);lf=np.log(ff);r=np.log1p(1/ff);out=np.zeros(n);J=np.zeros((n,n));ds=np.empty(5);sg=np.array([-1.,-1.,1.,1.,1.])
    for h in range(len(rate)):
        ids=idx[h];c=coef[h];pr=0.
        for a in range(5):pr+=c[a]*lf[ids[a]]
        rin=r[ids[0]]+r[ids[1]];rout=r[ids[2]]+c[3]*r[ids[3]]+c[4]*r[ids[4]];diff=rout-rin
        F=rate[h]*np.sign(diff)*np.exp(pr+max(rin,rout))*(-np.expm1(-abs(diff)))
        for a in range(5):out[ids[a]]+=sg[a]*c[a]*F
        if jac:
            fw=rate[h]*np.exp(pr+rout);bw=rate[h]*np.exp(pr+rin)
            for a in range(5):
                g=ff[ids[a]]
                if f[ids[a]]<=1e-90:ds[a]=0.
                elif a<2:ds[a]=c[a]*(fw/g-bw/(1+g)) if g<1e-3 else c[a]*(F/g+bw/(g*(1+g)))
                else:ds[a]=c[a]*(fw/(1+g)-bw/g) if g<1e-3 else c[a]*(F/g-fw/(g*(1+g)))
            for a in range(5):
                for b in range(5):J[ids[a],ids[b]]+=sg[a]*c[a]*ds[b]
    return out/vol,J/vol[:,None]

def exposure(factor,bias):
    sol,loc=background(vev_factor=factor,bias=bias); t=np.linspace(-np.log(101),0,20001);x,dx=sol.sol(t); H=np.sqrt([loc(a,b,c)[0] for a,b,c in zip(t,x,dx)])
    h=t>=-np.log(3);A4=trapezoid(x[h]**4/H[h],t[h]);cumulative=cumulative_trapezoid(np.exp(-4*t)*x**4/H,t,initial=0);C=cumulative[-1];u=cumulative/C;unique=np.r_[True,np.diff(u)>0]
    return C/A4,PchipInterpolator(u[unique],-t[unique]),A4

def solve(mesh,factor,bias,slices,nodes,rtol):
    e,f,tau,v0,ref=get_psd(mesh);vol,K,idx,coef,rate,S,lo,hi=geometry(e);N=vol@f;E=(vol*K)@f; assert np.max(abs((K[idx]*coef)@np.array([-1,-1,1,1,1]))/np.maximum(S,1e-100))<2e-15
    base,_,A0=exposure(1,0);ratio,logr,A=exposure(factor,bias);tau_new=tau*ratio/base;d0=10*(VREF/v0)**2;scale=tau_new/(128*np.pi**4*N);normal=np.maximum(f,1.);y=np.maximum(f,1e-38)/normal
    moments=[vol*K*K]+[vol*np.maximum(e[1:]**3-np.maximum(e[:-1],c)**3,0)/np.diff(e**3) for c in (2,6,20)];M=np.array(moments);m0=M@f
    p,w=roots_legendre(nodes);res=[];neg=0.; maxn=0.;maxe=0.;de=0.
    for j in range(slices):
        a,b=j/slices,(j+1)/slices; avg=np.zeros_like(rate)
        for u,wt in zip((a+b+(b-a)*p)/2,w):avg+=wt/2*pole(S,lo,hi,d0*np.exp(-2*logr(u)))
        ratej=rate*avg
        if j==0:
            eq=1/np.expm1(.3+80*K/K[-1]);r0,_=flux_jac(eq,idx,coef,ratej,vol,False);de=float(vol@abs(r0)/np.sum(ratej*eq[idx[:,0]]*eq[idx[:,1]]))
        def rhs(t,y):return flux_jac(y*normal,idx,coef,ratej,vol,False)[0]*scale/normal
        def jac(t,y):return flux_jac(y*normal,idx,coef,ratej,vol,True)[1]*scale*normal[None,:]/normal[:,None]
        s=solve_ivp(rhs,(a,b),y,method='BDF',jac=jac,rtol=rtol,atol=1e-13,max_step=(b-a)/8);assert s.success,s.message
        ff=s.y*normal[:,None];neg=max(neg,float(max(vol@np.maximum(-ff,0)/N)));maxn=max(maxn,float(max(abs(vol@ff/N-1))));maxe=max(maxe,float(max(abs((vol*K)@ff/E-1))));assert neg<1e-9
        y=s.y[:,-1];res.append(dict(slice=j+1,nfev=s.nfev,moment_changes=(M@(y*normal)/m0-1).tolist()));print(mesh,factor,j+1,res[-1],flush=True)
    return dict(mesh=mesh,n=len(f),events=len(rate),vev_factor=factor,bias=bias,slices=slices,nodes=nodes,rtol=rtol,reference_tail=ref,tau=tau_new,tau_ratio=ratio/base,A4_ratio=A/A0,sigma0=.22089528043464016*A0/A,relative_moments=res[-1]['moment_changes'],max_number_drift=maxn,max_energy_drift=maxe,max_negative_number_fraction=neg,detailed_balance_residual=de,steps=res,scope='Conditional post-production homogeneous NR Bose 2-to-2 scattering on reconstructed 220-bin PSD; same graph convention as prior archive. No early collision or in-medium calculation.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mesh',type=int,choices=(240,320),default=240);p.add_argument('--baseline',action='store_true');p.add_argument('--slices',type=int,default=8);p.add_argument('--nodes',type=int,default=16);p.add_argument('--rtol',type=float,default=1e-6);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=solve(a.mesh,1. if a.baseline else 10.,0. if a.baseline else BIAS,a.slices,a.nodes,a.rtol);a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print(json.dumps(r),flush=True)
```
<!-- END FILE -->

### direct_pole_bias_local_modes.py

<!-- FILE: direct_pole_bias_local_modes.py -->
```python
"""Frozen local collisionless dispersion scan: necessary sensitivity check,
not a time-integrated phase evolution or a proof of linear stability.
The static infrared response is regularized by finite positive gamma.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import brentq
from scipy.special import roots_legendre
from direct_pole_metric_pilot import background,potential,VE,M0,T0,COEF,C,HBAR,MPC

BIAS=.0032765520306085277

def response(y):
    y=np.asarray(y);out=np.empty_like(y);low=y<.003;u=y[low]**2;out[low]=u*(1/3+u*(-1/5+u/7));out[~low]=1-np.arctan(y[~low])/y[~low];return out

def scan(factor,bias,order=20,ns=400,nk=110):
    qedge=np.r_[0.,np.geomspace(1e-13,.1,15),.3,1.,2.,4.,8.,16.,32.,56.];u,w=roots_legendre(order);r0=np.diff(qedge)/2
    q=((qedge[1:]+qedge[:-1])[:,None]/2+r0[:,None]*u).ravel();dq=(r0[:,None]*w).ravel();f=1/np.expm1(q);fd=f*(1+f)
    sol,loc=background(vev_factor=factor,bias=bias);ve=VE*factor;g=M0/ve
    ts=np.linspace(-np.log(101),0,ns);ks=np.geomspace(1e-5,10,nk);x,dx=sol.sol(ts);rs=np.exp(-ts);grid=np.zeros((ns,nk));unc=np.zeros_like(grid);wave=np.zeros_like(grid);curv=np.zeros(ns)
    for i,(t,r,xx,dxx) in enumerate(zip(ts,rs,x,dx)):
        H2=loc(t,xx,dxx)[0];H=np.sqrt(H2);T=T0*r;m=M0*xx;e=np.hypot(q,m/T);v=q/e
        D=COEF*g*g*T*T*np.sum(dq*q**4*f/e**3);S=COEF*g*g*m*m*dq*v*fd/H2;V2=potential(xx,bias=bias)[2]/ve**2;curv[i]=(V2+D)/H2
        wave[i]=ks*r*C*HBAR/(MPC*H);w2=wave[i]**2
        unc[i]=np.maximum((np.sqrt(np.maximum(9-4*(w2+curv[i]),0))-3)/2,0)
        for j,wh in enumerate(wave[i]):
            def disp(gam):return gam*(gam+3)+wh**2+curv[i]-S@response(wh*v/gam)
            if disp(1e-5)<0:
                upper=max(1,np.sqrt(max(0.,S.sum()-wh**2-curv[i]))+1);grid[i,j]=brentq(disp,1e-5,upper,xtol=1e-9)
    idx=np.flatnonzero(abs(x)**4>=.5)[0];intg=trapezoid(grid,ts,axis=0);late=trapezoid(grid[idx:],ts[idx:],axis=0);ungu=trapezoid(unc,ts,axis=0);late_un=trapezoid(unc[idx:],ts[idx:],axis=0)
    valid=(wave>=10)&(grid>=5)&(grid>=3*np.abs(dx/np.maximum(abs(x),1e-30))[:,None]);fast=trapezoid(np.where(valid,grid,0.),ts,axis=0)
    return dict(vev_factor=factor,bias=bias,order=order,ns=ns,nk=nk,z_half_sample=float(rs[idx]-1),max_gamma_over_H=float(grid.max()),k_peak_integral=float(ks[np.argmax(intg)]),max_frozen_positive_integral_all=float(intg.max()),max_frozen_positive_integral_post_half=float(late.max()),max_unforced_post_half=float(late_un.max()),max_eta_increment_post_half=float(max(late-late_un)),max_fast_subhorizon_integral=float(max(fast)),min_curvature_after_half=float(min(curv[idx:])),k=ks.tolist(),all_integral=intg.tolist(),post_half_integral=late.tolist(),post_half_unforced=late_un.tolist(),fast_integral=fast.tolist(),scope='Positive frozen roots can overestimate or underestimate true growth across time-dependent phases. No metric feedback or actual time evolution in this diagnostic.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();rows=[]
    for f,b,n,k in [(1.,0.,480,100),(10.,BIAS,480,100),(10.,BIAS,720,160)]:
        row=scan(f,b,20 if n==480 else 28,n,k);rows.append(row);a.output.write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n');print({k:v for k,v in row.items() if not isinstance(v,list)},flush=True)
```
<!-- END FILE -->

### direct_pole_bias_parameter_window.py

<!-- FILE: direct_pole_bias_parameter_window.py -->
```python
"""Small prospective screen at fixed nominal z_on=35 and DeltaNeff=.1.
Vary X, with m_eta^2 proportional to X as in original thermal matching;
CW ratio scales with X, base vev scales with sqrt(X) to keep mu fixed.
The bias is retuned to half activation. No likelihood, full collision
solve, or perturbation solve is performed for these new screened rows.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import brentq
import direct_pole_metric_pilot as m

BASE=(m.X,m.VE,m.M0,m.CW); BASE_A4=4.95835323e32; BASE_RATIO=52.43035455726708; BASE_PROXY=.008746239106334854; BASE_LOOP=.7514722644010279; BASE_SIG=.22089528043464016

def screen(xfactor,vev_factor,target_z):
    x0,v0,m0,cw=BASE; m.X=x0*xfactor; m.VE=v0*np.sqrt(xfactor); m.M0=m0*np.sqrt(xfactor); m.CW=cw*xfactor
    try:
        target=-np.log1p(target_z);hit=.5**.25
        def mismatch(ls):
            sol,_=m.background(vev_factor=vev_factor,bias=np.exp(ls),end=target,stop_x=hit)
            return (target-sol.t[-1])*max(abs(sol.y[1,-1]),.001) if sol.t_events[0].size else abs(sol.y[0,-1])-hit
        lo,hi=np.log(1e-14),np.log(.8)
        if mismatch(lo)*mismatch(hi)>0:return dict(xfactor=xfactor,vev_factor=vev_factor,target_z=target_z,tuned=False)
        bias=float(np.exp(brentq(mismatch,lo,hi,xtol=1e-6)));sol,loc=m.background(vev_factor=vev_factor,bias=bias)
        t=np.linspace(-np.log(101),0,7001);x,dx=sol.sol(t);H2,dh,G,er,ep,rm,rr=np.array([loc(a,b,c) for a,b,c in zip(t,x,dx)]).T;H=np.sqrt(H2);r=np.exp(-t);T=m.T0*r;ve=m.VE*vev_factor
        U,V,V2=np.array([m.potential(a,bias=bias) for a in x]).T;q,dq,f=m.radial(6);e=np.hypot(q[None,:],(m.M0*x/T)[:,None]);direct=m.COEF*m.M0**2*T*T*(1/e**3@(dq*q**4*f));curv=(V2+direct)/(ve**2*H2)
        K=.5*ve**2*H2*dx**2;frac=(U+K+er)/(3*m.MPL**2*H2);halo=r<=3;A4=trapezoid(x[halo]**4/H[halo],t[halo]);cosmic=trapezoid(r**4*x**4/H,t);ratio=cosmic/A4;sig=BASE_SIG*BASE_A4/A4;loop=BASE_LOOP*np.sqrt(sig/BASE_SIG)/xfactor;proxy=BASE_PROXY*ratio/BASE_RATIO
        ix=np.flatnonzero(abs(x)**4>=.5)[0];root=trapezoid(np.maximum((np.sqrt(np.maximum(9-4*curv[ix:],0))-3)/2,0),t[ix:]);
        return dict(xfactor=xfactor,vev_factor=vev_factor,target_z=target_z,bias=bias,tuned=True,max_extra=float(max(frac)),max_loop_over_X=float(loop),cosmic_contact_proxy=float(proxy),sigma0=float(sig),positive_local_root_integral=float(root),min_curvature_after_half=float(min(curv[ix:])),max_abs_x=float(max(abs(x))),vev_over_Mpl=float(ve/m.MPL),passes=bool(max(frac)<.02 and loop<=1 and proxy<=.01),X=m.X,m_eta=m.M0,CW=m.CW)
    finally:m.X,m.VE,m.M0,m.CW=BASE

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--group',type=int,default=0);a=p.parse_args();cases=[(.8,10.,3.8),(.85,10.,3.8),(.9,10.,3.8),(.8,15.,3.8),(.85,15.,3.8),(.9,15.,3.8),(.8,10.,4.2),(.85,10.,4.2),(.9,10.,4.2),(.8,15.,4.2),(.85,15.,4.2),(.9,15.,4.2)];rows=[]
    for i,c in enumerate(cases):
        if i//6!=a.group:continue
        d=screen(*c);rows.append(d);a.output.write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n');print(d,flush=True)
```
<!-- END FILE -->

### direct_pole_metric_acoustic.py

<!-- FILE: direct_pole_metric_acoustic.py -->
```python
"""Independent early-transfer sensitivity: sudden photon/baryon decoupling.
Standard neutrinos and post-decoupling photons free-stream. Before z_dec,
photons and baryons share a tight-coupling velocity. No recombination width,
photon polarization, baryon pressure after decoupling, or reionization.
The late eta and scalar equations use the same stress convention as the
metric pilot, with phi evolved directly and the density constraint audited.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre
from direct_pole_metric_pilot import background,potential,radial,MPL,HBAR,C,MPC,VE,M0,T0,COEF,SEED,RHO0,OM,OR


def run(k=.03,order=5,angles=96,rtol=1e-7,zdec=1090.,zmatch=12.,amplitude=1e-5,baryon_fraction=.049/.315,nu_fraction=.409,vev_factor=1.,bias=0.,control="full"):
    ve=VE*vev_factor; assert control in ("full","lcdm","massless_eta")
    if control=="full":
        sol,loc=background(vev_factor=vev_factor,bias=bias); field=sol.sol
        m0=M0; coef=COEF
        pot=lambda x:potential(x,bias=bias)
    else:
        m0=0.; coef=COEF if control=="massless_eta" else 0.
        field=lambda t:np.zeros(2); pot=lambda x:(0.,0.,0.)
        def loc(t,x,dx):
            r=np.exp(-t); er=coef*(T0*r)**4*np.pi**4/15; ep=er/3
            rm=RHO0*OM*r**3; rr=RHO0*OR*r**4
            H2=(rm+rr+er+RHO0*(1-OM-OR))/(3*MPL*MPL)
            return H2,-(rm+4*(rr+er)/3)/(2*MPL*MPL*H2),0.,er,ep,rm,rr
    q,dq,f=radial(order); fd=f*(1+f); q2=coef*dq*q*q*fd; nweight=dq*q*q*fd/np.sum(dq*q*q*f)
    mu,wm=roots_legendre(angles); mu=(mu+1)/2; wm=wm/2; m1=wm*mu; P2=(3*mu*mu-1)/2; m2=wm*P2
    t0=-np.log1p(1e6); td=-np.log1p(zdec); tm=-np.log1p(zmatch)
    def env(t):
        x,dx=field(t); H2,dh,G,er,ep,rm,rr=loc(t,x,dx); T=T0*np.exp(-t); rho=3*MPL*MPL*H2; w=k*np.exp(-t)*C*HBAR/MPC/np.sqrt(H2)
        rb=baryon_fraction*rm; rc=rm-rb; rn=nu_fraction*rr; rg=rr-rn
        return x,dx,H2,dh,G,er,ep,rc,rb,rn,rg,T,rho,w
    def early_rhs(t,s,dec):
        dc,vc,db,vb,dg,vg,u,du=s[:8]; R,I,Rn,In,Rg,Ig=s[8:].reshape(6,angles)
        x,dx,H2,dh,G,er,ep,rc,rb,rn,rg,T,rho,w=env(t); U,V,V2=pot(x)
        b0=R@wm; b1=-I@m1; b2=-R@m2; bn0=Rn@wm; bn1=-In@m1; bn2=-Rn@m2
        bg0=Rg@wm if dec else dg/4; bg1=-Ig@m1 if dec else vg/3; bg2=-Rg@m2 if dec else 0.
        slip=12*(er*b2+rn*bn2+rg*bg2)/(rho*w*w)
        mom=(rc*vc+rb*vb+4*(er*b1+rn*bn1+rg*bg1)+ve**2*H2*w*dx*u)/rho
        d0=(rc*dc+rb*db+4*(er*b0+rn*bn0+rg*bg0)+ve**2*H2*dx*du+V*u)/rho
        psi=(-1.5/w**2*(d0+3*mom/w)-slip)/(1-1.5*ve**2*H2*dx*dx/rho/w**2); phid=-psi+1.5*mom/w
        dR=phid+w*mu*I; dI=-w*mu*(R+psi); dRn=phid+w*mu*In; dIn=-w*mu*(Rn+psi)
        if dec:
            dRg=phid+w*mu*Ig; dIg=-w*mu*(Rg+psi); dvg=0.; dvb=-vb+w*psi; ddg=0.
        else:
            dRg=np.zeros_like(Rg); dIg=np.zeros_like(Ig); inertia=3*rb/(4*rg)
            dvg=(-inertia*vg+w*dg/4)/(1+inertia)+w*psi; dvb=dvg; ddg=-4*w*vg/3+4*phid
        slipd=-2*slip-12*(er*(dR@m2)+rn*(dRn@m2)+rg*(dRg@m2))/(rho*w*w)
        ddu=-(3+dh)*du-(w*w+(V2+m0*m0*T*T/12)/(ve**2*H2))*u-2*m0*m0*x*T*T/12*b0/(ve**2*H2)+dx*(4*phid-slipd)-2*G/(ve**2*H2)*psi
        return np.r_[-w*vc+3*phid,-vc+w*psi,-w*vb+3*phid,dvb,ddg,dvg,du,ddu,dR,dI,dRn,dIn,dRg,dIg]
    x,dx,H2,dh,G,er,ep,rc,rb,rn,rg,T,rho,w=env(t0)
    v=w/2; b2=w*w/30; R=-5*b2*P2; I=-v*mu; slip=12*(er+rn)*b2/(rho*w*w); phi=1+slip
    mom=((rc+rb)*v+4*(er+rn+rg)*v/3)/rho; dr=(-2*w*w*phi/3-3*mom/w)*rho/(.75*(rc+rb)+rg+rn+er)
    initial=np.r_[.75*dr,v,.75*dr,v,dr,v,0.,0.,R+dr/4,I,R+dr/4,I,np.zeros(2*angles)]
    early=solve_ivp(lambda t,s:early_rhs(t,s,False),(t0,td),initial,rtol=rtol,atol=rtol*.005,method='DOP853',max_step=.025)
    assert early.success; jump=early.y[:,-1].copy(); jump[-2*angles:-angles]=jump[4]/4; jump[-angles:]=-jump[5]*mu
    # Stress (monopole/dipole and vanishing quadrupole) is continuous by construction.
    e2=solve_ivp(lambda t,s:early_rhs(t,s,True),(td,tm),jump,rtol=rtol,atol=rtol*.005,method='DOP853',max_step=.025)
    assert e2.success; si=e2.y[:,-1]; R,I,Rn,In,Rg,Ig=si[8:].reshape(6,angles)
    start=np.r_[si[:8],(q[:,None]*R).ravel(),(q[:,None]*I).ravel(),Rn,In,Rg,Ig,0.]
    def late_rhs(t,s,deriv=True):
        dc,vc,db,vb,dg,vg,u,du=s[:8]; R,I=s[8:8+2*q.size*angles].reshape(2,q.size,angles); Rn,In,Rg,Ig=s[-1-4*angles:-1].reshape(4,angles); phi=s[-1]
        x,dx,H2,dh,G,er,ep,rc,rb,rn,rg,T,rho,w=env(t); U,V,V2=pot(x)
        g0=R@wm; g1=-I@m1; g2=-R@m2; bn0=Rn@wm; bn1=-In@m1; bn2=-Rn@m2; bg0=Rg@wm; bg1=-Ig@m1; bg2=-Rg@m2
        y=m0*x/T; e=np.hypot(q,y); velo=q/e; de=m0*m0*x/(T*T*e)
        direct=coef*T*T*m0*m0*np.sum(dq*q**4*f/e**3); J=coef*T**4*np.sum(dq*q*q*f*de)
        drho=T**4*np.sum(q2*e*g0)+J*u; eta_mom=T**4*np.sum(q2*q*g1)
        slip=(3*T**4*np.sum(q2*q*q/e*g2)+12*(rn*bn2+rg*bg2))/(rho*w*w)
        mom=(rc*vc+rb*vb+4*(rn*bn1+rg*bg1)+eta_mom+ve**2*H2*w*dx*u)/rho
        d0=(rc*dc+rb*db+4*(rn*bn0+rg*bg0)+drho+ve**2*H2*dx*du+V*u)/rho
        phi_c=(-1.5/w**2*(d0+3*mom/w)-slip)/(1-1.5*ve**2*H2*dx*dx/rho/w**2)+slip
        psi=phi-slip; phid=-psi+1.5*mom/w
        dR=q[:,None]*phid+w*mu[None,:]*velo[:,None]*I; dI=-w*mu[None,:]*(velo[:,None]*R+(e*psi+de*u)[:,None])
        dRn=phid+w*mu*In; dIn=-w*mu*(Rn+psi); dRg=phid+w*mu*Ig; dIg=-w*mu*(Rg+psi)
        edot=(m0/T)**2*x*(x+dx)/e
        slipd=-2*slip+3*T**4/(rho*w*w)*np.sum(q2*q*q/e*(-dR@m2-g2*edot/e))-12*(rn*(dRn@m2)+rg*(dRg@m2))/(rho*w*w)
        ddu=-(3+dh)*du-(w*w+(V2+direct)/(ve**2*H2))*u-T**4*np.sum(q2*de*g0)/(ve**2*H2)+dx*(4*phid-slipd)-2*(V+J)/(ve**2*H2)*psi
        dn=nweight@g0; vn=nweight@(velo*g1)
        if not deriv: return np.array([dn,dn+3*vn/w,(rc*dc+rb*db)/(rc+rb),u/max(abs(x),.2),phi,phi_c,w])
        return np.r_[-w*vc+3*phid,-vc+w*psi,-w*vb+3*phid,-vb+w*psi,0.,0.,du,ddu,dR.ravel(),dI.ravel(),dRn,dIn,dRg,dIg,phid]
    start[-1]=late_rhs(tm,start,False)[5]
    def event(t,s): return .1-abs(amplitude*late_rhs(t,s,False)[0]) if control=="full" else 1.
    event.terminal=True; event.direction=-1
    def chi_event(t,s): return .1-abs(amplitude*s[6])/max(abs(field(t)[0]),.2)
    chi_event.terminal=True; chi_event.direction=-1
    ls=solve_ivp(late_rhs,(tm,0),start,rtol=rtol,atol=rtol*.005,method='DOP853',max_step=.003,t_eval=np.linspace(tm,0,1601),events=[event,chi_event])
    assert ls.success; ts=ls.t; vals=np.array([late_rhs(t,s,False) for t,s in zip(ts,ls.y.T)])
    for tt,yy in zip(ls.t_events,ls.y_events):
        if tt.size and tt[-1]>ts[-1]: ts=np.r_[ts,tt[-1]]; vals=np.vstack([vals,late_rhs(tt[-1],yy[-1],False)])
    return dict(k=k,control=control,vev_factor=vev_factor,bias=bias,order=order,angles=angles,rtol=rtol,zdec=zdec,zmatch=zmatch,baryon_fraction=baryon_fraction,nu_fraction=nu_fraction,amplitude=amplitude,
                z_at_stop=float(np.expm1(-ts[-1])),reached_eta_threshold=bool(ls.t_events[0].size),reached_chi_threshold=bool(ls.t_events[1].size),
                max_eta=float(amplitude*max(abs(vals[:,0]))),max_chi_rel=float(amplitude*max(abs(vals[:,3]))),eta_rest_at_stop=float(amplitude*vals[-1,1]),max_matter=float(amplitude*max(abs(vals[:,2]))), matter_at_stop=float(amplitude*vals[-1,2]), phi_at_stop=float(amplitude*vals[-1,4]),
                max_density_constraint_rel=float(max(abs(vals[:,4]-vals[:,5])/np.maximum.reduce([abs(vals[:,4]),abs(vals[:,5]),np.full(len(vals),1e-10)]))),
                match_eta_number_per_psi=float(3*(R@wm)),match_matter_per_psi=float((1-baryon_fraction)*si[0]+baryon_fraction*si[2]),
                scope='Sudden decoupling, ideal tight coupling before z_dec, collisionless photons thereafter; no visibility function, polarization or full recombination. This is a sensitivity diagnostic.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--k',type=float,default=.03);p.add_argument('--radial',type=int,default=5);p.add_argument('--angles',type=int,default=96);p.add_argument('--rtol',type=float,default=1e-7);p.add_argument('--zdec',type=float,default=1090.);p.add_argument('--zmatch',type=float,default=12.);p.add_argument('--vev-factor',type=float,default=1.);p.add_argument('--bias',type=float,default=0.);p.add_argument('--control',choices=('full','lcdm','massless_eta'),default='full');p.add_argument('--output',required=True,type=Path);a=p.parse_args(); r=run(a.k,a.radial,a.angles,a.rtol,a.zdec,zmatch=a.zmatch,vev_factor=a.vev_factor,bias=a.bias,control=a.control);a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print(json.dumps(r,indent=2))
```
<!-- END FILE -->

### direct_pole_metric_pilot.py

<!-- FILE: direct_pole_metric_pilot.py -->
```python
"""Independent Newtonian-gauge pilot, late activation benchmark.
All standard matter is pressureless. Standard radiation is split into a
perfect fluid and an optional collisionless massless component. Eta is
collisionless (radial and angular quadrature); no recombination or baryon drag. This is a linear diagnostic, not a CLASS fit.
Reference background: direct-subthreshold-friedmann.json at 14999f92.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre

MPL=2.435e27; HBAR=6.582119569e-16; C=2.99792458e10; MPC=3.085677581e24
H0=.674*100*1e5/MPC*HBAR; OM=.315; OR=9.2e-5; RHO0=3*MPL*MPL*H0*H0
X=2.036548001081445e-10; VE=9.049809269161514e24; M0=.012664880357500833
T0=.0039033442679353545/36; CW=.2; SEED=1e-8; NETA=1; COEF=NETA/(2*np.pi**2)

def radial(order):
    edges=np.array([0.,.001,.01,.1,.4,1.,2.,4.,8.,16.,32.,56.])
    u,w=roots_legendre(order); h=np.diff(edges)[:,None]/2
    q=((edges[1:]+edges[:-1])[:,None]/2+h*u).ravel(); dq=(h*w).ravel()
    f=1/np.expm1(q); return q,dq,f

def potential(x,shape=0.,bias=0.):
    xx=x*x; lg=np.log(max(xx,1e-300))
    base=(X*(.25*(xx-1)**2+CW*(xx*xx*lg-1.5*xx*xx+2*xx-.5)),
            X*x*(xx-1+4*CW*(1+xx*(lg-1))),
            X*(3*xx-1+4*CW*(1+xx*(3*lg-1))))
    if bias != 0.: base=tuple(b+X*bias*v for b,v in zip(base,((x-1)**2,2*(x-1),2.)))
    if shape == 0.: return base
    # Even deformation vanishing at the origin and both vacuum minima.
    f=xx**2*(xx-1)**2/(1+xx**2)
    fp=4*x**3*(xx-1)*(xx**3+2*xx-1)/(1+xx**2)**2
    fpp=4*xx*(3*xx**6-xx**5+9*xx**4+9*xx**2-15*xx+3)/(1+xx**2)**3
    return tuple(b+X*shape*v for b,v in zip(base,(f,fp,fpp)))

def background(seed=SEED,shape=0.,vev_factor=1.,end=0.,stop_x=None,bias=0.):
    ve=VE*vev_factor
    q,dq,f=radial(40); weights=COEF*dq*q*q*f
    ys=np.geomspace(1e-10,1e5,1800); e=np.hypot(ys[:,None],q)
    fun=[PchipInterpolator(np.log(ys),np.log(v)) for v in [e@weights,(1/e)@(weights*q*q)/3,(1/e)@weights]]
    def local(t,x,dx):
        r=np.exp(-t); T=T0*r; yy=M0*x/T; logy=np.log(np.clip(abs(yy),ys[0],ys[-1])); er,ep,iv=[np.exp(a(logy)) for a in fun]
        er*=T**4; ep*=T**4; U,V,V2=potential(x,shape,bias); J=M0*M0*x*T*T*iv
        rm=RHO0*OM*r**3; rr=RHO0*OR*r**4
        H2=(rm+rr+RHO0*(1-OM-OR)+U+er)/(3*MPL*MPL-.5*ve*ve*dx*dx)
        K=.5*ve*ve*H2*dx*dx
        dh=-(rm+4*rr/3+er+ep+2*K)/(2*MPL*MPL*H2)
        return H2,dh,V+J,er,ep,rm,rr
    def rhs(t,s):
        x,dx=s; H2,dh,G,*_=local(t,x,dx)
        return dx,-(3+dh)*dx-G/(ve*ve*H2)
    def event(t,s): return abs(s[0])-stop_x
    event.terminal=True; event.direction=1
    sol=solve_ivp(rhs,(-np.log(1e14),end),[seed,0.],rtol=3e-10,atol=3e-20,method='DOP853',max_step=.025,dense_output=True,events=event if stop_x is not None else None)
    assert sol.success,sol.message
    return sol,local

def run(k,order,angles,rtol,metric=True,seed=SEED,zmatch=12.,zstart=1e6,amplitude=1e-5,stop_at=1.,free_rad_fraction=0.,shape=0.,evolve_phi=False,chi_threshold=None,compact=False,vev_factor=1.,bias=0.,initial_mode="adiabatic"):
    ve=VE*vev_factor
    assert 0 <= free_rad_fraction <= 1 and stop_at > 0 and amplitude > 0
    assert not evolve_phi or metric
    sol,loc=background(seed,shape,vev_factor,bias=bias); q,dq,f=radial(order); fd=f*(1+f)
    q2=COEF*dq*q*q*fd; number=np.sum(dq*q*q*f); nweight=dq*q*q*fd/number
    mu,wm=roots_legendre(angles); mu=(mu+1)/2; wm=wm/2; P2=(3*mu*mu-1)/2
    m1=wm*mu; m2=wm*P2
    t0=-np.log1p(zstart); tm=-np.log1p(zmatch)
    def env(t):
        x,dx=sol.sol(t); H2,dh,G,er,ep,rm,rr=loc(t,x,dx); T=T0*np.exp(-t); rho=3*MPL*MPL*H2; w=k*np.exp(-t)*C*HBAR/MPC/np.sqrt(H2)
        return x,dx,H2,dh,G,er,ep,rm,rr,T,rho,w
    # Until matching, eta is ultrarelativistic. Its perturbation factorizes
    # as -q f'_0 times one angular brightness; evaluate max m/T below.
    def early_rhs(t,s):
        dc,vc,dr,vr,u,du=s[:6]; R,I,Rn,In=s[6:].reshape(4,angles); x,dx,H2,dh,G,er,ep,rm,rr,T,rho,w=env(t)
        rn=rr*free_rad_fraction; rg=rr-rn
        b0=R@wm; b1=-I@m1; b2=-R@m2
        bn0=Rn@wm; bn1=-In@m1; bn2=-Rn@m2
        shear=8*(er*b2+rn*bn2)/3/rho; slip=4.5*shear/w**2
        mom=(rm*vc+4*rg*vr/3+4*er*b1+4*rn*bn1+ve**2*H2*w*dx*u)/rho
        U,V,V2=potential(x,shape,bias); d0=(rm*dc+rg*dr+4*er*b0+4*rn*bn0+ve**2*H2*dx*du+V*u)/rho
        psi=(-1.5/w**2*(d0+3*mom/w)-slip)/(1-1.5*ve**2*H2*dx*dx/rho/w**2); phi=psi+slip; phid=-psi+1.5*mom/w
        dR=phid+w*mu*I; dI=-w*mu*(R+psi)
        dRn=phid+w*mu*In; dIn=-w*mu*(Rn+psi)
        slipd= -2*slip+12/(rho*w*w)*(-er*(dR@m2)-rn*(dRn@m2))
        # eta direct mass curvature is retained at zero mass; angular
        # distribution feedback O(m^2) retained in relativistic limit.
        direct=M0*M0*T*T/12
        fb=2*M0*M0*x*T*T/12*b0
        ddu=-(3+dh)*du-(w*w+(V2+direct)/(ve**2*H2))*u-fb/(ve**2*H2)+dx*(4*phid-slipd)-2*G/(ve**2*H2)*psi
        return np.r_[-w*vc+3*phid,-vc+w*psi,-4*w*vr/3+4*phid,w*(dr/4+psi),du,ddu,dR,dI,dRn,dIn]
    x,dx,H2,dh,G,er,ep,rm,rr,T,rho,w=env(t0)
    psi=1.; v=w*psi/2; b1=v/3; b2=w*w*psi/30; R=-5*b2*P2; I=-3*b1*mu
    slip=12*(er+rr*free_rad_fraction)/rho*b2/w**2; phi=psi+slip
    # Enforce the density constraint instead of relying on leading-order
    # superhorizon density values, which would suffer cancellation at small w.
    mom=(rm*v+4*(rr+er)*v/3)/rho
    dR=(-2*w*w*phi/3-3*mom/w)*rho/(.75*rm+rr+er)
    initial=np.r_[.75*dR,v,dR,v,0.,0.,R+dR/4,I,R+dR/4,I]
    if initial_mode != "adiabatic":
        # Two compensated trial modes: initial eta number shift or absolute
        # delta x; vanishing initial total density and momentum (phi=psi=0).
        # The early eta is ultrarelativistic, as in early_rhs.
        assert initial_mode in ("eta_iso", "chi_iso")
        u=float(initial_mode=="chi_iso"); b=float(initial_mode=="eta_iso")/3
        U,V,V2=potential(x,shape,bias)
        dr=-(4*er*b+V*u)/rr; vr=-ve**2*H2*w*dx*u/(4*rr/3)
        initial=np.r_[0.,0.,dr,vr,u,0.,np.full(angles,b),np.zeros(angles),np.full(angles,dr/4),-vr*mu]
    early=solve_ivp(early_rhs,(t0,tm),initial,rtol=rtol,atol=rtol*.005,method='DOP853',max_step=.025,dense_output=True)
    assert early.success,early.message
    tcheck=np.linspace(t0,tm,400); max_y=float(np.max(abs(sol.sol(tcheck)[0])*M0*np.exp(tcheck)/T0))
    # At match expand brightness into the full eta distribution; subsequent
    # evolution uses the exact signed mass gradient and relativistic stress.
    si=early.y[:,-1]; R,I,Rn,In=si[6:].reshape(4,angles)
    start=np.r_[si[:6],(q[:,None]*R).ravel(),(q[:,None]*I).ravel(),Rn,In]
    def evaluate(t,s,deriv=True):
        dc,vc,dr,vr,u,du=s[:6]; R,I=s[6:6+2*q.size*angles].reshape(2,q.size,angles); Rn,In=s[6+2*q.size*angles:6+2*q.size*angles+2*angles].reshape(2,angles)
        x,dx,H2,dh,G,er,ep,rm,rr,T,rho,w=env(t); U,V,V2=potential(x,shape,bias)
        rn=rr*free_rad_fraction; rg=rr-rn; bn0=Rn@wm; bn1=-In@m1; bn2=-Rn@m2
        y=M0*x/T; e=np.hypot(q,y); velo=q/e; de=M0*M0*x/(T*T*e)
        g0=R@wm; g1=-I@m1; g2=-R@m2
        direct=COEF*T*T*M0*M0*np.sum(dq*q**4*f/e**3)
        J=COEF*T**4*np.sum(dq*q*q*f*de)
        # delta rho at fixed momentum plus direct varying-mass energy weight
        drho=T**4*np.sum(q2*e*g0)+J*u
        eta_mom=T**4*np.sum(q2*q*g1)
        shear=(2*T**4/3*np.sum(q2*q*q/e*g2)+8*rn*bn2/3)/rho; slip=4.5*shear/w**2
        mom=(rm*vc+4*rg*vr/3+4*rn*bn1+eta_mom+ve**2*H2*w*dx*u)/rho
        d0=(rm*dc+rg*dr+4*rn*bn0+drho+ve**2*H2*dx*du+V*u)/rho
        psi_c=(-1.5/w**2*(d0+3*mom/w)-slip)/(1-1.5*ve**2*H2*dx*dx/rho/w**2)
        phi_c=psi_c+slip
        if metric:
            phi=s[-1] if evolve_phi else phi_c
            psi=phi-slip; phid=-psi+1.5*mom/w
        else: phi=psi=phid=0.
        dR=q[:,None]*phid+w*mu[None,:]*velo[:,None]*I
        dI=-w*mu[None,:]*(velo[:,None]*R+(e*psi+de*u)[:,None])
        g2d=-dR@m2
        dRn=phid+w*mu*In; dIn=-w*mu*(Rn+psi)
        edot=(M0/T)**2*x*(x+dx)/e
        slipd=-2*slip+3*T**4/(rho*w*w)*np.sum(q2*q*q/e*(g2d-g2*edot/e))-12*rn/(rho*w*w)*(dRn@m2) if metric else 0.
        fb=T**4*np.sum(q2*de*g0)
        ddu=-(3+dh)*du-(w*w+(V2+direct)/(ve**2*H2))*u-fb/(ve**2*H2)+dx*(4*phid-slipd)-2*(V+J)/(ve**2*H2)*psi
        dn=nweight@g0
        vn=nweight@(velo*g1)
        if not deriv: return np.array([x,dx,dn,dc,phi,psi,phid,slip,drho/rho,w,u,dn+3*vn/w,er/rho,phi_c])
        rhs=np.r_[-w*vc+3*phid,-vc+w*psi,-4*w*vr/3+4*phid,w*(dr/4+psi),du,ddu,dR.ravel(),dI.ravel(),dRn,dIn]
        return np.r_[rhs,phid] if evolve_phi else rhs
    if evolve_phi:
        start=np.r_[start,0.]
        start[-1]=evaluate(tm,start,False)[-1]
    def event(t,s): return stop_at-abs(amplitude*evaluate(t,s,False)[2])
    event.terminal=True; event.direction=-1
    def chi_event(t,s): return chi_threshold-amplitude*abs(s[4])/max(abs(sol.sol(t)[0]),.2)
    chi_event.terminal=True; chi_event.direction=-1
    events=[event] if chi_threshold is None else [event,chi_event]
    late=solve_ivp(evaluate,(tm,0),start,rtol=rtol,atol=rtol*.005,method='DOP853',max_step=.003,events=events,dense_output=not compact,t_eval=np.linspace(tm,0,1601) if compact else None)
    assert late.success,late.message
    if compact:
        ts=late.t; vals=np.array([evaluate(t,y,False) for t,y in zip(ts,late.y.T)])
        for tt,yy in zip(late.t_events,late.y_events):
            if tt.size and tt[-1]>ts[-1]:
                ts=np.r_[ts,tt[-1]]; vals=np.vstack([vals,evaluate(tt[-1],yy[-1],False)])
        err=scale=None
    else:
        ts=np.linspace(tm,late.t[-1],4001); vals=np.concatenate([np.array([evaluate(t,s,False) for t,s in zip(chunk,late.sol(chunk).T)]) for chunk in np.array_split(ts,40)])
        # Independent finite-difference check of the momentum constraint.
        eps=2e-6; idx=np.linspace(2,len(ts)-3,300,dtype=int); tsel=ts[idx]
        fd_phi=np.array([(evaluate(t+eps,late.sol(t+eps),False)[4]-evaluate(t-eps,late.sol(t-eps),False)[4])/(2*eps) for t in tsel])
        err=np.abs(fd_phi-vals[idx,6]); scale=np.maximum.reduce([np.abs(fd_phi),np.abs(vals[idx,6]),np.abs(vals[idx,4]),np.full(idx.size,1e-10)])
    report=dict(k=k,radial_order=order,angles=angles,rtol=rtol,metric=metric,evolve_phi=evolve_phi,compact=compact,seed=seed,shape=shape,vev_factor=vev_factor,bias=bias,zmatch=zmatch,zstart=zstart,initial_mode=initial_mode,initial_amplitude=amplitude,psi_initial=amplitude if initial_mode=="adiabatic" else 0.,free_standard_radiation_fraction=free_rad_fraction,eta_number_threshold=stop_at,
                early_max_eta_m_over_T=max_y,initial_radiation_delta_per_psi=float(initial[2]),
                match_delta_m_per_psi=float(si[0]),match_eta_delta_n_per_psi=float(3*(R@wm)),
                z_at_stop=float(np.expm1(-ts[-1])),reached_eta_threshold=bool(late.t_events[0].size),chi_threshold=chi_threshold,reached_chi_threshold=bool(chi_threshold is not None and late.t_events[1].size),
                max_eta_number_fraction=float(np.max(abs(vals[:,2]))*amplitude),max_matter_fraction=float(np.max(abs(vals[:,3]))*amplitude),
                max_phi=float(np.max(abs(vals[:,4]))*amplitude),max_delta_x=float(np.max(abs(vals[:,10]))*amplitude),max_delta_x_over_x_floor=float(amplitude*np.max(abs(vals[:,10])/np.maximum(abs(vals[:,0]),.2))),
                at_stop=dict(k_over_aH=float(vals[-1,9]),eta_number=float(vals[-1,2]*amplitude),eta_number_rest_frame=float(vals[-1,11]*amplitude),eta_background_fraction=float(vals[-1,12]),eta_density_over_total=float(vals[-1,8]*amplitude),delta_x_over_x=float(amplitude*vals[-1,10]/vals[-1,0])),
                max_constraint_derivative_rel=float(np.max(err/scale)) if metric and not compact else None,rms_constraint_derivative_rel=float(np.sqrt(np.mean((err/scale)**2))) if metric and not compact else None,
                max_density_constraint_rel=float(np.max(np.abs(vals[:,4]-vals[:,13])/np.maximum.reduce([abs(vals[:,4]),abs(vals[:,13]),np.full(len(vals),1e-10)]))) if evolve_phi else None,
                early_nfev=early.nfev,late_nfev=late.nfev,
                background_today_x=float(sol.sol(0)[0]),background_max_x=float(np.max(sol.sol(np.linspace(-np.log(101),0,16000))[0])),
                scope='Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like or compensated trial mode; no primordial spectrum imposed. Stop at eta number contrast threshold.')
    return report

def main():
    p=argparse.ArgumentParser(); p.add_argument('--k',type=float,default=.03);p.add_argument('--radial',type=int,default=5);p.add_argument('--angles',type=int,default=48);p.add_argument('--rtol',type=float,default=2e-7);p.add_argument('--zmatch',type=float,default=12.);p.add_argument('--seed',type=float,default=SEED);p.add_argument('--amplitude',type=float,default=1e-5);p.add_argument('--stop-at',type=float,default=1.);p.add_argument('--no-late-metric',action='store_true');p.add_argument('--free-radiation',type=float,default=0.);p.add_argument('--zstart',type=float,default=1e6);p.add_argument('--shape',type=float,default=0.);p.add_argument('--evolve-phi',action='store_true');p.add_argument('--chi-threshold',type=float);p.add_argument('--compact',action='store_true');p.add_argument('--vev-factor',type=float,default=1.);p.add_argument('--bias',type=float,default=0.);p.add_argument('--initial-mode',choices=("adiabatic","eta_iso","chi_iso"),default="adiabatic");p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=run(a.k,a.radial,a.angles,a.rtol,not a.no_late_metric,a.seed,a.zmatch,zstart=a.zstart,amplitude=a.amplitude,stop_at=a.stop_at,free_rad_fraction=a.free_radiation,shape=a.shape,evolve_phi=a.evolve_phi,chi_threshold=a.chi_threshold,compact=a.compact,vev_factor=a.vev_factor,bias=a.bias,initial_mode=a.initial_mode)
    a.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print(json.dumps(r,indent=2),flush=True)
if __name__=='__main__':main()
```
<!-- END FILE -->

### direct_pole_slow_trigger.py

<!-- FILE: direct_pole_slow_trigger.py -->
```python
"""Retune the initial trigger displacement to preserve half activation.
The potential energy and eta mass/temperature are unchanged; vev changes
both the field inertia and the eta coupling. Contact exposure is a proxy,
not the full Bose collision operator. The early proxy is not a lower bound.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import trapezoid
from direct_pole_metric_pilot import background,potential,radial,VE,M0,T0,COEF,MPL

TARGET_Z=4.3445; TARGET=-np.log1p(TARGET_Z); HIT=.5**.25
BASE_RATIO=52.43035455726708; BASE_PROXY=.008746239106334854

def tune(factor,bias_mode=False,target_z=TARGET_Z):
    target=-np.log1p(target_z)
    def mismatch(ls):
        seed=1e-8 if bias_mode else np.exp(ls); bias=np.exp(ls) if bias_mode else 0.
        sol,_=background(seed=seed,vev_factor=factor,end=target,stop_x=HIT,bias=bias)
        if sol.t_events[0].size: return (target-sol.t[-1])*max(abs(sol.y[1,-1]),.001)
        return abs(sol.y[0,-1])-HIT
    lo,hi=np.log(1e-14),np.log(.8)
    if mismatch(lo)*mismatch(hi)>0: return dict(vev_factor=factor,tuned=False)
    tuned=float(np.exp(brentq(mismatch,lo,hi,xtol=5e-6)))
    seed=1e-8 if bias_mode else tuned; bias=tuned if bias_mode else 0.
    sol,loc=background(seed=seed,vev_factor=factor,bias=bias); ve=VE*factor
    t=np.r_[np.linspace(-np.log(1e5+1),-np.log(101),6000)[:-1],np.linspace(-np.log(101),0,14000)]
    x,dx=sol.sol(t); H2,dh,G,er,ep,rm,rr=np.array([loc(a,b,c) for a,b,c in zip(t,x,dx)]).T
    U,V,V2=np.array([potential(a,bias=bias) for a in x]).T
    T=T0*np.exp(-t); q,dq,f=radial(8); e=np.hypot(q[None,:],(M0*x/T)[:,None]); direct=COEF*M0*M0*T*T*(1/e**3 @ (dq*q**4*f))
    c=(V2+direct)/(ve*ve*H2); late=t>=-np.log(101); halo=t>=-np.log(3); pre=t<=-np.log(101)+1e-13
    after=np.flatnonzero(abs(x)**4>=.5); assert after.size
    positive=trapezoid(np.maximum((np.sqrt(np.maximum(9-4*c[after[0]:],0))-3)/2,0),t[after[0]:])
    A4=trapezoid(x[halo]**4/np.sqrt(H2[halo]),t[halo]); weight=np.exp(-4*t)*x**4/np.sqrt(H2)
    early=trapezoid(weight[pre],t[pre])/A4; late_ratio=trapezoid(weight[late],t[late])/A4
    # Representative pair cross-section suppression only, v_rel=v_peak,0/a.
    # Not a rigorous bound on the quantum collision operator or full PSD.
    v0=17.834903110469433; vref=np.sqrt(2)*44.1e5; pole=(10/(10+(v0*np.exp(-t)/vref)**2))**2
    early_pair=trapezoid((weight*pole)[pre],t[pre])/A4
    extra=(.5*ve*ve*H2*dx*dx+U+er)/(3*MPL*MPL*H2)
    return dict(vev_factor=factor,tuned=True,seed=seed,bias=bias,target_z=target_z,z_half=float(np.expm1(-t[after[0]])),max_extra_zle100=float(max(extra[late])),
                positive_local_root_integral=float(positive),minimum_curvature_after_half=float(min(c[after[0]:])),
                late_contact_proxy=float(BASE_PROXY*late_ratio/BASE_RATIO),early_contact_proxy_100_to_1e5=float(BASE_PROXY*early/BASE_RATIO),
                early_representative_pair_proxy=float(BASE_PROXY*early_pair/BASE_RATIO),early_x_z1e5=float(x[0]),
                max_abs_x=float(max(abs(x))),today_x=float(x[-1]),activation_integral=float(A4),
                scope='Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--factor',type=float,nargs='+',required=True);p.add_argument('--bias-mode',action='store_true');p.add_argument('--target-z',type=float,default=TARGET_Z);p.add_argument('--output',required=True,type=Path);a=p.parse_args();rows=[]
    for fac in a.factor:
        r=tune(fac,a.bias_mode,a.target_z);rows.append(r);a.output.write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n');print(json.dumps(r),flush=True)
```
<!-- END FILE -->

### direct_pole_trigger_shape.py

<!-- FILE: direct_pole_trigger_shape.py -->
```python
"""Exploratory potential shape screen, not a microscopic trigger model."""
import json,argparse
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from direct_pole_metric_pilot import background,potential,radial,M0,T0,X,VE,MPL,COEF

def screen(shape):
    xx=np.r_[0.,np.geomspace(1e-8,.99,4000),np.linspace(.99,1.01,1200),np.geomspace(1.01,100,4000)]
    uu,v,vv=np.array([potential(x,shape) for x in xx]).T/X
    # Numerical finite field monotonicity/positivity screen, not a proof.
    valid=bool(np.min(uu)>=-1e-12 and np.max(v[xx<1])<=1e-12 and np.min(v[xx>1])>=-1e-12)
    if not valid: return dict(shape=shape,finite_field_potential_screen=False,min_U_over_X=float(min(uu)))
    sol,loc=background(shape=shape)
    t=np.linspace(-np.log(101),0.,18000); x,dx=sol.sol(t)
    H2,dh,G,er,ep,rm,rr=np.array([loc(a,b,c) for a,b,c in zip(t,x,dx)]).T
    U,V,V2=np.array([potential(a,shape) for a in x]).T
    q,dq,f=radial(8); T=T0*np.exp(-t); e=np.hypot(q[None,:],(M0*x/T)[:,None])
    direct=COEF*M0*M0*T*T*(1/e**3 @ (dq*q**4*f))
    c=(V2+direct)/(VE**2*H2)
    i=np.flatnonzero(abs(x)**4>=.5)[0]; tlate=t>=t[i]; halo=t>=-np.log(3)
    # Local frozen curvature-only screen: not accumulated phase growth.
    proxy=trapezoid(np.maximum((np.sqrt(np.maximum(9-4*c[tlate],0))-3)/2,0),t[tlate])
    A4=trapezoid(x[halo]**4/np.sqrt(H2[halo]),t[halo]); cosmic=trapezoid(np.exp(-4*t)*x**4/np.sqrt(H2),t)
    extra=(.5*VE*VE*H2*dx*dx+U+er)/(3*MPL*MPL*H2)
    return dict(shape=shape,finite_field_potential_screen=True,z_half=float(np.expm1(-t[i])),max_extra=float(max(extra)),min_curvature_after_half=float(min(c[tlate])),positive_local_root_integral=float(proxy),activation_integral=float(A4),cosmic_to_activation=float(cosmic/A4),max_x=float(max(x)),min_x_after_half=float(min(x[tlate])),today_x=float(x[-1]))

def main():
    p=argparse.ArgumentParser(); p.add_argument('--output',required=True,type=Path); a=p.parse_args()
    rows=[screen(s) for s in (-.4,-.3,-.2,-.14,-.1,0.,.1,.2,.3,.4)]
    ref=next(row for row in rows if row['shape']==0.)
    for row in rows:
        if row['finite_field_potential_screen']:
            row['retuned_cosmic_proxy']=0.008746239106334854*row['cosmic_to_activation']/ref['cosmic_to_activation']
    a.output.write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
```
<!-- END FILE -->

### elastic-base240.json

<!-- FILE: elastic-base240.json -->
```json
{
  "mesh": 240,
  "n": 220,
  "events": 3482680,
  "vev_factor": 1.0,
  "bias": 0.0,
  "slices": 1,
  "nodes": 16,
  "rtol": 1e-06,
  "reference_tail": -0.008102174905421533,
  "tau": 1.5194994680886315e-10,
  "tau_ratio": 1.0,
  "A4_ratio": 1.0,
  "sigma0": 0.22089528043464016,
  "relative_moments": [
    -1.1392580778757377e-07,
    -0.008102173869313112,
    0.0021997976187422275,
    0.0013832632287753555
  ],
  "max_number_drift": 6.661338147750939e-16,
  "max_energy_drift": 2.220446049250313e-16,
  "max_negative_number_fraction": 1.7508978152768644e-23,
  "detailed_balance_residual": 5.91658271853033e-16,
  "steps": [
    {
      "slice": 1,
      "nfev": 98,
      "moment_changes": [
        -1.1392580778757377e-07,
        -0.008102173869313112,
        0.0021997976187422275,
        0.0013832632287753555
      ]
    }
  ],
  "scope": "Conditional post-production homogeneous NR Bose 2-to-2 scattering on reconstructed 220-bin PSD; same graph convention as prior archive. No early collision or in-medium calculation."
}
```
<!-- END FILE -->

### elastic-bias240-sliced.json

<!-- FILE: elastic-bias240-sliced.json -->
```json
{
  "mesh": 240,
  "n": 220,
  "events": 3482680,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "slices": 4,
  "nodes": 12,
  "rtol": 3e-07,
  "reference_tail": -0.008102174905421533,
  "tau": 1.4558174594712668e-10,
  "tau_ratio": 0.9580901408952318,
  "A4_ratio": 0.9933716419965258,
  "sigma0": 0.22236922325533096,
  "relative_moments": [
    -1.0893339963402582e-07,
    -0.007786545690242175,
    0.0021096985163995274,
    0.001325090393938133
  ],
  "max_number_drift": 6.661338147750939e-16,
  "max_energy_drift": 2.220446049250313e-16,
  "max_negative_number_fraction": 1.095370321503643e-16,
  "detailed_balance_residual": 5.915873276518954e-16,
  "steps": [
    {
      "slice": 1,
      "nfev": 88,
      "moment_changes": [
        -2.706682322717313e-08,
        -0.0020571355373100175,
        0.0005365819592473287,
        0.00033040946833229334
      ]
    },
    {
      "slice": 2,
      "nfev": 40,
      "moment_changes": [
        -5.4354069467521526e-08,
        -0.004037371162492587,
        0.0010669837892101075,
        0.0006613961495112175
      ]
    },
    {
      "slice": 3,
      "nfev": 38,
      "moment_changes": [
        -8.16430145711422e-08,
        -0.005945481439681877,
        0.001591317498263356,
        0.0009929572234579354
      ]
    },
    {
      "slice": 4,
      "nfev": 38,
      "moment_changes": [
        -1.0893339963402582e-07,
        -0.007786545690242175,
        0.0021096985163995274,
        0.001325090393938133
      ]
    }
  ],
  "scope": "Conditional post-production homogeneous NR Bose 2-to-2 scattering on reconstructed 220-bin PSD; same graph convention as prior archive. No early collision or in-medium calculation."
}
```
<!-- END FILE -->

### elastic-bias240.json

<!-- FILE: elastic-bias240.json -->
```json
{
  "mesh": 240,
  "n": 220,
  "events": 3482680,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "slices": 1,
  "nodes": 16,
  "rtol": 1e-06,
  "reference_tail": -0.008102174905421533,
  "tau": 1.4558174594712668e-10,
  "tau_ratio": 0.9580901408952318,
  "A4_ratio": 0.9933716419965258,
  "sigma0": 0.22236922325533096,
  "relative_moments": [
    -1.0915287007584595e-07,
    -0.007786537942998439,
    0.0021096983715429563,
    0.0013250903901160793
  ],
  "max_number_drift": 6.661338147750939e-16,
  "max_energy_drift": 2.220446049250313e-16,
  "max_negative_number_fraction": 1.603986570373627e-23,
  "detailed_balance_residual": 5.916717067748623e-16,
  "steps": [
    {
      "slice": 1,
      "nfev": 96,
      "moment_changes": [
        -1.0915287007584595e-07,
        -0.007786537942998439,
        0.0021096983715429563,
        0.0013250903901160793
      ]
    }
  ],
  "scope": "Conditional post-production homogeneous NR Bose 2-to-2 scattering on reconstructed 220-bin PSD; same graph convention as prior archive. No early collision or in-medium calculation."
}
```
<!-- END FILE -->

### elastic-bias320.json

<!-- FILE: elastic-bias320.json -->
```json
{
  "mesh": 320,
  "n": 220,
  "events": 3479014,
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "slices": 1,
  "nodes": 16,
  "rtol": 1e-06,
  "reference_tail": -0.00810603014148481,
  "tau": 1.4376317459358914e-10,
  "tau_ratio": 0.9580901408952318,
  "A4_ratio": 0.9933716419965258,
  "sigma0": 0.22236922325533096,
  "relative_moments": [
    -1.0846607534364949e-07,
    -0.007790656708083765,
    0.0017775419692096506,
    0.0014280187809658784
  ],
  "max_number_drift": 6.661338147750939e-16,
  "max_energy_drift": 2.220446049250313e-16,
  "max_negative_number_fraction": 3.532990660663646e-24,
  "detailed_balance_residual": 6.332898541647532e-16,
  "steps": [
    {
      "slice": 1,
      "nfev": 106,
      "moment_changes": [
        -1.0846607534364949e-07,
        -0.007790656708083765,
        0.0017775419692096506,
        0.0014280187809658784
      ]
    }
  ],
  "scope": "Conditional post-production homogeneous NR Bose 2-to-2 scattering on reconstructed 220-bin PSD; same graph convention as prior archive. No early collision or in-medium calculation."
}
```
<!-- END FILE -->

### evolved-phi-fine.json

<!-- FILE: evolved-phi-fine.json -->
```json
{
  "k": 0.03,
  "radial_order": 6,
  "angles": 96,
  "rtol": 3e-08,
  "metric": true,
  "evolve_phi": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4086214687368,
  "match_eta_delta_n_per_psi": -1.1235193692035854,
  "z_at_stop": 3.390754766439937,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.10000000000000007,
  "max_matter_fraction": 0.03133711477301986,
  "max_phi": 3.8164583555616655e-06,
  "max_delta_x": 0.020241695394288838,
  "at_stop": {
    "k_over_aH": 110.96188997106198,
    "eta_number": -0.10000000000000007,
    "eta_number_rest_frame": -0.09905977563464329,
    "eta_background_fraction": 9.636801396595778e-05,
    "eta_density_over_total": -6.75205958580017e-06,
    "delta_x_over_x": 0.03195172269846626
  },
  "max_constraint_derivative_rel": 1.950674366154076e-06,
  "rms_constraint_derivative_rel": 1.453455905375788e-07,
  "max_density_constraint_rel": 1.2474709264088631e-08,
  "early_nfev": 8849,
  "late_nfev": 6575,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### evolved-phi.json

<!-- FILE: evolved-phi.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 80,
  "rtol": 8e-08,
  "metric": true,
  "evolve_phi": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4087498615208,
  "match_eta_delta_n_per_psi": -1.1235194997239533,
  "z_at_stop": 3.390754814484032,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.10000000000000013,
  "max_matter_fraction": 0.03133711819687248,
  "max_phi": 3.816458839989201e-06,
  "max_delta_x": 0.02024169896477308,
  "at_stop": {
    "k_over_aH": 110.96188943772837,
    "eta_number": -0.10000000000000013,
    "eta_number_rest_frame": -0.09905977634162481,
    "eta_background_fraction": 9.636804438065252e-05,
    "eta_density_over_total": -6.7520692096929e-06,
    "delta_x_over_x": 0.0319516397406544
  },
  "max_constraint_derivative_rel": 3.3721834371126275e-06,
  "rms_constraint_derivative_rel": 2.767919744776597e-07,
  "max_density_constraint_rel": 8.489893106268443e-08,
  "early_nfev": 8339,
  "late_nfev": 6224,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### free-rad-early.json

<!-- FILE: free-rad-early.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "seed": 1e-08,
  "zmatch": 10.0,
  "zstart": 2000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.48782898058577e-10,
  "initial_radiation_delta_per_psi": -2.000037080334355,
  "match_delta_m_per_psi": -1264.2748261694717,
  "match_eta_delta_n_per_psi": -1.1127804282289473,
  "z_at_stop": 3.390772891942069,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.10000000000000041,
  "max_matter_fraction": 0.03134197743292831,
  "max_phi": 3.817074704968816e-06,
  "max_delta_x": 0.02024142480949959,
  "max_constraint_derivative_rel": 0.00013699425704715004,
  "rms_constraint_derivative_rel": 9.980707490451143e-06,
  "early_nfev": 10721,
  "late_nfev": 5378,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### free-rad-fine.json

<!-- FILE: free-rad-fine.json -->
```json
{
  "k": 0.03,
  "radial_order": 6,
  "angles": 96,
  "rtol": 8e-08,
  "metric": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4087498443564,
  "match_eta_delta_n_per_psi": -1.1235194997733842,
  "z_at_stop": 3.390754788177674,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.09999999999999834,
  "max_matter_fraction": 0.03133711837982009,
  "max_phi": 3.816458862008217e-06,
  "max_delta_x": 0.020241698154782172,
  "at_stop": {
    "k_over_aH": 110.9618897297531,
    "eta_number": -0.09999999999999834,
    "eta_number_rest_frame": -0.09905977593529014,
    "eta_background_fraction": 9.636802772717908e-05,
    "eta_density_over_total": -6.752063212405174e-06,
    "delta_x_over_x": 0.03195169339202129
  },
  "max_constraint_derivative_rel": 2.926569354710466e-05,
  "rms_constraint_derivative_rel": 1.9423698642742106e-06,
  "early_nfev": 8339,
  "late_nfev": 6260,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### free-rad.json

<!-- FILE: free-rad.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "seed": 1e-08,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4088035686345,
  "match_eta_delta_n_per_psi": -1.123519554436656,
  "z_at_stop": 3.39075482557096,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.09999999999999963,
  "max_matter_fraction": 0.03133711969199764,
  "max_phi": 3.816459033525467e-06,
  "max_delta_x": 0.02024170044668394,
  "max_constraint_derivative_rel": 0.0001287236889651407,
  "rms_constraint_derivative_rel": 7.700072416169424e-06,
  "early_nfev": 8234,
  "late_nfev": 6155,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### growth-k03-full.json

<!-- FILE: growth-k03-full.json -->
```json
{
  "k": 0.03,
  "control": "full",
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "order": 5,
  "angles": 160,
  "rtol": 1e-07,
  "zdec": 1090.0,
  "zmatch": 100.0,
  "baryon_fraction": 0.15555555555555556,
  "nu_fraction": 0.409,
  "amplitude": 1e-05,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "reached_chi_threshold": false,
  "max_eta": 0.0082609768128561,
  "max_chi_rel": 4.361213124599871e-07,
  "eta_rest_at_stop": -0.008259614516103598,
  "max_matter": 0.08811976162003239,
  "matter_at_stop": -0.08811976162003239,
  "phi_at_stop": 2.33819862926446e-06,
  "max_density_constraint_rel": 3.9751448223307225e-13,
  "match_eta_number_per_psi": -0.8856598814067413,
  "match_matter_per_psi": -115.3291625476508,
  "scope": "Sudden decoupling, ideal tight coupling before z_dec, collisionless photons thereafter; no visibility function, polarization or full recombination. This is a sensitivity diagnostic."
}
```
<!-- END FILE -->

### growth-k03-lcdm.json

<!-- FILE: growth-k03-lcdm.json -->
```json
{
  "k": 0.03,
  "control": "lcdm",
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "order": 5,
  "angles": 160,
  "rtol": 1e-07,
  "zdec": 1090.0,
  "zmatch": 100.0,
  "baryon_fraction": 0.15555555555555556,
  "nu_fraction": 0.409,
  "amplitude": 1e-05,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "reached_chi_threshold": false,
  "max_eta": 9.73848820758697e-06,
  "max_chi_rel": 0.0,
  "eta_rest_at_stop": -7.249573537475955e-06,
  "max_matter": 0.08987456357018311,
  "matter_at_stop": -0.08987456357018311,
  "phi_at_stop": 2.384711977797961e-06,
  "max_density_constraint_rel": 2.376117345770304e-15,
  "match_eta_number_per_psi": -0.8903578451086226,
  "match_matter_per_psi": -115.71811074282473,
  "scope": "Sudden decoupling, ideal tight coupling before z_dec, collisionless photons thereafter; no visibility function, polarization or full recombination. This is a sensitivity diagnostic."
}
```
<!-- END FILE -->

### growth-k03-massless_eta.json

<!-- FILE: growth-k03-massless_eta.json -->
```json
{
  "k": 0.03,
  "control": "massless_eta",
  "vev_factor": 10.0,
  "bias": 0.0032765520306085277,
  "order": 5,
  "angles": 160,
  "rtol": 1e-07,
  "zdec": 1090.0,
  "zmatch": 100.0,
  "baryon_fraction": 0.15555555555555556,
  "nu_fraction": 0.409,
  "amplitude": 1e-05,
  "z_at_stop": -0.0,
  "reached_eta_threshold": false,
  "reached_chi_threshold": false,
  "max_eta": 9.708027602880782e-06,
  "max_chi_rel": 0.0,
  "eta_rest_at_stop": -7.2210516432152686e-06,
  "max_matter": 0.08955027195139473,
  "matter_at_stop": -0.08955027195139473,
  "phi_at_stop": 2.3761072960460936e-06,
  "max_density_constraint_rel": 3.9848207374434116e-13,
  "match_eta_number_per_psi": -0.8856604040049789,
  "match_matter_per_psi": -115.32921298320977,
  "scope": "Sudden decoupling, ideal tight coupling before z_dec, collisionless photons thereafter; no visibility function, polarization or full recombination. This is a sensitivity diagnostic."
}
```
<!-- END FILE -->

### k02.json

<!-- FILE: k02.json -->
```json
{
  "k": 0.02,
  "radial_order": 4,
  "angles": 64,
  "rtol": 2e-07,
  "metric": true,
  "seed": 1e-08,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.0,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000056740619312,
  "match_delta_m_per_psi": -592.5502781701796,
  "match_eta_delta_n_per_psi": -1.3580304981509188,
  "z_at_stop": 3.2907865067045137,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.10000000000000019,
  "max_matter_fraction": 0.01773128228227702,
  "max_phi": 5.106205010772511e-06,
  "max_delta_x": 0.054571242692681585,
  "max_constraint_derivative_rel": 0.00043107852111758127,
  "rms_constraint_derivative_rel": 2.8213184104806327e-05,
  "early_nfev": 8588,
  "late_nfev": 6242,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### k04.json

<!-- FILE: k04.json -->
```json
{
  "k": 0.04,
  "radial_order": 4,
  "angles": 96,
  "rtol": 2e-07,
  "metric": true,
  "seed": 1e-08,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.0,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.0002269624772473,
  "match_delta_m_per_psi": -1391.755120804984,
  "match_eta_delta_n_per_psi": -0.7992572374988032,
  "z_at_stop": 3.163705955775518,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.09999999999999959,
  "max_matter_fraction": 0.04291603529927072,
  "max_phi": 2.713355071635976e-06,
  "max_delta_x": 0.011469202840788263,
  "max_constraint_derivative_rel": 3.869083195128205e-05,
  "rms_constraint_derivative_rel": 3.504300495103282e-06,
  "early_nfev": 7919,
  "late_nfev": 6512,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### pilot-1.json

<!-- FILE: pilot-1.json -->
```json
{
  "k": 0.03,
  "radial_order": 4,
  "angles": 48,
  "rtol": 3e-07,
  "metric": true,
  "seed": 1e-08,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.0,
  "eta_number_threshold": 1.0,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.0001276663934515,
  "match_delta_m_per_psi": -1000.091787106264,
  "match_eta_delta_n_per_psi": -1.0498422738243895,
  "z_at_stop": 3.1306240584753073,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 1.0000000000000027,
  "max_matter_fraction": 0.03108014136438697,
  "max_phi": 5.390328154045022e-06,
  "max_delta_x": 0.3033819953477027,
  "max_constraint_derivative_rel": 0.00045112946198142825,
  "rms_constraint_derivative_rel": 4.018547489367681e-05,
  "early_nfev": 7814,
  "late_nfev": 6497,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### pilot-2.json

<!-- FILE: pilot-2.json -->
```json
{
  "k": 0.03,
  "radial_order": 6,
  "angles": 80,
  "rtol": 8e-08,
  "metric": true,
  "seed": 1e-08,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.0,
  "eta_number_threshold": 1.0,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.0001276663934515,
  "match_delta_m_per_psi": -1000.0911431215327,
  "match_eta_delta_n_per_psi": -1.0498416158417327,
  "z_at_stop": 3.131620155099254,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.9999999999999966,
  "max_matter_fraction": 0.03107284945780168,
  "max_phi": 5.385572498166782e-06,
  "max_delta_x": 0.30633483591542066,
  "max_constraint_derivative_rel": 0.0001508672673018676,
  "rms_constraint_derivative_rel": 9.557823455266002e-06,
  "early_nfev": 8339,
  "late_nfev": 6917,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### psd_reference.py

<!-- FILE: psd_reference.py -->
```python
"""Compact transcription of public prior benchmark input at 14999f92.
See kinetics/results/direct-pole-evolution-{240,320}.json. Grid reconstructed
from geometric endpoints and nearest original index; relative face error
<=2.2e-13 (240), <=8e-15 (320) against the archived values.
"""
import numpy as np

DATA={
240:(0.00001506234557978815,1.0964998669026034,1.5194994680886315e-10,17.834902767879523,845992330.2269042,-0.008102174905421533,
'''3.23479450900695e-19,2.86737007096245e-14,9.87997900194025e-13,1.46465861938221e-11,1.72153812734187e-10,1.64706326958764e-9,1.29010725309488e-8,8.57655113599411e-8,4.86507610802761e-7,0.00000238384244304437,0.0000102177102965560,0.0000387739343535104,0.000286337586216272,0.00112814719744310,0.00289683118689333,0.00688868062442572,0.0152893126755936,0.0319204349415787,0.0632273885429313,0.119947924682214,0.220344697169392,0.397169632360139,0.713769836253527,1.97769197128673,5.10446106466832,11.6881055620021,31.0456363751531,98.9149218140978,385.285271393151,1830.94204272447,10366.9816890932,66863.9137622580,461131.794950962,3165789.88486826,70700228.3253159,477288270.752775,1596906756.84427,3959913277.23624,7219096614.58250,9810758617.68927,10250312551.1056,8742104537.25922,6425790566.39999,4283223486.17699,2735659912.20911,1420127279.78946,814303613.459981,594978133.714366,449047933.681795,344373699.145609,264580755.401852,201721617.768319,151861424.169498,112752980.440966,82710876.8892533,60065380.7430276,36358914.5894962,22206982.5659612,15817802.6475299,11222748.2290013,7932441.70426200,5576942.19221970,3895407.12471318,2701653.57685877,1856718.26634728,1262940.80034396,848872.045065024,451885.502764979,235880.064672908,148506.134009526,91567.1609374851,55263.2930073157,32637.5972250003,18890.7048385460,10745.4190392617,6027.08774087637,3347.68128899872,1849.91729053591,762.440532842799,316.122184882362,178.677962581485,102.630494072788,60.1487403935279,36.1154136451968,22.3061238208926,14.2231968854141,9.38794209944558,6.42245219329092,4.55143171240935,2.86770714365648,1.94391472078553,1.53340128966890,1.22854535752592,0.995686519299282,0.813464067522531,0.667980204949797,0.549990519092996,0.453139677421784,0.372988717912688,0.306291720471319,0.224260383567604,0.165595679425245,0.133555774784540,0.107141751856610,0.0855302001580482,0.0680228700146107,0.0539779144787877,0.0428324540673738,0.0340689075424029,0.0272133115778634,0.0218343711252463,0.0155217934474271,0.0109346539633095,0.00822868421062388,0.00584965541127513,0.00385886035075166,0.00232845695263726,0.00127232551430926,0.000626251197797036,0.000277208360598522,0.000110397466563871,0.0000395798188360959,0.00000760261519740105,9.49429039886151e-7,2.15773276491395e-7,4.28361799511795e-8,7.33689263300493e-9,1.06812734405028e-9,1.30067758026674e-10,1.32610757235812e-11,1.08613725693399e-12,7.34708441825913e-14,8.01196603076191e-15,2.92712569585953e-15,9.62282716451073e-17,0,0,0,7.79966320667010e-23,0,3.48036754063803e-25'''),
320:(0.000015253260877550325,1.0715351883028994,1.5005182545690112e-10,17.61183942070332,877560814.2531368,-0.00810603014148481,
'''6.50164131613495e-15,3.89216149647155e-14,2.53395566641218e-13,9.13755621456808e-12,9.74912755090447e-11,1.86373639011258e-9,1.35357395209542e-8,5.70964604200464e-8,5.25911011884654e-7,0.00000252077333077956,0.0000151938001815391,0.0000562177411088663,0.000241572266639497,0.000719187998042160,0.00235788287135752,0.00585758967890091,0.0155100718000084,0.0333047720600962,0.0557677992280325,0.121066569643294,0.228800313765155,0.464390644233927,0.860306856153328,1.82123109703144,3.69156034006790,9.76396679371297,25.5696419697877,107.980117131500,423.405833436328,1351.35124551708,12294.4517998694,75524.6997987034,901558.380209074,5719392.91832259,54612623.9428277,258211742.315574,1260951064.60826,3406673948.82170,7042412397.75002,9789269554.88670,10210338269.4927,8588439318.98449,6148958390.86101,3835093052.95870,2325724849.39301,1431622188.65569,932981478.223123,643410381.284633,466744036.752,343626314.310490,256335491.573044,209033257.278117,150983423.824896,108587140.901261,76038679.6184656,53263442.7527046,36596985.1974647,25270475.5774765,17189698.8688673,11761843.7646515,7920245.48831159,5357012.89127919,4094657.21290232,2699604.83237203,1779089.07711689,1145397.21877013,733844.076896360,456868.100830985,281494.953739435,167893.664421647,98430.8405057613,55762.1068706887,30924.9758072774,20541.5569278730,10913.6904964944,5700.07008781781,2945.69605118198,1510.40128671955,777.891128791799,402.001597086505,211.183515290362,112.682667490440,61.7439830738375,34.8304469244932,24.2602037762630,14.5607904010196,9.14781241081574,5.98883875121900,4.11677405085656,2.92824637320430,2.17063513681663,1.64213607513542,1.27947518249420,1.00388288544460,0.802992303753624,0.692435711469428,0.553431320215728,0.447069495052379,0.357104654041510,0.287235890044,0.227556692989583,0.181143212666036,0.141714156678159,0.111326725617765,0.0859535020645363,0.0667636392564478,0.0561162880737024,0.0430119846731950,0.0333878154410370,0.0257919708978495,0.0202111676507273,0.0156371236339048,0.0120456191218030,0.00879431603143199,0.006124530825798,0.00382301584571110,0.00217237512962343,0.00138105590163676,0.000628761229326432,0.000248192297218991,0.0000879844734375806,0.0000259980439529947,0.00000722291076674258,0.00000156774305572957,3.35909123675293e-7,5.06329994211040e-8,8.04351150024181e-9,7.68937291432e-10,1.56297384452077e-10,1.54767884832538e-11,6.83700074638533e-13,4.60315308289076e-14,1.93583798256500e-15,2.80171316235817e-16,0,2.98273643953e-20,1.47340195939254e-21,0,0,7.90191681717920e-27,1.00528539424640e-27''')}

def get_psd(mesh):
    e0,ratio,tau,v0,number,ref,s=DATA[mesh]
    edges=e0*ratio**np.rint(np.linspace(0,mesh,221)); f=np.zeros(220); f[70:211]=np.fromstring(s,sep=',')
    assert abs((np.diff(edges**3)/(6*np.pi**2))@f/number-1)<2e-10
    return edges,f,tau,v0,ref
```
<!-- END FILE -->

### reproduce.sh

<!-- FILE: reproduce.sh -->
```bash
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# Install exact dependencies into an isolated environment if needed:
# python -m venv .venv; source .venv/bin/activate; python -m pip install -r requirements.txt
mkdir -p rerun
python -m pytest -q test_metric_research.py
python direct_pole_metric_pilot.py --k .03 --radial 6 --angles 96 --rtol 3e-8 --free-radiation .409 --evolve-phi --stop-at .1 --output rerun/reference-metric.json
python direct_pole_metric_acoustic.py --vev-factor 10 --bias .0032765520306085277 --zmatch 100 --k .03 --angles 200 --radial 6 --rtol 5e-8 --output rerun/new-acoustic.json
python direct_pole_bias_elastic.py --baseline --mesh 240 --slices 1 --output rerun/elastic-base.json
python direct_pole_bias_elastic.py --mesh 240 --slices 4 --nodes 12 --rtol 3e-7 --output rerun/elastic-biased240.json
python direct_pole_bias_elastic.py --mesh 320 --slices 1 --output rerun/elastic-biased320.json
python direct_pole_bias_local_modes.py --output rerun/local-modes.json
python direct_pole_bias_early_bound.py --output rerun/early-bound.json
python direct_pole_bias_bao.py --output rerun/bao.json
```
<!-- END FILE -->

### requirements.txt

<!-- FILE: requirements.txt -->
```text
numpy==2.3.5
scipy==1.17.0
numba==0.65.1
pytest==9.0.2
```
<!-- END FILE -->

### seed-double.json

<!-- FILE: seed-double.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "seed": 2e-08,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 1.4997928337700945e-09,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4088035686345,
  "match_eta_delta_n_per_psi": -1.123519554436656,
  "z_at_stop": 3.4226621878905314,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.09999999999999945,
  "max_matter_fraction": 0.03111816936356059,
  "max_phi": 3.80294404285835e-06,
  "max_delta_x": 0.023332324099811196,
  "max_constraint_derivative_rel": 8.703314597618916e-05,
  "rms_constraint_derivative_rel": 6.545905720765539e-06,
  "early_nfev": 8234,
  "late_nfev": 6023,
  "background_today_x": 1.0140956462376973,
  "background_max_x": 1.301819098222254,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### seed-half.json

<!-- FILE: seed-half.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "seed": 5e-09,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 3.749482084411928e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4088035686345,
  "match_eta_delta_n_per_psi": -1.123519554436656,
  "z_at_stop": 3.3598719295742296,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.09999999999999935,
  "max_matter_fraction": 0.03155186624785635,
  "max_phi": 3.830753199293453e-06,
  "max_delta_x": 0.02199141601286066,
  "max_constraint_derivative_rel": 7.926126368931971e-05,
  "rms_constraint_derivative_rel": 6.304480742858167e-06,
  "early_nfev": 8234,
  "late_nfev": 6176,
  "background_today_x": 1.0253692342617273,
  "background_max_x": 1.3019854335083367,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### shape--0.1-k03.json

<!-- FILE: shape--0.1-k03.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "seed": 1e-08,
  "shape": -0.1,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4088035686345,
  "match_eta_delta_n_per_psi": -1.123519554436656,
  "z_at_stop": 3.3461197427620055,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.1000000000000001,
  "max_matter_fraction": 0.03165189103299151,
  "max_phi": 3.7725852560959138e-06,
  "max_delta_x": 0.02804978482706251,
  "max_constraint_derivative_rel": 2.848973850949564e-05,
  "rms_constraint_derivative_rel": 2.1709791998850246e-06,
  "early_nfev": 8234,
  "late_nfev": 6161,
  "background_today_x": 1.003521354079112,
  "background_max_x": 1.3412714179244782,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### shape--0.14-k03.json

<!-- FILE: shape--0.14-k03.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "seed": 1e-08,
  "shape": -0.14,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4088035686345,
  "match_eta_delta_n_per_psi": -1.123519554436656,
  "z_at_stop": 3.3183177751531066,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.09999999999999984,
  "max_matter_fraction": 0.031850885003802516,
  "max_phi": 3.791481914309556e-06,
  "max_delta_x": 0.022677007693098855,
  "max_constraint_derivative_rel": 4.629744703606218e-05,
  "rms_constraint_derivative_rel": 3.234316420681474e-06,
  "early_nfev": 8234,
  "late_nfev": 6248,
  "background_today_x": 0.9992064377602506,
  "background_max_x": 1.3630220208513277,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### shape-0.1-k03.json

<!-- FILE: shape-0.1-k03.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "seed": 1e-08,
  "shape": 0.1,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4088035686345,
  "match_eta_delta_n_per_psi": -1.123519554436656,
  "z_at_stop": 3.4275925558871543,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.10000000000000034,
  "max_matter_fraction": 0.03108173521110416,
  "max_phi": 3.87206327449011e-06,
  "max_delta_x": 0.029773639999693534,
  "max_constraint_derivative_rel": 3.222472868334931e-05,
  "rms_constraint_derivative_rel": 2.8463471562294935e-06,
  "early_nfev": 8234,
  "late_nfev": 6203,
  "background_today_x": 0.991799676517834,
  "background_max_x": 1.2749895740728854,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### shape-0.2-k03.json

<!-- FILE: shape-0.2-k03.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "seed": 1e-08,
  "shape": 0.2,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4088035686345,
  "match_eta_delta_n_per_psi": -1.123519554436656,
  "z_at_stop": 3.4539488227787687,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.10000000000000041,
  "max_matter_fraction": 0.03090154588644654,
  "max_phi": 3.880943731947848e-06,
  "max_delta_x": 0.04336816663883203,
  "max_constraint_derivative_rel": 9.336624109250267e-05,
  "rms_constraint_derivative_rel": 7.2802042559016395e-06,
  "early_nfev": 8234,
  "late_nfev": 6227,
  "background_today_x": 0.985299875709109,
  "background_max_x": 1.2548536333606537,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### shape-screen.json

<!-- FILE: shape-screen.json -->
```json
[
  {
    "shape": -0.4,
    "finite_field_potential_screen": false,
    "min_U_over_X": -0.06382242559949684
  },
  {
    "shape": -0.3,
    "finite_field_potential_screen": true,
    "z_half": 4.353173233556498,
    "max_extra": 0.01666081019081464,
    "min_curvature_after_half": -20891.94252945557,
    "positive_local_root_integral": 6.746545485151842,
    "activation_integral": 5.351371388381981e+32,
    "cosmic_to_activation": 85.10545839777812,
    "max_x": 1.561196632953058,
    "min_x_after_half": 0.33349269488629774,
    "today_x": 0.996210658947026,
    "retuned_cosmic_proxy": 0.014196979873332372
  },
  {
    "shape": -0.2,
    "finite_field_potential_screen": true,
    "z_half": 4.350428724866024,
    "max_extra": 0.017240123176980118,
    "min_curvature_after_half": -19586.12547116548,
    "positive_local_root_integral": 8.662209386798809,
    "activation_integral": 5.075581613198397e+32,
    "cosmic_to_activation": 64.14004380346715,
    "max_x": 1.407421472693072,
    "min_x_after_half": 0.32193569485992807,
    "today_x": 1.0324951060016196,
    "retuned_cosmic_proxy": 0.01069960644235544
  },
  {
    "shape": -0.14,
    "finite_field_potential_screen": true,
    "z_half": 4.349056998219829,
    "max_extra": 0.01758785169466212,
    "min_curvature_after_half": -18666.16153202656,
    "positive_local_root_integral": 9.737555381635753,
    "activation_integral": 5.012973851732103e+32,
    "cosmic_to_activation": 58.97052249585672,
    "max_x": 1.363003074932408,
    "min_x_after_half": 0.319661024413116,
    "today_x": 0.9992064377602506,
    "retuned_cosmic_proxy": 0.009837245891803209
  },
  {
    "shape": -0.1,
    "finite_field_potential_screen": true,
    "z_half": 4.3476856232527386,
    "max_extra": 0.017835353549930703,
    "min_curvature_after_half": -18204.50123540569,
    "positive_local_root_integral": 10.553550616441807,
    "activation_integral": 4.989625865195312e+32,
    "cosmic_to_activation": 56.54885307346978,
    "max_x": 1.3413009099865518,
    "min_x_after_half": 0.3189275338382211,
    "today_x": 1.003521354079112,
    "retuned_cosmic_proxy": 0.009433271896518445
  },
  {
    "shape": 0.0,
    "finite_field_potential_screen": true,
    "z_half": 4.344943927995254,
    "max_extra": 0.018486238804454765,
    "min_curvature_after_half": -17452.394500529168,
    "positive_local_root_integral": 13.306105793786578,
    "activation_integral": 4.958353227179941e+32,
    "cosmic_to_activation": 52.43035455726708,
    "max_x": 1.3017822162984705,
    "min_x_after_half": 0.31885561955769715,
    "today_x": 0.9770158387549878,
    "retuned_cosmic_proxy": 0.008746239106334854
  },
  {
    "shape": 0.1,
    "finite_field_potential_screen": true,
    "z_half": 4.340834020449015,
    "max_extra": 0.01913217314621902,
    "min_curvature_after_half": -18740.95815569216,
    "positive_local_root_integral": 17.19106652523951,
    "activation_integral": 4.936272841892157e+32,
    "cosmic_to_activation": 49.944625595481156,
    "max_x": 1.2750148659120346,
    "min_x_after_half": 0.32043293549349555,
    "today_x": 0.991799676517834,
    "retuned_cosmic_proxy": 0.008331578934056698
  },
  {
    "shape": 0.2,
    "finite_field_potential_screen": true,
    "z_half": 4.338095837927363,
    "max_extra": 0.019732575043419198,
    "min_curvature_after_half": -23421.844475903083,
    "positive_local_root_integral": 21.534891917095692,
    "activation_integral": 4.931520532692382e+32,
    "cosmic_to_activation": 48.206684333637924,
    "max_x": 1.2548724627747154,
    "min_x_after_half": 0.32326055199567294,
    "today_x": 0.985299875709109,
    "retuned_cosmic_proxy": 0.008041661958342878
  },
  {
    "shape": 0.3,
    "finite_field_potential_screen": true,
    "z_half": 4.333991196108401,
    "max_extra": 0.020268864267429184,
    "min_curvature_after_half": -30175.30105283628,
    "positive_local_root_integral": 25.66206734017962,
    "activation_integral": 4.927137458246844e+32,
    "cosmic_to_activation": 47.03578402015731,
    "max_x": 1.2390943813269666,
    "min_x_after_half": 0.3272492004364436,
    "today_x": 1.0183212819937026,
    "retuned_cosmic_proxy": 0.007846336670198998
  },
  {
    "shape": 0.4,
    "finite_field_potential_screen": true,
    "z_half": 4.329889710486887,
    "max_extra": 0.020730844054058593,
    "min_curvature_after_half": -37508.73739101137,
    "positive_local_root_integral": 29.163036292830565,
    "activation_integral": 4.925089393322681e+32,
    "cosmic_to_activation": 46.230363433218486,
    "max_x": 1.2262312361299195,
    "min_x_after_half": 0.33249904668289454,
    "today_x": 1.011898739097619,
    "retuned_cosmic_proxy": 0.00771197936718215
  }
]
```
<!-- END FILE -->

### slow-1.json

<!-- FILE: slow-1.json -->
```json
[
  {
    "vev_factor": 1.25,
    "tuned": true,
    "seed": 3.7448034398118554e-09,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.018235864909957027,
    "positive_local_root_integral": 10.119561918819562,
    "minimum_curvature_after_half": -11023.850405730998,
    "late_contact_proxy": 0.008862519389152444,
    "early_contact_proxy_100_to_1e5": 2.3045158078559035e-26,
    "early_representative_pair_proxy": 2.2856162123308828e-26,
    "early_x_z1e5": 3.744413444653837e-09,
    "max_x": 3.744413444653837e-09,
    "today_x": -1.0065532498761696,
    "activation_integral": 4.957079946100083e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 1.5,
    "tuned": true,
    "seed": 1.883896789802231e-07,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.017994982001772394,
    "positive_local_root_integral": 8.020948910350338,
    "minimum_curvature_after_half": -7510.633133394172,
    "late_contact_proxy": 0.008965455469665968,
    "early_contact_proxy_100_to_1e5": 1.475760129003028e-19,
    "early_representative_pair_proxy": 1.463663534038262e-19,
    "early_x_z1e5": 1.883760542172953e-07,
    "max_x": 1.2983209329082603,
    "today_x": 0.9977584510573932,
    "activation_integral": 4.9617567272499934e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  }
]
```
<!-- END FILE -->

### slow-2.json

<!-- FILE: slow-2.json -->
```json
[
  {
    "vev_factor": 2.0,
    "tuned": true,
    "seed": 1.9397784187650756e-05,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.017547452200784686,
    "positive_local_root_integral": 5.440290075212942,
    "minimum_curvature_after_half": -4097.4000420926695,
    "late_contact_proxy": 0.009175369958161294,
    "early_contact_proxy_100_to_1e5": 1.660443256149645e-11,
    "early_representative_pair_proxy": 1.646840132778569e-11,
    "early_x_z1e5": 1.939699505388466e-05,
    "max_x": 1.939699505388466e-05,
    "today_x": -0.9669978931982461,
    "activation_integral": 4.960749798078133e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 3.0,
    "tuned": true,
    "seed": 0.0023935739573714174,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.016756836025883848,
    "positive_local_root_integral": 2.9444289975524724,
    "minimum_curvature_after_half": -1731.1814046014956,
    "late_contact_proxy": 0.009549145854671,
    "early_contact_proxy_100_to_1e5": 0.0038517945149076603,
    "early_representative_pair_proxy": 0.0038202514154730732,
    "early_x_z1e5": 0.002393530679546539,
    "max_x": 1.2879733450613933,
    "today_x": 1.0069236438567974,
    "activation_integral": 4.960633753919806e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  }
]
```
<!-- END FILE -->

### slow-3.json

<!-- FILE: slow-3.json -->
```json
[
  {
    "vev_factor": 4.0,
    "tuned": true,
    "seed": 0.013242740446847728,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.016064630862038614,
    "positive_local_root_integral": 1.7564778275003479,
    "minimum_curvature_after_half": -914.3447058782585,
    "late_contact_proxy": 0.00984312145993535,
    "early_contact_proxy_100_to_1e5": 3.59589584379369,
    "early_representative_pair_proxy": 3.56645259232538,
    "early_x_z1e5": 0.013242605762128616,
    "max_x": 0.013242605762128616,
    "today_x": -1.0240745767545896,
    "activation_integral": 4.979739916185609e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 6.0,
    "tuned": true,
    "seed": 0.1942185364441907,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.014851082440954448,
    "positive_local_root_integral": 0.7109165312212711,
    "minimum_curvature_after_half": -310.58254550143135,
    "late_contact_proxy": 0.015129748453342168,
    "early_contact_proxy_100_to_1e5": 166413.99499614356,
    "early_representative_pair_proxy": 165051.53755720003,
    "early_x_z1e5": 0.19421765859170906,
    "max_x": 1.267561888825082,
    "today_x": 1.0287391804074124,
    "activation_integral": 4.9789744303545555e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  }
]
```
<!-- END FILE -->

### slow-4.json

<!-- FILE: slow-4.json -->
```json
[
  {
    "vev_factor": 2.25,
    "tuned": true,
    "seed": 0.0006928845827668006,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.017337836412755565,
    "positive_local_root_integral": 4.598352482579945,
    "minimum_curvature_after_half": -3142.8826228895596,
    "late_contact_proxy": 0.009270518036571,
    "early_contact_proxy_100_to_1e5": 2.702779504197751e-05,
    "early_representative_pair_proxy": 2.680640374033134e-05,
    "early_x_z1e5": 0.000692862310962468,
    "max_abs_x": 1.2931648570386305,
    "today_x": -0.98945156957048,
    "activation_integral": 4.9623736594259775e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 2.5,
    "tuned": true,
    "seed": 0.0003732530303159969,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.017139356713777188,
    "positive_local_root_integral": 3.9129793269516195,
    "minimum_curvature_after_half": -2524.335538938604,
    "late_contact_proxy": 0.009347264252610956,
    "early_contact_proxy_100_to_1e5": 2.272021724281845e-06,
    "early_representative_pair_proxy": 2.2534130437570313e-06,
    "early_x_z1e5": 0.0003732433121635388,
    "max_abs_x": 1.291448646492892,
    "today_x": 0.9891751207870625,
    "activation_integral": 4.971929115558076e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  }
]
```
<!-- END FILE -->

### slow-5.json

<!-- FILE: slow-5.json -->
```json
[
  {
    "vev_factor": 2.75,
    "tuned": true,
    "seed": 0.00075405672727723,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.016946958371745077,
    "positive_local_root_integral": 3.401658413316363,
    "minimum_curvature_after_half": -2076.873741989655,
    "late_contact_proxy": 0.009440536040966922,
    "early_contact_proxy_100_to_1e5": 3.785723338934946e-05,
    "early_representative_pair_proxy": 3.75471939312841e-05,
    "early_x_z1e5": 0.0007540405017227031,
    "max_abs_x": 1.289666015516822,
    "today_x": 1.028448074162698,
    "activation_integral": 4.9709749297170015e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  },
  {
    "vev_factor": 3.25,
    "tuned": true,
    "seed": 0.021490181828789316,
    "z_half": 4.344498797161065,
    "max_extra_zle100": 0.016578594565798327,
    "positive_local_root_integral": 2.5855090027156216,
    "minimum_curvature_after_half": -1458.523303946255,
    "late_contact_proxy": 0.009639032451592807,
    "early_contact_proxy_100_to_1e5": 25.039286135412212,
    "early_representative_pair_proxy": 24.834244454856993,
    "early_x_z1e5": 0.021489850749619465,
    "max_abs_x": 1.2862364587556776,
    "today_x": 0.9727040836410384,
    "activation_integral": 4.9588549682447394e+32,
    "scope": "Fixed potential height and eta bath, tuned seed; contact and representative pair proxies only, no high-redshift production or in-medium scattering solve."
  }
]
```
<!-- END FILE -->

### slow2-k015.json

<!-- FILE: slow2-k015.json -->
```json
{
  "k": 0.015,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 1.9397784187650756e-05,
  "shape": 0.0,
  "vev_factor": 2.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 8.0242384686952e-06,
  "initial_radiation_delta_per_psi": -2.000037024540537,
  "match_delta_m_per_psi": -421.95774328731216,
  "match_eta_delta_n_per_psi": -1.7615186712636821,
  "z_at_stop": 2.431592605269815,
  "reached_eta_threshold": true,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.09999999999999994,
  "max_matter_fraction": 0.015671417559170855,
  "max_phi": 6.36834123195464e-06,
  "max_delta_x": 0.020876534178235744,
  "max_delta_x_over_x_floor": 0.02446059643893275,
  "at_stop": {
    "k_over_aH": 61.94047169594386,
    "eta_number": -0.09999999999999994,
    "eta_number_rest_frame": -0.09881147439128231,
    "eta_background_fraction": 0.00017134085098589825,
    "eta_density_over_total": -1.5158328259783543e-05,
    "delta_x_over_x": 0.011764167973600475
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 7.173190861818551e-08,
  "early_nfev": 9791,
  "late_nfev": 6887,
  "background_today_x": -0.9669978931982461,
  "background_max_x": 8.939887910071359e-07,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### slow25-k012.json

<!-- FILE: slow25-k012.json -->
```json
{
  "k": 0.012,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 0.0003732530303159969,
  "shape": 0.0,
  "vev_factor": 2.5,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 0.00017157885573922251,
  "initial_radiation_delta_per_psi": -2.000023695705944,
  "match_delta_m_per_psi": -301.18133991905955,
  "match_eta_delta_n_per_psi": -1.9573829473265358,
  "z_at_stop": 1.9098000142674538,
  "reached_eta_threshold": true,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.09999999999999999,
  "max_matter_fraction": 0.013088117957685424,
  "max_phi": 7.168525554819441e-06,
  "max_delta_x": 0.014826644691805781,
  "max_delta_x_over_x_floor": 0.016809506524065362,
  "at_stop": {
    "k_over_aH": 52.975758653423725,
    "eta_number": -0.09999999999999999,
    "eta_number_rest_frame": -0.09894556337979322,
    "eta_background_fraction": 0.00015356348238778012,
    "eta_density_over_total": -1.3483956401535e-05,
    "delta_x_over_x": 0.012388598594485527
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 1.3414798295939293e-07,
  "early_nfev": 10706,
  "late_nfev": 7706,
  "background_today_x": 0.9891751207870625,
  "background_max_x": 1.291417309077793,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### slow275-k011.json

<!-- FILE: slow275-k011.json -->
```json
{
  "k": 0.011,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 0.00075405672727723,
  "shape": 0.0,
  "vev_factor": 2.75,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 0.000441602313990983,
  "initial_radiation_delta_per_psi": -2.000019910975133,
  "match_delta_m_per_psi": -262.9038573302812,
  "match_eta_delta_n_per_psi": -1.9749541325989124,
  "z_at_stop": 1.6454048411025681,
  "reached_eta_threshold": true,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.10000000000000003,
  "max_matter_fraction": 0.01249460747457958,
  "max_phi": 7.410274239067635e-06,
  "max_delta_x": 0.011867233667373736,
  "max_delta_x_over_x_floor": 0.013151427807546766,
  "at_stop": {
    "k_over_aH": 50.27492863862153,
    "eta_number": -0.10000000000000003,
    "eta_number_rest_frame": -0.09906762782208114,
    "eta_background_fraction": 0.00014241342162873397,
    "eta_density_over_total": -1.257371470932804e-05,
    "delta_x_over_x": 0.011881577524082403
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 7.739569760635078e-08,
  "early_nfev": 11126,
  "late_nfev": 8186,
  "background_today_x": 1.028448074162698,
  "background_max_x": 1.289710205023639,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### slow3-k01.json

<!-- FILE: slow3-k01.json -->
```json
{
  "k": 0.01,
  "radial_order": 5,
  "angles": 48,
  "rtol": 1e-07,
  "metric": true,
  "evolve_phi": true,
  "compact": true,
  "seed": 0.0023935739573714174,
  "shape": 0.0,
  "vev_factor": 3.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 0.000910193902311648,
  "initial_radiation_delta_per_psi": -2.00001645535135,
  "match_delta_m_per_psi": -225.93090681172757,
  "match_eta_delta_n_per_psi": -2.120738960492483,
  "z_at_stop": 1.360679317739265,
  "reached_eta_threshold": true,
  "chi_threshold": 0.1,
  "reached_chi_threshold": false,
  "max_eta_number_fraction": 0.1,
  "max_matter_fraction": 0.011929827406289659,
  "max_phi": 7.682648464150224e-06,
  "max_delta_x": 0.00875721963938936,
  "max_delta_x_over_x_floor": 0.009541386268411455,
  "at_stop": {
    "k_over_aH": 47.40421858857535,
    "eta_number": -0.1,
    "eta_number_rest_frame": -0.09929624053206647,
    "eta_background_fraction": 0.00012706243402503038,
    "eta_density_over_total": -1.2190997572573246e-05,
    "delta_x_over_x": 0.00416912982758165
  },
  "max_constraint_derivative_rel": null,
  "rms_constraint_derivative_rel": null,
  "max_density_constraint_rel": 1.3961358212278067e-07,
  "early_nfev": 11621,
  "late_nfev": 8744,
  "background_today_x": 1.0069236438567974,
  "background_max_x": 1.2879582950210953,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### small-amplitude-fine.json

<!-- FILE: small-amplitude-fine.json -->
```json
{
  "k": 0.03,
  "radial_order": 6,
  "angles": 96,
  "rtol": 3e-08,
  "metric": true,
  "seed": 1e-08,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-08,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4086214687368,
  "match_eta_delta_n_per_psi": -1.1235193692035854,
  "z_at_stop": 2.6036512072233498,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.09999999999999964,
  "max_matter_fraction": 3.982761937224497e-05,
  "max_phi": 1.660317446116828e-07,
  "max_delta_x": 0.026317600137544823,
  "max_constraint_derivative_rel": 0.0006187794500819836,
  "rms_constraint_derivative_rel": 3.747932853710604e-05,
  "early_nfev": 8849,
  "late_nfev": 9347,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### small-amplitude.json

<!-- FILE: small-amplitude.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "seed": 1e-08,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-08,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4088035686345,
  "match_eta_delta_n_per_psi": -1.123519554436656,
  "z_at_stop": 2.6036511918578427,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.09999999999999967,
  "max_matter_fraction": 3.9827625711291214e-05,
  "max_phi": 1.6603170341429427e-07,
  "max_delta_x": 0.026317597176063034,
  "max_constraint_derivative_rel": 0.004917826168396373,
  "rms_constraint_derivative_rel": 0.00035897708180482687,
  "early_nfev": 8234,
  "late_nfev": 8555,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### test_metric_research.py

<!-- FILE: test_metric_research.py -->
```python
"""Unit and archived-regression tests for the independent research pilot."""
import json
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import zeta
from direct_pole_metric_pilot import potential,radial,X,COEF,M0,T0
from direct_pole_bias_elastic import geometry,pole,flux_jac
from direct_pole_bias_early_bound import matching,integrate
from direct_pole_bias_bao import evaluate,COV
from psd_reference import get_psd

P=Path(__file__).parent

def read(name):return json.loads((P/name).read_text())

@pytest.mark.parametrize('shape,bias',[(0,0),(-.14,0),(.2,0),(0,.0032765520306085277)])
def test_potential_derivatives(shape,bias):
    for x in (-1.3,-.81,-.21,.21,.81,1.3):
        eps=2e-6;u,g,c=potential(x,shape,bias);um,gm,_=potential(x-eps,shape,bias);up,gp,_=potential(x+eps,shape,bias)
        assert abs((up-um)/(2*eps)-g)/X<1e-8
        assert abs((gp-gm)/(2*eps)-c)/X<1e-8
    assert abs(potential(1,shape,bias)[0])<1e-25
    assert abs(potential(1,shape,bias)[1])<1e-25

def test_eta_massless_and_mass_derivatives():
    q,dq,f=radial(32);n=COEF*np.sum(dq*q*q*f);e=COEF*np.sum(dq*q**3*f)
    assert np.isclose(n,zeta(3)/np.pi**2,rtol=1e-12)
    assert np.isclose(e,np.pi**2/30,rtol=1e-12)
    for x in (-.9,.25,1.1):
        temp=.001; y=M0*x/temp;E=np.hypot(q,y);dv=M0*M0*x/(temp*temp*E)
        dr=COEF*temp**4*np.sum(dq*q*q*f*dv);d2=COEF*temp*temp*M0*M0*np.sum(dq*q**4*f/E**3)
        rho=lambda a:COEF*temp**4*np.sum(dq*q*q*f*np.hypot(q,M0*a/temp))
        grad=lambda a:COEF*temp*temp*M0*M0*a*np.sum(dq*q*q*f/np.hypot(q,M0*a/temp))
        h=1e-6;assert np.isclose((rho(x+h)-rho(x-h))/(2*h),dr,rtol=2e-8);assert np.isclose((grad(x+h)-grad(x-h))/(2*h),d2,rtol=2e-7)

@pytest.mark.parametrize('mesh',[240,320])
def test_psd(mesh):
    e,f,tau,v0,ref=get_psd(mesh);assert len(e)==221 and np.all(np.diff(e)>0) and f.min()==0 and f.max()>1e10
    assert np.isclose(ref,-.0081,rtol=.001)

def test_graph_pole_and_conservation():
    e=np.geomspace(.01,5,18);v,K,idx,c,rate,S,l,h=geometry(e);sg=np.array([-1,-1,1,1,1]);assert np.max(abs((K[idx]*c)@sg)/S)<2e-15
    for b in (.01,3.,1e8):
        a=pole(S,l,h,b);assert np.all((a>0)&(a<=1+1e-14))
        for i in np.linspace(0,len(S)-1,20,dtype=int):
            exact=quad(lambda P:(b/(b+S[i]-P*P))**2,l[i],h[i],epsabs=1e-13)[0]/(h[i]-l[i]);assert np.isclose(a[i],exact,rtol=5e-10,atol=1e-13)
    r=rate*pole(S,l,h,3.);eq=1/np.expm1(.3+.2*K);d,J=flux_jac(eq,idx,c,r,v,True);assert v@abs(d)<1e-12*np.sum(r*eq[idx[:,0]]*eq[idx[:,1]])
    f=eq*(1+.2*np.sin(np.arange(len(eq))));d,J=flux_jac(f,idx,c,r,v,True);assert abs(v@d)<1e-12*(v@abs(d));assert abs((v*K)@d)<1e-12*((v*K)@abs(d));u=np.random.default_rng(57).normal(size=len(f))*f;eps=1e-6;fd=(flux_jac(f+eps*u,idx,c,r,v,False)[0]-flux_jac(f-eps*u,idx,c,r,v,False)[0])/(2*eps);assert np.linalg.norm(fd-J@u)/np.linalg.norm(fd)<1e-7

def test_baseline_and_new_elastic_archives():
    b=read('elastic-base240.json');assert abs(b['relative_moments'][1]-b['reference_tail'])<3e-8
    for name in ('elastic-bias240.json','elastic-bias320.json','elastic-bias240-sliced.json'):
        d=read(name);assert -.01<d['relative_moments'][1]<0;assert max(d['max_number_drift'],d['max_energy_drift'])<2e-14;assert d['detailed_balance_residual']<2e-14
    assert abs(read('elastic-bias240.json')['relative_moments'][1]-read('elastic-bias240-sliced.json')['relative_moments'][1])<2e-8

def test_metric_thresholds_and_constraints():
    a=read('evolved-phi-fine.json');assert a['reached_eta_threshold'] and abs(a['z_at_stop']-3.39075)<2e-4 and a['max_density_constraint_rel']<2e-8
    for p in P.glob('bias10-early-k*.json'):
        d=json.loads(p.read_text());assert not(d['reached_eta_threshold'] or d['reached_chi_threshold']);assert d['max_density_constraint_rel']<1e-6;assert d['early_max_eta_m_over_T']<2e-5
    a=read('growth-k03-full.json');assert a['max_eta']<.01 and a['max_density_constraint_rel']<1e-8
    rows=read('bias-local-modes.json');assert rows[0]['max_frozen_positive_integral_post_half']>13;assert max(x['max_frozen_positive_integral_post_half'] for x in rows[1:])<.23;assert max(x['max_fast_subhorizon_integral'] for x in rows[1:])==0

def test_early_bound_and_matching():
    a=matching();assert a['field_cap_over_Mpl']==.1 and a['lambda_A']<2e-21 and a['constant_A4_finite_loop_ratio']<.7
    d=integrate(1e-8,240,6000);assert abs(d['combined_amplitude_rge101_L1_bound']-read('bias-early-bound.json')['rows'][1]['combined_amplitude_rge101_L1_bound'])<3e-7
    assert d['combined_amplitude_rge101_L1_bound']<.001

def test_bao_baseline_and_covariance():
    assert np.all(np.linalg.eigvalsh(COV)>0);d=evaluate(.315,'lcdm');assert abs(d['chi2']-14.308052789821941)<1e-6;assert abs(d['ruler_Mpc']-148.63535474133934)<1e-5
    out=read('bias-bao.json')['results'];assert abs(out['biased']['profiled_omega_m']['chi2']-out['lcdm']['profiled_omega_m']['chi2'])<.1
```
<!-- END FILE -->

### threshold-point1-new.json

<!-- FILE: threshold-point1-new.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 64,
  "rtol": 1e-07,
  "metric": true,
  "seed": 1e-08,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "free_standard_radiation_fraction": 0.0,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.0001276663934515,
  "match_delta_m_per_psi": -1000.0911987884457,
  "match_eta_delta_n_per_psi": -1.0498416727301105,
  "z_at_stop": 3.383480490238939,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.10000000000000006,
  "max_matter_fraction": 0.029325979378868158,
  "max_phi": 3.5657398325544344e-06,
  "max_delta_x": 0.02446129650965613,
  "max_constraint_derivative_rel": 3.231000561477097e-05,
  "rms_constraint_derivative_rel": 2.6639147260412317e-06,
  "early_nfev": 8234,
  "late_nfev": 6170,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### threshold-point1.json

<!-- FILE: threshold-point1.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 64,
  "rtol": 1e-07,
  "metric": true,
  "seed": 1e-08,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-05,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.0001276663934515,
  "match_delta_m_per_psi": -1000.0911994274165,
  "match_eta_delta_n_per_psi": -1.0498416733452398,
  "z_at_stop": 3.383480486957721,
  "nonlinear_eta": true,
  "max_eta_number_fraction": 0.09999999999999981,
  "max_matter_fraction": 0.029325979419042848,
  "max_phi": 3.565739867053653e-06,
  "max_delta_x": 0.024461304705956706,
  "max_constraint_derivative_rel": 4.785544100016013e-05,
  "rms_constraint_derivative_rel": 3.67950085589524e-06,
  "early_nfev": 8234,
  "late_nfev": 6158,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Ideal standard radiation fluid, pressureless standard matter, no recombination, standard neutrino hierarchy or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### tiny-amplitude-fine.json

<!-- FILE: tiny-amplitude-fine.json -->
```json
{
  "k": 0.03,
  "radial_order": 6,
  "angles": 96,
  "rtol": 3e-08,
  "metric": true,
  "seed": 1e-08,
  "shape": 0.0,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-10,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4086214687368,
  "match_eta_delta_n_per_psi": -1.1235193692035854,
  "z_at_stop": 2.049008223337972,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.09999999999999987,
  "max_matter_fraction": 3.7434040552288556e-06,
  "max_phi": 1.5376989396438848e-07,
  "max_delta_x": 0.014775259040164359,
  "at_stop": {
    "k_over_aH": 129.96281070075426,
    "eta_number": -0.09999999999999987,
    "eta_number_rest_frame": -0.09954483663582922,
    "eta_background_fraction": 0.00012276742078269975,
    "eta_density_over_total": -1.0169217261668595e-05,
    "delta_x_over_x": 0.017575174174265067
  },
  "max_constraint_derivative_rel": 0.004781522323095467,
  "rms_constraint_derivative_rel": 0.0003583035061228184,
  "early_nfev": 8849,
  "late_nfev": 12272,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### tiny-amplitude.json

<!-- FILE: tiny-amplitude.json -->
```json
{
  "k": 0.03,
  "radial_order": 5,
  "angles": 80,
  "rtol": 1e-07,
  "metric": true,
  "seed": 1e-08,
  "zmatch": 12.0,
  "zstart": 1000000.0,
  "psi_initial": 1e-10,
  "free_standard_radiation_fraction": 0.409,
  "eta_number_threshold": 0.1,
  "early_max_eta_m_over_T": 7.498964168829178e-10,
  "initial_radiation_delta_per_psi": -2.000148098162149,
  "match_delta_m_per_psi": -1070.4088035686345,
  "match_eta_delta_n_per_psi": -1.123519554436656,
  "z_at_stop": 2.049008082567796,
  "reached_eta_threshold": true,
  "max_eta_number_fraction": 0.09999999999999988,
  "max_matter_fraction": 3.743404464434018e-06,
  "max_phi": 1.537696037065621e-07,
  "max_delta_x": 0.014775241224624395,
  "max_constraint_derivative_rel": 0.02653476386746583,
  "rms_constraint_derivative_rel": 0.0015477743047628847,
  "early_nfev": 8234,
  "late_nfev": 11165,
  "background_today_x": 0.9770158387549878,
  "background_max_x": 1.3018933502539558,
  "scope": "Pressureless standard matter, radiation fluid plus optional massless collisionless radiation, no recombination, baryon-photon drag or nonlinear evolution. Initial adiabatic-like mode normalized to early Newtonian psi. Stop at eta number contrast threshold."
}
```
<!-- END FILE -->

### SHA256SUMS

<!-- FILE: SHA256SUMS -->
```text
84f3e6f42cf3e4672a3611c5c51eeaecabe3c5acec00e1d30640e2e1e79fb058  acoustic.json
6e32d1f82567d8fbf53066afbfa1cb70adcc352e45720b1cc8cd23600e7f4ccb  band-.003.json
efe99f93c6ed17d69a1d479ef40b51b1d743c5f47b6c458136a7b39e66b5afae  band-.006.json
d0c918bce834e1818cfe98c9a8927d8503f226399903bc4e0bcfd6bf15dbecbb  band-.01.json
fa7bd5219aaa506e8b2a8e0b1263e706d053f27736d718d3192dcaa33eb5976c  band-.015.json
32c09c7e3b6212f5ae1c282083df58e088cd1e9754351a0f90df1f0b4e6db2db  band-.02.json
925e39bf43ed4a428fd8715dba86fd03fde91b059c5f823b16a06e1f90202aeb  band-.03.json
5a2c939e05423409cdefc726d4c8995a2a35ee183493874f1a4e64cc16cdf89f  band-.04.json
f220df17f4b8de18f26aac14b1d426070a2e29610ce752b8facd0bdb4cb6783b  band-.05.json
01a09c9a6783bd464e3cd3b072cfccc80d3795ee61f5ac96bb4df123e3981e09  band-.06.json
aa88c76640a65301aa4feea489e6f9ff5dd74ca729b601c1587977b4cfb5ba5f  bias-1.json
cd409c41c987bcbbf383b5928ff2e8b2ab1191df4070b022af338093134f8f2e  bias-2.json
d33eb7280b706e5280235b491334d958c1e19d7c0583f726c7914cd9e514e758  bias-3.json
bd542ce6c7877c6765c66a8e479d91fb0f8e11bed71f036b0342f565e1f66f86  bias-bao.json
ef51023745f286c3f7a376d946f8fccf7c033d3b01a9bf5a0bc3d2c3c3567c50  bias-early-bound.json
95fab31dd5242a1eb50b8a3cef8ef5f1ab256481dc397b62e51baac67dcef226  bias-late32.json
ad3763c794d7853920eb6953b67f6859daf9a1aec0f5679dbd0e0da11e6bc202  bias-late35.json
476d4ab675269118e52ae95a5db0b082891738556750642b670ff07715bbe5c4  bias-late38.json
11735f8bfb8028d8e3ec82304f224d8052414aff9c6c0c9d673c9cad148f3630  bias-local-modes.json
16eb1f22bb000896c476c365fdc5de8674e6042b92d600a7662582fea2039fb5  bias-seed-wall.json
997bbee3d2e92d6fdc3aa7264c243b945513581e50f8424bb5f8b76e14febcc4  bias-window-0.json
38fb2a4df58c5c36e810f630d6e17400efc44ca0bc2e363b1f54735eca497bdb  bias10-acoustic-fine-k.03.json
b25ddab00aaa1f51140e31088c66cd7a8aeb3ae2f21c3277e5c7fb11329ee3c6  bias10-acoustic-fine-k.05.json
9a84304a4a2739fe4e805b7cfc171da898e8f6b49e16e0d39f2d1455152a0f32  bias10-acoustic-fine-k.06.json
5752431918cf091dae0c0b08280d50e70a21ec22dd944e11b51bc333882bd748  bias10-acoustic-k03.json
b6b34f393ef0b0d0539d674b1f49e44eee6fe7b05f14c2038cd73d29440f72fa  bias10-chi_iso-k003.json
2e6a6dd0dc6f0c1f982526ac801744d52bbbb3e98d4876c39b75e26cf91c5a59  bias10-early-k.0003.json
ba441e69e43646237e61b44ea99582e33016e954c1a2461340c5071cb9bee1ca  bias10-early-k.001.json
67c576999a56c369ee9d2c65f10b3867c4569cf5583437e087a5879063c0bac1  bias10-early-k.003.json
0f85fb13894582826f8d043b077c5db9c4cc90e1fadb1a8e6b31b6ba445bc42d  bias10-early-k.006.json
efa1f9681506b3e6f7f5abca49fa489c353950f2af567f8a9869482dfd0e207e  bias10-early-k.01.json
15a9e521e695c0351a88ad6f5259f1ae5dba51b400ed0f8ccfb20d6929877b67  bias10-early-k.02.json
0bd42a31c57b7a7161258cddbe1ebc752d88eefa31c37a2167610ad6e142b714  bias10-eta_iso-k003.json
7c7a380550539a9502d8253fff78c5c6a4f9c48049ddb6a1f1f5ee139086becf  bias10-k.001.json
7ba9fe1bdf58c593af4c688dcc9340c2a7eb8d3c2dcb8861cdb7c145f3a62ac3  bias10-k.002.json
ae9883c17bd51dfbb9fd819901441151f338f660a198f802b4cfb9b4d2ee6a09  bias10-k.003.json
d000c54943ad263cd5f982e5d69fc2ea0c35b7cedbb7717bd7b2f645466ac6be  bias10-k.006.json
07c53f344295819ee5e563434de5c1f249231ab0fcdf6f221f846ecc8f973d81  bias10-k.01.json
498ea19d3c254acc6643ccf5e810095820ad713c73ba0a4905919b56513e681c  control-no-late-metric.json
9add24639604f45109e9922e25d964279a586d8d9be8d815bda5afb96b69aa1d  direct_pole_bias_background_checks.py
16f7f54775a39241b7f17f95c14fa1b14e3b2e4d206018527c0c1c306d3ce462  direct_pole_bias_bao.py
e40c0e7968919404cdf83ff6bc4e531f563df3b64f3bbb495182f4bfb5030ebb  direct_pole_bias_early_bound.py
0d73f6a9db0e6f5e63bed0cff2b471f2677900613e82386bba008dd19e99a9e6  direct_pole_bias_elastic.py
361f78c8dfde4b321f89ad7e02b02d122697c835c1bac79e2e5a0a11a604b2ef  direct_pole_bias_local_modes.py
2beea9f18fe384333ede52a038ffc318b2a14b642a2ee28cd3a30b65de3f66ae  direct_pole_bias_parameter_window.py
34c74a2bdebee8b34fbb4268f470948ea63a30b3aa01ee3aaee31dbbc3ddb65b  direct_pole_metric_acoustic.py
6396d68172dbbfdca0748620ee338f30ece28d22170088f0079f44b61439d7b7  direct_pole_metric_pilot.py
f4d73c83de8ca1f62f8738834400ce569a0b3b060e114d1918f0edee100b9773  direct_pole_slow_trigger.py
1ebfe5662c7ad363d7e131d134a5300577b98ad6a992f459cc23e7ad8c209128  direct_pole_trigger_shape.py
4489d28caf2ef81a29f37ce06e3f547558bf5cae56571a62f86e16d79ebfdc67  elastic-base240.json
895588bb3e1490ceb78b9dcccc2ca80c8abe377bf8799b8699cf87b0173f8a98  elastic-bias240-sliced.json
2225e3b1345d8a5fb42122ee5e8051b00c7a781d3673a1290b246c1c86e30ce4  elastic-bias240.json
ccbc5e0ec6b111fb81f5f54556a62fdfbd9dbaf7e970b9ed183318fe0367e9c7  elastic-bias320.json
e1ae55b9ea1457342f63b7bb244d9a776d8581f02183326bf3a652098e8f330e  evolved-phi-fine.json
e72fd4a779a03bd4b4d6c13bbb95df656f288fe101d516f4e994ccf0780c1fcd  evolved-phi.json
70271fcfd81ec107f3d32c72bca7b8c84f53281b0df2a0a61c1f697875620265  free-rad-early.json
775e158f5ec26e2d34a92d0f510575afa4231aa0df12f3904eee62fc294d7ee1  free-rad-fine.json
88b764abb0d9d7af0d5f50b859c2e7b4af58ea686433d0c29e9586b32b9aa9db  free-rad.json
b1cf4bf13dec4ecbde12df2590a93a42e9c6b40544c204f0e39dad117d0102f0  growth-k03-full.json
bd4649687c785e0b3116106cf9c61e77ea2677d1cd01cd8044d5f71042d54f07  growth-k03-lcdm.json
dee6ba9b88489346d04cfdb1a0db4e06df0fd9731f813bb23e478196f2917932  growth-k03-massless_eta.json
a36788bc0404b0141c1533f89528a4df7bf0c5277f2092f4244dfeaf75407761  k02.json
871440446b99274aa3409dd72bb452d5a38dc0057821379c59852764c4080fcf  k04.json
4cb1c0a61e74834e1053ea36bf1630891d9374b266e6a105d238a0c836f77660  pilot-1.json
b8bdb169fcd5a4d09a09388eb07e30dd13b95d4c6069ec7417203daebf347e46  pilot-2.json
53ff3ce4bf70057919d99db138ad4a416f64e5601f1cf70bf8ae65a0691d5583  psd_reference.py
07b23feb9390d907c4a0144512a430d6480b54a5282bfc08407b83353d238455  reproduce.sh
a1f530e13f46bbed2c482046c8fb5ea82f21202b991cdb29d3dd64ddd0cf286f  requirements.txt
11fe4ea2eeea2479a5a6cff9e6fe86f601309e79aef3f6474135039408bfd364  seed-double.json
d82edbff47cc411d50af771979ecb9793b0e23c7dfdf2cfddad7745fe05dd527  seed-half.json
190a8b5b85afc4375280c6495f163db4af35ad0cab7690e168626d20f242877e  shape--0.1-k03.json
d37bcb4132b847f31bc4004f134ee86286134931146a13e3a04d42e9d91c2559  shape--0.14-k03.json
edb0e77020e184eaf0ff4b3c97db0ea806bc453e88ae1bd7086fd63582fd7253  shape-0.1-k03.json
a7eb998414a3e4ef38eadd2a9b152d6f8b21699a810e57f63da765e645827ba2  shape-0.2-k03.json
13faadc716b937d27f927d0fa754721d700e114531340f5ee93e87e42a24c4ff  shape-screen.json
f6bfcc3c3d2b5dcc269ff5f8b69be170768611e42739d03f829d0b2c7ddbba32  slow-1.json
3ac18e9dd04d3978c08348f20c8efb2ea6edebe66bfa59f81775fcf2c52349bb  slow-2.json
24a30841d56306248402b771beafadd68820a6c73c926cf6b97e395496e382b1  slow-3.json
4d23d72d17bc2dbe334303dcbe31a7e320efce8ceab1791fb2325d555dd25881  slow-4.json
ce6a544a5ab004c2e0364aed3f17e605e92e3a970d4054db44a96d15545c3b9b  slow-5.json
bd71eeee58ec15c8c39d806bb52237c0fcde1a6b396477c22516ebf70c45fe7d  slow2-k015.json
002f294147e12da7c5c144b25904cde44a6811f0360629718914d8cc2a2d3ba3  slow25-k012.json
211ef6ea0a9317e8d9dd4dcb6ea7056fd7cb12d96b572eaadad2bade8856c151  slow275-k011.json
16e95c2a81672fc7bbf41169ffcd6c8aefea536e89846222808521ff728225fb  slow3-k01.json
7883fe1e8b18671744dcf0a5d45b3eac12b962bc97cf3104ee780cc29aa06af2  small-amplitude-fine.json
48e0442a7b7f3e4bcea5974205a7cd78130926001b51b98d05d9332b8524925c  small-amplitude.json
9726ed97d9f3cbda32fa2d1d0b15c03c9104acc312a92201eeff53a97e5019bd  test_metric_research.py
6542946e0dc832e92d51e2d7a873653b6a128c568ddcc5368f35e56995d56768  threshold-point1-new.json
7e267c714e6f7eacbd29165b09a37f1672b1b91931e89f69c5ad5b34261a65c4  threshold-point1.json
3550536c691ca76488011e608019e81674c8530f029564227d991f28da8e25ef  tiny-amplitude-fine.json
d16ef0adcf9dde8f81300ac3e58ba027d615e90e16a9f135ff58e4b90ed265ff  tiny-amplitude.json
```
<!-- END FILE -->
