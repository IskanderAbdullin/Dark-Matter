# PLAN — карта исследований Dark-Matter

**Живой маршрут, не замена [GOAL.md](GOAL.md).** Научная цель охватывает космологическое происхождение бозонной тёмной материи и галактические распределения/гало. Нельзя считать всё исследование законченным из-за успеха отдельной модели.

| Если задача | Открыть |
|---|---|
| Определить новый научный этап, критерии и очередность | [Подробный изменяемый план](plans/RESEARCH_ROADMAP_DETAIL.md) |
| Перейти к кинетике, микрофизике и freeze-out | [Кинетическое ТЗ](kinetics/CONTINUOUS-RESEARCH-SPEC.md), [специальная цель](kinetics/research-control/GOAL.md) |
| Продолжить фактический расчёт | [STATE](kinetics/research-control/STATE.json) → [NEXT_ACTION](kinetics/research-control/NEXT_ACTION.json); открыть нужную часть [EVIDENCE](kinetics/research-control/EVIDENCE.md) |
| Проверить FRW/космологическую наблюдаемость | [cosmology/README](cosmology/README.md) |
| Исследовать гало, профили, кривые вращения | [PROJECT_MAP](PROJECT_MAP.md) → ветка галактических исследований |
| Проверить историю и воспроизводимость | [archive/CHRONOLOGY](archive/CHRONOLOGY.md), [archive/README](archive/README.md) |

**Правило:** PLAN/подробные научные планы можно менять при появлении доказательств, сохраняя причину, альтернативы, научные ворота и следы Git. GOAL и разрешения не менять без пользователя. Исследование не прекращается на первом CI или коммите: следующий релевантный шаг выбирать самостоятельно в пределах текущего запуска.

**Научный приоритет P0 после 2026-11-24:** минимальный real-φ Yukawa + одинаковый векторный U1 заряд χ1/χ2 сохраняет стабильного тяжёлого χ2 (20 МэВ). Его симметричный тепловой toy-реликт при y=.001 сильно перепроизводится даже без φ; см. [постоянный отрицательный физический аудит](kinetics/PHI-BOSE-STABLE-CHI2-THERMAL-RELIC-GATE-2026-11-24.md), [условный векторный SM portal](kinetics/PHI-BOSE-VECTOR-PORTAL-PARTIAL-SM-GATE-2026-11-24.md). **Проверять суммарную Ωχ2+Ωφ и реальные начальные условия в полной FD/BE+SM+FRW модели**, не игнорировать стабильного партнёра. Правило остановки только при объективно полной блокировке — в [AGENTS](AGENTS.md).

**Новый gate P0 (2026-11-24):** реальные CLEO R_had(9.996GeV) 3.497 ± статистическая/систематические ошибки измерены (но отличаются на 4MeV от mV10); точный Noether dark U1 ограничивает тяжёлую χ2 асимметрию вне MB; электронный A'10GeV νtransport/H≈1.56e−7 при χ2 химическом x8 тогда как точный парный Yukawa Γ/H≈4.97. **Нельзя** называть прежний toy Ωχ2pair h²≈880 настоящим relic указанного vector SM модели, поскольку общий Tχ2,φ,SM не выведен. Следующий шаг: отдельная dark/SM entropy+FRW, сумма Ωheavy+Ωφ и наблюдения. [Первичный исследовательский отчёт](kinetics/PHI-BOSE-CLEO-HADRONIC-U1-ASYMMETRY-AND-THERMAL-GATES-2026-11-24.md).

**P0 update (2026-11-24, 21:34 UTC):** [Variational entropy-energy/count gate](kinetics/PHI-BOSE-ENTROPY-ENERGY-GATE-2026-11-24.md) shows initial thermal + isolated benchmark requires at least ~0.89 neutrino energy units during BBN independently of later dark thermalization; scalar thermal freezeout Ω≈0.848; exact entropy-only today bound Ω≥0.0837 does not alone exclude. Prioritize actual initial entropy/heat export, lower-entropy nonthermal histories, and time-dependent BBN; count gate conditional on number nondecrease, not proven. No universal exclusion.

**Cold-initial-history P0 (2026-11-24, 21:36 UTC):** [low-reheat lepton injection/one-pass count gate](kinetics/PHI-BOSE-PORTAL-LOW-REHEAT-GATE-2026-11-24.md): at TR=20MeV, even converting all injected heavy counts to φ in one pass yields Ωφ≈1.34e-5; 99% scalar mass fraction needs heavy survival ≤8.69e-8. This is not a full freeze-in result: evaluate physical φ multiplication, backreaction, hadrons and other initial production mechanisms rather than assume a cold sector suffices.

**P0 cold PSD (2026-11-24, 21:39 UTC):** [thermal free-streaming and arbitrary cold-PSD entropy gate](kinetics/PHI-BOSE-THERMAL-STREAM-ENTROPY-GATE-2026-11-24.md). Ordinary all-DM 86eV BE spectrum equality free stream 6.65 Mpc. Chosen 5.7keV FD rms diagnostic requires q_rms<=1.14e-6 eV; ideal Bose entropy upper bound at that moment and abundance, including <=1% heavy energy, s_total<=.07939 cm^-3. Thus initial high-entropy isolated branch also fails this cold-moment condition; actual LSS likelihood and physical nonthermal production remain open.


**P0 nonthermal source (2026-11-24, 21:44 UTC):** [momentum-resolved portal gain gate](kinetics/PHI-BOSE-PORTAL-INITIAL-GAIN-GATE-2026-11-24.md). Tiny MB initial inversion plus spontaneous parent depletion yields reduced gain exponent upper <3.47e-5 on 42 sampled points. No nonlinear or FD claim; prioritize actual FD angular source and depletion over assumed fixed mother bath.


**FD pole source (2026-11-24):** [exact angular quantum occupation correction at vacuum pole](kinetics/PHI-BOSE-PORTAL-FD-POLE-GATE-2026-11-24.md). Eight sampled positive MB pole inversions survive; FD total rate suppressed at T5 by 34%. Finite width matters for differential inversion, so no full FD or burst claim.


**Full FD vacuum BW source (2026-11-24):** [finite-width electron injection and reduced gain](kinetics/PHI-BOSE-PORTAL-FD-BW-GATE-2026-11-24.md). Max reduced exponent 2.77e-5 on 8 samples; negative inversion in one cold .688 GeV point. Next actual stable chi2 depletion and full bath, not more pole-only diagnostics.


**Vacuum stable partner depletion (2026-11-24):** [arbitrary-momentum envelope](kinetics/PHI-BOSE-PORTAL-VACUUM-DEPLETION-ENVELOPE-2026-11-24.md). Exact σvM<=y⁴/(8πM1²). Conditional low-TR MB fixed injection + maximal vacuum scalar annihilation still leaves >=86.1% at TR20MeV. Need actual stimulated modes and additional channels; no full exclusion.


**Low-reheat cold count (2026-11-24):** [kinematic/source-tail envelope](kinetics/PHI-BOSE-PORTAL-COLD-TAIL-BOUND-2026-11-24.md). At TR20MeV single-pass redshift+decay+terminal pair with arbitrary stimulation gives cold Ω<=2.27e-14 vs .060 needed for chosen rms reference. Not a full kinetic/LSS no-go; test redistribution/recycling or alternative initial production.


**Repeated production energy gate (2026-11-24):** [number/energy upper bound](kinetics/PHI-BOSE-PORTAL-RECYCLING-ENERGY-BOUND-2026-11-24.md). Includes fermion conversion, repeated inverse pairs and lepton elastic heat under defined low-TR vacuum network. Arbitrary dark momentum redistribution/stimulation still gives Ωφ<=4.387e-4 at TR20, <=.03844 at TR30. Missing scalar self number-changing and full visible thermal channels remain; prioritize these and distinct initial histories.


**Scalar radiative completion (2026-11-24):** [CW matching and halo/cold field diagnostic](kinetics/PHI-BOSE-SCALAR-QUANTUM-COMPLETION-2026-11-24.md). Renormalized scalar mass and quartic independent. At Q=.1 δm²/mφ²~8.27e4 (scheme-dependent), loop-size positive illustrative λ makes all-DM harmonic extrapolation inconsistent above ~6.97MeV. Need actual interacting background and scalar number-changing; no claim λ fixed.


**Coherent alternative (2026-11-24):** [exact mixed-potential action and Lamé resonance](kinetics/PHI-BOSE-COHERENT-RESONANCE-2026-11-24.md). At TR20, conditional all-DM homogeneous positive-λ examples μ/H>5.86e13, but resonant q below the chosen cold rms diagnostic. Cannot extrapolate to nonlinear final PSD. Distinct coherent initial history not covered by particle-injection bounds.

**P0 update (2026-11-25):** [Nonlinear classical fragmentation](kinetics/PHI-BOSE-QUANTUM-SEED-2026-11-25.md) from conditional finite-band vacuum covariance reaches >10% fluctuations after 77/42 periods, negligible Hdt, with no converged late UV spectrum. [Swept 4→2 comparison](kinetics/PHI-BOSE-NARROW-SWEEP-2026-11-25.md), [chemical/count endpoints](kinetics/PHI-BOSE-COHERENT-ENDPOINTS-2026-11-25.md), [fermion mass crossing](kinetics/PHI-BOSE-COHERENT-FERMION-CROSSING-2026-11-25.md) prevent extrapolating homogeneous branch to a relic. Prioritize physical initial production plus FRW scalar number-changing/redistribution and stable partner abundance, not further fixed-background resonance scans.

**Halo P0 update (2026-11-25):** [Global projected Antonov form](continuation-04/papers/15-global-projected-antonov-critical-halo-2026-11-25.md) proves radial A−B≥A/4 and nonradial sign modulo translation even at critical equilibrium (conditional linear energy statement, not critical nonlinear theorem). [Initial C22 escape](continuation-04/papers/16-initial-bose-collision-escape-2026-11-25.md) and [deep-state flux](continuation-04/papers/17-initial-bose-deep-state-flux-2026-11-25.md) quantify collisional nonstationarity, including independent soft coefficient .185/.367 per Gyr but tiny absolute deep count. Next construct conservative orbit-averaged quantum transport with escaping bound/unbound populations, self-consistent potential and eventually GP exchange; derive actual evolution before claiming condensates or survival. Compare observation only with physical evolved profiles. Both scientific paths remain active; CLASS-ready=false.
