# Этап: 2026-11-22 — независимый entropy maximum и неравновесная эволюция ν/antiν/Z′

- **ID / UTC:** DM-KINETICS-20261122-NU-ENTROPY-LONGTIME; точная длительность не зафиксирована, не выдумывается
- **Направление / статус:** kinetics, аналитическое решение и локально проверенный численный пилот, новый GitHub Scientific CI после push обязателен
- **Ветка / стартовый HEAD:** research/kinetic-equilibrium-audit-20261121, [baa501e5394c](https://github.com/IskanderAbdullin/Dark-Matter/commit/baa501e5394c0fe52604025e2c7e2ebb4f534cd7), итоговый commit — Git-запись этого stage
- **Источники правил:** [GOAL](../../../../GOAL.md), [PLAN](../../../../PLAN.md), [kinetics GOAL](../../../../kinetics/research-control/GOAL.md), [AGENTS](../../../../AGENTS.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json)

## Предварительная гипотеза и независимый эталон

После доказанного отрицательного результата исходной integer ν сетки — два искусственных k=Emax conserved modes, nullity5 вместо трёх — проверить, что исправленный half-offset discrete event graph с nullity3 действительно эволюционирует от **отдельных неравновесных** fν(k), fantiν(k), bV(p) к физически ожидаемому для ЭТОЙ сетки максимуму квантовой энтропии, определяемому ровно C=Nν+Nantiν+2NV, L=Nν−Nantiν, E. Критерии: три constraints независимого root до 2e−11, стационарность FD/BE в оригинальном микроскопическом event source до 3e−13, рост S, соответствие long-time mode occupations и S с указанием абсолютных невязок, сохранение нерешённых q/dt/cutoff. Альтернативы: extra approximate slow charges, ошибка формы энтропии/дегенерации, первая схема ошибается равновесием, long-time splitting затягивает релаксацию.

## Метод, физические ограничения и первоисточники

На одинаковом энергетическом half-offset пространстве k_i=(i+1/2)h, E_j=jh с exact cell angle overlap и clipped vector density-of-states из [предыдущей геометрической проверки](2026-11-22-vector-neutrino-radial-endpoints.md), максимизировать дискретную квантовую энтропию с линейными C,L,E, не используя исходники collision/time в constraint root. Внутреннее решение: T>0, μν/T, μantiν/T; FD fν=expit(μν/T−k/T), fanti=expit(μanti/T−k/T), BE bV=1/expm1((E−μν−μanti)/T), µV=µν+µanti<min E. Отдельно проверить fixed point quantum collision и долгую динамику с теми же начальным C,L,E. Quantum framework: [Froustey/Pitrou/Volpe](https://arxiv.org/abs/2008.01074) — физический контекст, ни код, ни лицензируемые данные не копировали.

**Строгий scope:** только broken-phase Z′↔ν+antiν, gV=3, gν=ganti=2, mV=.4, Emax4, Γ0=1 условная единица. Начальные спектры синтетические fν=expit(−k+.15e^−k), fanti=expit(−k−.1e^−k/2), bV=.75BE(E). Нет физического χ SM Γelastic, actual lepton weak ν↔eγ thermalization, Z′→dark channels, χ1χ2φ Bose, H/FRW или 86-eV DM model.

## Данные и воспроизводимость

- [Независимый thermodynamic root/time сравнение](../../../../kinetics/vector_neutrino_constrained_equilibrium.py)
- [Регрессионные тесты](../../../../tests/test_vector_neutrino_constrained_equilibrium.py)
- [Полный научный аудит](../../../../kinetics/VECTOR-NEUTRINO-CONSTRAINED-EQUILIBRIUM-AUDIT.md)
- [Версионированный численный JSON](../../../../kinetics/results/vector-neutrino-constrained-equilibrium-2026-11-22.json)
- CLI: python -m pytest -q tests/test_vector_neutrino_constrained_equilibrium.py
- CLI: python -m kinetics.vector_neutrino_constrained_equilibrium --h 0.2 --tau 500 --dt 1

Python NumPy/SciPy, без случайности, seed N/A. Сетки h=.4,.2,.1 для равновесия; τ=100/500/1000, dt=2/1/.5/.25 для временного negative-control. Входы, единицы и фактические численные данные находятся в указанном JSON, CI URL будет добавлен после фактического завершения.

## Фактический результат и независимость

Получено независимо от collision и времени:
h=.4: T=1.03058018, μν/T=−.10959360, μanti/T=−.16466833, μV/T=−.27426193, ΔSmax=.00108105814.
h=.2: T=1.02842113, μν/T=−.10581032, μanti/T=−.16095289, μV/T=−.26676321, ΔSmax=.00113072711.
h=.1: T=1.02762325, μν/T=−.10440679, μanti/T=−.15957578, μV/T=−.26398257, ΔSmax=.00115047138.
Все constrained C/L/E невязки ≤4.5e−16. Поскольку L≠0, μν≠μanti; это реальная редуцированная дискретная charge, не произвольное FD-closure. В microscopic event sweep эти состояния остаются fixed points ~1e−13.

Подлинная **неравновесная эволюция без thermal reset** в том же графе при h=.2 и Γ0τ=500, dt1 даёт max |b/f−FD/BE_eq|≈.00103721 и оставшийся положительный entropy gap≈1.8630e−6; при Γ0τ=1000, dt1 разница≈.00033923, gap≈3.9929e−7. Дискретные C/L/E сохраняются до машинной ошибки, S растёт. Это свидетельство *ограниченной дискретной* релаксации, не континуального физического eigenvalue.

## Отрицательные данные и научная граница

При h=.2, τ=100 результат низкоэнергетического vector bV для dt=2,1,.5,.25: .75305223,.74397501,.74667793,.74874261. Сходимость **не монотонна**, time-order/first-order splitting bias остаётся. Стационарные T и µ меняются по h, а предшествующая независимая weak ↔ lattice rate ошибка на h=.05 number ~−0.142%, cutoff tails и threshold неисправлены. Zeroν mode, Γ/H, H/SM, 86 eV φ/Ω/BEC/CLASS отсутствуют. Нет оснований утверждать подтверждение физической бозонной DM модели.

## Решение по плану и точный NEXT_ACTION

[PLAN](../../../../PLAN.md) сохраняет цель; отдельная независимая constrained-equilibrium проверка становится важным численным gate перед long-time cosmology. Полученный дискретный максимум нельзя выдавать за thermal cosmological prediction. Проверить реальный новый GitHub Actions, затем построить высокопорядковый bounded/reversible timestepper либо оценить dt Lie-splitting со сменой event order и Richardson только после проверки асимптотики, сравнить full spectra на 3+ h/dt/Emax. Затем полнофизические ν/e±γ, dark χ1/χ2/φ и FRW. **CLASS-ready=false.**

## Поправки

Продолжает, а не переписывает [аудит spurious k=Emax invariants](2026-11-22-vector-neutrino-invariants.md); сохраняет два прежних отрицательных frozen modes и все численные ошибки finite-grid.
