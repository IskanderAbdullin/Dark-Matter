# Этап: 2026-11-22 — нетермальные χ1/χ2–ℓ физические Pauli 2→2 слабые моменты

- **ID / UTC:** DM-KINETICS-20261122-NONTHERMAL-CHI-LEPTON-WEAK; точный UTC интервал локального прогона не фиксировался
- **Направление, статус:** kinetics, independently derived local 5D pilot; GitHub Scientific CI нового коммита необходимо проверить, не считать успех до фактического run
- **Стартовый SHA/ветка:** [b820abf](https://github.com/IskanderAbdullin/Dark-Matter/commit/b820abf487b1ffd4c85ec38c13239cb9e3cdd1a9), research/kinetic-equilibrium-audit-20261121; итоговый SHA — commit этого этапа
- **Контекст:** [GOAL](../../../../GOAL.md) → [PLAN](../../../../PLAN.md) → [PROJECT_MAP](../../../../PROJECT_MAP.md) → [специальное ТЗ](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md) → [научные запреты](../../../../kinetics/RESEARCH_GUARDRAILS.md) → [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json)

## Вопрос, гипотеза и альтернативы

Важнейший пропуск в текущей физике: численные finite-mV 5D χ+ℓ Pauli интегралы до сих пор принимали заранее назначенные FD χ/ℓ распределения и не работали с реальной неравновесной χ2 от Yukawa decay. Сформулирован ограниченный начальный gate: построить **тот же физический on-shell 2→2 weak source** с произвольными внутренними `fχ(E),fℓ(ω)`, finite mℓ и mV, с полной gain/loss, обратной реакцией и entropy; независимо сравнить FD-subset с прежним опубликованным в репозитории [full_heat](../../../../kinetics/finite_mediator_5d_adaptive_reference.py), сохранить ограничения по сильному `Cχ`. Альтернативы: явный paired event graph или Fokker–Planck, другая gauge/UV completion, физический 2→2 collision через BEST. Наличие weak heat не принимается за физическую gamma.

## Микроскопическая модель, математический метод и источники

Dirac χ1 и χ2 с одинаковым ненулевым vectorlike gauge charge (необходимость для *одной real neutral φ* показана [ранее](../../../../kinetics/REAL-PHI-YUKAWA-GAUGE-CHARGE-AUDIT.md)); лептон ℓ либо massless ν, либо массивный charged lepton, полный finite-mass t-channel vector propagator. Матрица:
`|M|²_avg=2G²[(s-M²-m²)²+(u-M²-m²)²+2(M²+m²)t]×(mV²/(mV²-t))²`.
`F=fχ(E)fℓ(ω)(1-fχ(E'))(1-fℓ(ω'))`, `R=fχ(E')fℓ(ω')(1-fχ(E))(1-fℓ(ω))`.
Окончательные энергии получаются из Lorentz boost CM scattering в произвольном азимуте; `E+ω=E'+ω'`. Weak `Qχ∝∫(E'-E)(F-R)|M|²`, `Qℓ=-Qχ`, entropy `∝∫(F-R)log(F/R)|M|²≥0`. При общем T и *независимых* μχ и μℓ 상세ный balance точен, даже без µχ=µℓ. Полный `Cχ(p),Cℓ(k)` пока не вычислен; spin degeneracies прежние 4×4 — только **согласованный численный reference convention**, не универсальный ν flavor counting.

Ссылки: [базовый 5D angular reference](../../../../kinetics/finite_mediator_forward_event_reference.py), [микрофизический Yukawa portal](../../../../kinetics/MICROPHYSICS-MODEL.md), [допустимость charges](../../../../kinetics/REAL-PHI-YUKAWA-GAUGE-CHARGE-AUDIT.md). Код старых моделей не удалялся и исторические отрицательные случаи не перезаписаны.

## Первичные локальные результаты, ошибки и воспроизведение

- [Код с arbitrary isotropic Fermi occupation callbacks](../../../../kinetics/nonthermal_chi_lepton_5d_weak.py)
- [Сравнительные тесты](../../../../tests/test_nonthermal_chi_lepton_5d_weak.py)
- [Научный отчёт](../../../../kinetics/NONTHERMAL-CHI-LEPTON-5D-WEAK-AUDIT.md)
- [Машинный JSON](../../../../kinetics/results/nonthermal-chi-lepton-5d-weak-2026-11-22.json)
- `python -m pytest -q tests/test_nonthermal_chi_lepton_5d_weak.py`
- `python -m kinetics.nonthermal_chi_lepton_5d_weak`

Локальный независимый bootstrap использовал NumPy/SciPy, deterministic `nr=4,nμ=4,nφ=12,ncos=12`, `cut_kinetic_T=12`, `Tlepton=1,mV=.1,G=.001`, `Mχ1=.1,Mχ2=.02,mℓ=0,.5` в общей условной системе энергетических единиц. RNG/seed нет; лицензируемые PDFs не копировались. Для FD `fχ=FD(T=.9)`, `fℓ=FD(T=1)` при χ1 получено Qχ=+6.935904667×10⁻¹⁴, S=+7.706560741×10⁻¹⁵ (mℓ=0); Qχ=+6.088215658×10⁻¹⁴, S=+6.764684064×10⁻¹⁵ (mℓ=.5). При общем T=1 Qχ≈−2.8×10⁻²⁹ или −3.5×10⁻²⁹, entropy≈10⁻⁴². При синтетических нетермальных `fχ=expit(-E/.9+.15exp(-E/.9))` и `fℓ=expit(-ω+.05exp(-ω))` Qχ1 для mℓ=0/+mℓ=.5 =7.725499230×10⁻¹⁴/6.839775255×10⁻¹⁴; Qχ2=7.799808434×10⁻¹⁴/6.889894767×10⁻¹⁴. Все entropy положительны, `Qχ+Qℓ=0` аналитически.

Новый сравнительный CI должен воспроизвести **специализацию** произвольного fχ/fℓ к FD и проверить вычисления с [существующей программой](../../../../kinetics/finite_mediator_5d_adaptive_reference.py) при совершенно одинаковых сетках и cutoff, без ослабления точности. Данный dataset **не** независимая физическая проверка χ2 abundance и не physical gamma. Числа nrad4 **не доказаны radial/UV/azimuth сходящимися**.

## Отрицательные данные и неизменённая цель

Строго отсечены Fermi-граничные 0,1 — без регуляризованного entropy `xlogy` такая точка `fail closed`, не скрытая замена. Нет mode-resolved `Cχ1(p),Cχ2(p),Cℓ(k)`, coupled finite ν/e±γ background или actual fχi после χ1↔χ2+φ; нет ν flavor/gauge thermal masses, mediator on-shell population and cosmology, нет physical Гamma_el, freeze-out и абсолютной Ω. Gauge-consistent φ real Yukawa **требует** обе χ scattering; формулы для двух видов есть, но самосогласованного source пока нет. Нельзя выдавать arbitrary fχ2 при одинаковой нормировке за его настоящую abundance или доказывать массовую 86eV Ω/BEC/CLASS.

## Вывод, следующий исполнимый этап и поправки

Получен сквозной *неравновесный weak 2→2 phase-space functional* с конечными массами и полной Pauli, необходимый для будущего полноценного оператора. Следующее P0 — вывести **парные strong `Cχ(p)` и `Cℓ(k)`** из этого же event integral, добившись number/energy/entropy weak-strong equality с независимыми физическими мерами и хотя бы 3 разрешениями; затем одновременно χ1+χ2, ν/e±γ, Z′ и φ Bose FRW. Только после этого обсуждать физическую Γelastic/H и late 86eV f(q). Обновить [STATE](../../../../kinetics/research-control/STATE.json), [EVIDENCE](../../../../kinetics/research-control/EVIDENCE.md), [BACKLOG](../../../../kinetics/research-control/BACKLOG.md), подробный [план](../../../../plans/RESEARCH_ROADMAP_DETAIL.md) только если меняются научные ворота. `CLASS-ready=false`.

**Связь:** это НЕ заменяет ранее выполненные 5D thermal benchmarks, [χ2 conditional thermal heat](../../../../kinetics/CHI2-LEPTON-HEAT-GATE.md) и отрицательный UV selection audit. Внешние временные CI artifacts не являются достаточным постоянным архивом.

## Реальный CI после публикации (append-only)

- На [b4d6549](https://github.com/IskanderAbdullin/Dark-Matter/commit/b4d6549504f60f8dbf87658db2f3df0e7aaa4ca7) успешно завершён [Scientific verification #45016136831](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45016136831). Это scoped weak 5D Pauli functional, а не физический strong Cχ/Cℓ или 86eV dark matter.
