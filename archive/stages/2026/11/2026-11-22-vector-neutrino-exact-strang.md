# Этап: 2026-11-22 — точный локальный Riccati и симметричный квантовый Strang

- **ID / UTC:** DM-KINETICS-20261122-EXACT-MICROEVENT-STRANG. Дата 2026-11-22, точные UTC отметки начала и окончания эксперимента не регистрировались.
- **Направление и статус:** kinetics; локально независимо вычислен и проверен scoped finite-event experiment, новый GitHub Scientific CI до повышения статуса проверить.
- **Рабочая ветка, стартовый HEAD:** `research/kinetic-equilibrium-audit-20261121`, [92caea7ad2ae](https://github.com/IskanderAbdullin/Dark-Matter/commit/92caea7ad2ae4fa6f04d384e9fc37ea34b7ec2c7). Итоговый SHA — Git commit, содержащий этот stage.
- **Общая цель:** [GOAL](../../../../GOAL.md) → [PLAN](../../../../PLAN.md) → [подробный маршрут](../../../../plans/RESEARCH_ROADMAP_DETAIL.md) → [kinetics GOAL](../../../../kinetics/research-control/GOAL.md) → [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Исследовательский вопрос, гипотеза и альтернативы

После [независимого равновесного контроля](2026-11-22-vector-neutrino-constrained-equilibrium.md) установлено, что каждый элементарный backward-Euler шаг микроскопически безопасен, но сохраняет крупное и иногда немонотонное `dt` расщепление. До испытания сформулирован критерий: заменить **каждый микропоток** на точный аналитический поток конечномерного quantum Pauli/Bose 1↔2, сделать палиндромную Strang композицию и проверить методом независимого `solve_ivp`, 3+ разрешениями по времени и инвариантами. Альтернативы: общий implicit coupled solver, IMEX, изменение grid/cutoff, отложенная физическая модель. Позитивность, conservation и H-теорема **не** принимаются за доказательство order/continuum accuracy.

## Физическая модель, схема и первичные источники

Изолированная и условно нормированная broken-phase 1↔2 `Z′↔ν+antiν`, gV=3, gν=gantiν=2, конечная lattice `E_j=kν_i+kantiν_l` с положительными exact daughter-angle overlap и clipped radial density weights. Γ0=1 означает **произвольную единицу времени**, не физический `Γ/H`. Для каждого события с интегральными весами `wV,wν,wantiν` и extent `x`: `b=B−x/wV`, `f=F+x/wν`, `a=A+x/wantiν`. Точный quantum flux `J=b(1-f)(1-a)−(1+b)fa` — квадратный полином по x. Уравнение `dx/dt=W J(x)` интегрируется аналитической Riccati формулой; тогда forward `dt/2` по всем событиям и backward `dt/2` дают Strang композицию **точных микрореакций**. Скалярный случай independently проверен `solve_ivp`. Каждый элементарный поток сохраняет C=Nν+Nantiν+2NV, L=Nν−Nantiν, E; entropy не убывает внутри физической области. No silent clipping.

Общая методическая литература: [Strang splitting в кинетике](https://www.sciencedirect.com/science/article/pii/S0898122111007206), [детерминированные консервативные Boltzmann методы](https://www.sciencedirect.com/science/article/pii/S0021999103000652). Не копировали сторонний исходный код, формулы и solver выведены для собственной reaction network.

## Входы, воспроизведение и архивированные файлы

- [Исполняемая реализация](../../../../kinetics/vector_neutrino_exact_event_strang.py).
- [Регрессионные тесты](../../../../tests/test_vector_neutrino_exact_event_strang.py).
- [Научный отчёт](../../../../kinetics/VECTOR-NEUTRINO-EXACT-EVENT-STRANG-AUDIT.md).
- [Машинные численные результаты](../../../../kinetics/results/vector-neutrino-exact-strang-2026-11-22.json).

Команды для нового агента: `python -m pytest -q tests/test_vector_neutrino_exact_event_strang.py`; `python -m kinetics.vector_neutrino_exact_event_strang --h 0.2 --tau 4 --dt 0.0625`. Независимый путь: `scipy.integrate.solve_ivp` по физическому `J(x)` с `rtol=2e-11`, `atol=2e-17`; 6 event × 3 dt. Runtime Python NumPy/SciPy/pytest, без RNG/seed. Синтетический `T=1`, `mV=.4,Emax=4,h=.4,.2,.1`, начальные `fν=expit(−k+.15e^−k)`, `fantiν=expit(−k−.10e^−k/2)`, `b=.75BE(E)`. `τ=4`, `dtref=1/256`. Это НЕ лабораторные, космологические или наблюдательные данные.

## Фактические результаты и численные ошибки

Абсолютная разница exact event extent vs независимого `solve_ivp` для 18 событийных экспериментов ≤3.20×10⁻¹⁶. Для h=.2, 207 событий, максимум абсолютного спектрального отличия всех трёх распределений от Strang dtref=1/256: `dt=.0625→.00727514`, `.03125→.00209798`, `.015625→.00052151`, `.0078125→.00010547`. Отношения первых трёх ~3.47 и 4.02, что согласуется с асимптотическим порядком 2 **на данном конечном графе и этом диапазоне dt**. Forward/reverse full event при dt=.0625: `.01257733`/`.01487057` ошибки, против Strang `.00727514`.

В диапазоне h=.4 и `dt=.125,.0625,.03125` спектральные ошибки `.00174420,.00045490,.00011370`. C,L,E сохраняются до ~10⁻¹⁵, entropy растёт на каждом шаге, fermion Pauli и vector Bose границы не нарушаются. Отдельные тесты химической fixed point и empty-vector inverse production сохранены.

## Отрицательные результаты и область применимости

При **h=.1**, 808 событий, `dt=.125,.0625,.03125` ошибки vs `dtref=1/256` остаются `.0232120,.0141676,.0072225`: наблюдаемые отношения **1.64 и 1.96**, а не 4. Stiffness/неасимптотический dt режим и finite cutoff не решены. Ранее установленный angular/threshold/number/cutoff error сохраняется независимо от time-order. На успех этого discrete 1↔2 solver нельзя ссылаться как на доказательство полной релятивистской χ1/χ2-SM Pauli 2→2, Γ_el(T), физической 86-eV φ Bose, FRW, ΩDM, BEC, frozen f(q) и CLASS-ready. Текущий статус `in_progress`, `CLASS-ready=false`.

## Следующий исполнимый шаг и checkpoint

Проверить фактический CI нового HEAD; при failure сохранить traceback и исправить исходник/тест **не ослабляя физические критерии**. Затем измерить жёсткий h=.1 Strang на dt≤1/256 против существенно более мелкого эталона, добавить 3+ h/Emax/cutoff и shape проверки со spectrum-level нормами, рассмотреть adaptive substepping. Лишь затем возвращаться к физическим слабым ν–e±γ, χ1/χ2 Pauli 2→2 и real φ Bose/FRW. Результаты и следующий шаг обновляются в [STATE](../../../../kinetics/research-control/STATE.json), [EVIDENCE](../../../../kinetics/research-control/EVIDENCE.md), [BACKLOG](../../../../kinetics/research-control/BACKLOG.md). Итоговый SHA доступен в Git истории.

## Поправки и связь с прежними этапами

Этап продолжает [прежний неравновесный граф](2026-11-22-vector-neutrino-reversible-event-graph.md), [точную angular geometry](2026-11-22-vector-neutrino-shell-overlap.md), [исправление Emax](2026-11-22-vector-neutrino-radial-endpoints.md) и [нулевую моду stoichiometry](2026-11-22-vector-neutrino-invariants.md). Не исправляет прошлые h/Emax отрицательные данные задним числом. Внешние Actions artifacts не являются единственным архивом.
