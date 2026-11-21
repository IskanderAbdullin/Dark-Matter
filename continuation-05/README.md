# Теоретическое исследование V — положительные Bose cutoffs, устойчивость и открытая кинетика

Дата: 21 ноября 2026. Прямой исходный источник: Abdullin & Popov, Boson dark matter halos with a dominant noncondensed component, JCAP 11 (2021) 055, https://doi.org/10.1088/1475-7516/2021/11/055.

Это продолжение предыдущих исследований continuation-01–04. Цель: изучить зависимость реального размера и массы гало от регулярности энергетической границы, различить collisionless stationarity, quantum collision equilibrium и actual tidal escape.

## Что получено

1. [15 — Positive Bose Taylor remainder hierarchy](15-general-higher-order-lowering.md): строго положительные F_k(E) с df/dE<0, regularity C^(k−1), внешний показатель плотности n=k+3/2; k3 может дать ~10^9R, k4–6 не дали внешнего нуля до 10^10R в нескольких численных сценариях.
2. [16 — Log-radius Emden phase plane](16-phase-plane-and-outer-tail.md): автономная система в фазовой плоскости, след и определитель матрицы Якоби, связь со стандартным порогом n=5.
3. [17 — Existing classical stability criteria](17-classical-linear-stability.md): Doremus–Feix–Baumann/Antonov/Pérez–Aly уже дают достаточные условия линейной collisionless stability для отдельных bounded isotropic f(E); точный класс возмущений и regularity нужно сверить, singular critical center не покрыт автоматически.
4. [18 — Singular subcritical-to-critical limit](18-singular-critical-limit.md): fmax растёт как 1/delta, центральная quadratic curvature как 1/sqrt(delta), в то время как внешний радиус стремится к finite value.
5. [19 — Closed Bose detailed balance](19-closed-collision-equilibrium-no-go.md): жесткий lowered F_k(E) НЕ является стационарным распределением для полной закрытой эргодической 2→2 Uehling–Uhlenbeck collision dynamics; открытый tidal halo может быть квазистационарен.
6. [20 — Open energy-space flux](20-open-flux-energy-edge.md): finite diffusion с nonzero steady escape flux выбирает линейное падение f на границе; C1 quadratic cutoff совместим с этой игрушечной моделью только при нулевом граничном потоке. Получено exact Riccati solution, numerical match 3.2e-14 и finite interior pole.
7. [21 — Rein–Rendall compact-support theorem audit](21-mathematical-compactness-prior-art.md): существующие результаты по классическим isotropic polytropes при k=7/2 и границы переноса на Bose F_k.

## Код и проверки без GitHub Actions

- [reproduce_cutoff_ladder.py](reproduce_cutoff_ladder.py): positive cutoffs k1–6, численные решения Poisson, tail slopes.
- [critical_limit.py](critical_limit.py): bounded→singular critical transition.
- [reproduce_detailed_balance.py](reproduce_detailed_balance.py): непосредственная проверка 2→2 energy- and momentum-conserving scattering witness.
- [reproduce_open_flux.py](reproduce_open_flux.py): Riccati exact solution vs independently integrated ODE.
- [results/reproduced_profiles.json](results/reproduced_profiles.json): основные численные профили.
- [results/critical_limit.json](results/critical_limit.json): десять точек subcritical→critical с растущей fmax и сходящимся радиусом.
- [results/detailed_balance.json](results/detailed_balance.json), [results/open_flux.json](results/open_flux.json): collision/flux воспроизведение.
- [RESEARCH-LOG.md](RESEARCH-LOG.md): журнал, пределы доказательности, неразрешённые задачи.

Для запуска: Python 3.10+, установить зависимости из requirements.txt; выполнить приведённые Python-программы. Полные расширенные программы с 24+ радиальными примерами и численными наборами также помещены в независимый ZIP-архив, доступный в сообщении ChatGPT. GitHub содержит компактные воспроизводящие скрипты. Без CI, без запуска роботов, без вмешательства в main.

**Главный итог:** положительное и даже гладкое f(E) само по себе не означает компактное, термализованное или физически воспроизводимое гало. Напротив, искусственное повышение гладкости границы может создать огромный внешний хвост. Локальная открытая модель с конечным escape flux выбирает линейный, а не quadratic cutoff, но её результат нельзя экстраполировать на реальный галактический transport без полноценного quantum collision kernel и orbit averaging.

**Остаётся:** единая физическая модель рождения и поздней кинетики 86-eV бозона, GP+Boltzmann+Poisson с действительным escape, применение строго сформулированных stability theorems, cosmological transfer, статистическая проверка галактик. Ни одна новая заметка не считается рецензированной публикацией.
