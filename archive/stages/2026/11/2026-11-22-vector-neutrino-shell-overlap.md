# Этап: 2026-11-22 — точная геометрия угловых ячеек ν от Z′

- **ID / UTC:** DM-KINETICS-20261122-ZPRIME-NU-CELL-OVERLAP; дата 2026-11-22, точная продолжительность не измерялась
- **Направление, статус:** kinetics / независимая аналитика + локальная численная проверка; реальный CI нового исходника обязателен до повышения статуса
- **Ветка, исходный commit:** research/kinetic-equilibrium-audit-20261121, [974d6474](https://github.com/IskanderAbdullin/Dark-Matter/commit/974d6474cc6cb1427b5b8ff1bfccfe101de25954); итоговый SHA зафиксирован коммитом с этим этапом
- **Цель:** [GOAL](../../../../GOAL.md) → [PLAN](../../../../PLAN.md) → [kinetics GOAL](../../../../kinetics/research-control/GOAL.md) → [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json)

## Проверяемый вопрос, критерии и альтернативы

Предыдущий conservative event graph сохраняет дискретную энергию, Pauli/Bose и энтропию, но на грубой энергии ошибается вплоть до ~21% в числе. Гипотеза: физическая угловая ошибка частично связана с присвоением полной ширины ячейки h по одному попаданию центра в angular support. Независимый контроль: на *такой же* полуцелочисленной сетке сравнить center-hit и аналитическую длину пересечения ячейки с точной кинематической областью. Критерии до проверки: ∑_i Delta k_i=p_j на каждом векторном shell до округления, коэффициенты W>0, сохраняются eventwise C/L/E и H, отдельные относительные N/E/S к независимому p,cosθ weak benchmark 144×144; отрицательные остатки сохранить. Альтернативы: ν grid иной ориентации, conservative energy remap, midpoint+adaptive angular, coupled higher-order implicit.

## Физика и метод

Вектор E_j=jh, p_j=sqrt((jh)^2−mV²); полуцелые daughter energies kν,i=(i+1/2)h и kanti,l=(l+1/2)h, j=i+l+1. Никакого энергии remap: E=kν+kanti точно. Реальная допустимая область дочерней энергии K=[(E−p)/2,(E+p)/2], ячейка Ci=[ih,(i+1)h]. Точная длина d=length(K∩Ci)≥0, W_i = wV Γ0 (mV/E_j)d/p_j. Так как Ci разбивают [0,E_j], ∑d=p_j. Нет posthoc нормировки или отрицательных весов. Условия Pauli/Bose и implicit extent взяты из [первого microreversible пилота](2026-11-22-vector-neutrino-reversible-event-graph.md). Оговорка: это узкий Z′↔ν+antiν, нет e±γ, χ1/χ2/φ и расширения; Γ0=1 условная единица. Источник физического quantum collision: [Froustey/Pitrou/Volpe 2020](https://arxiv.org/abs/2008.01074). Чужой код не переносили.

## Постоянные исходные программы, данные и команды

[Скрипт](../../../../kinetics/vector_neutrino_shell_cell_overlap.py), [четыре теста](../../../../tests/test_vector_neutrino_shell_cell_overlap.py), [научный отчёт](../../../../kinetics/VECTOR-NEUTRINO-SHELL-CELL-OVERLAP-AUDIT.md), [JSON](../../../../kinetics/results/vector-neutrino-shell-overlap-2026-11-22.json). Запуск: python -m pytest -q tests/test_vector_neutrino_shell_cell_overlap.py; python -m kinetics.vector_neutrino_shell_cell_overlap.

Случайность отсутствует. mV=.4, Emax=4, h=.4,.2,.1,.05, T=1, Γ0=1, gV=3, gnu=ganti=2, synthetic fν/fantiν и .75BE как в предыдущем этапе, независимый continuum 144×144.

## Реальные результаты и независимая проверка

На каждой vector shell exact cell-overlap normalization ∑W/[wV Γ0(mV/E)]−1 имеет max abs <4e−15 для всех h. На той же половинной ν-сетке center-hit демонстрирует max ошибку 15.47%,34.16%,13.40%,21.27% соответственно h=.4,.2,.1,.05. Это аналитически прослеживаемая проблема геометрической квадратуры, не уравнения Bose.

При h=.2 center-hit relative dNV +5.2037%, overlap +0.0824%; relative dEV +5.0406%→+2.2292%; entropy +5.3046%→+0.8059%. При h=.05 number +0.4030%→+0.1640%, energy +1.0860%→+0.5922%, entropy +0.6363%→+0.2816%. Непрерывный reference dNV=.00580880057694,dEV=.01157798731686,dS=.00144991988903. Позитивный ограниченный вывод: точное нормирование угла устраняет существенный alias на этой сетке.

## Отрицательные результаты и ограничения

**Несмотря на ∑W exact, полная слабая энергия НЕ сходится по доказанному критерию:** h=.2 ошибка E +2.23%, h=.05 остаётся +0.592%. Векторный E/p midpoint, cutoff, низкий k, тепловые масcы и отдельные спектральные формы не доказаны continuum accurate. Нельзя сравнивать первую integer-ν lattice и эту shifted lattice как чистый эксперимент по angular weights: правильное сравнение — center-hit и overlap в ЭТОМ этапе при идентичной сетке. Старый dt-bias не исправлен. Настоящая темная материя χ1/χ2/φ, H, 86eV Ω/BEC/CLASS-ready отсутствует.

## Вывод, пересмотр плана и следующий этап

Получена независимая проверяемая поправка **к физически нормированному угловому shell** с сохранением позитивных весов и микроскопических инвариантов, но остаточная энергия и entropy bias требует улучшения энергетической меры, отсечения и интегратора времени. [PLAN](../../../../PLAN.md) дополнен коротким решением о приоритете continuum error перед космологическими интерпретациями. При новом запуске сверить реальные CI обоих новых этапов; затем расширить независимые h, angular, radial, cutoff и dt тесты, сравнить в одинаковой физической области с непрерывным weak solver; только после этого e±γν, χ1/χ2/φ и FRW. CLASS-ready=false.

## Поправки

Не заменяет [архив первого графа](2026-11-22-vector-neutrino-reversible-event-graph.md). Все отрицательные ошибки сохранены; отчёты/исходные коды/JSON версионированы в Git, не только в Actions.

## CI-проверка (append-only)

- На commit 174e9b4ed1fc2cc8b11f02caab7635cfda8ad65d фактически завершён успешный [Scientific verification #45011770137](https://github.com/IskanderAbdullin/Dark-Matter/actions/runs/45011770137). Подтверждена регрессия угловой клеточной схемы, но не утверждается континуальная точность физической φ/χ DM модели.
