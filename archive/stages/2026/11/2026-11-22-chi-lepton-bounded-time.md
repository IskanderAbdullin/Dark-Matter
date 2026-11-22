# Этап 2026-11-22: нелинейное Pauli-bounded time-evolution χ–ℓ 2→2 из полной 5D микрофизики

- **ID / UTC:** DM-KINETICS-20261122-CHI-LEPTON-BOUNDED-TIME; точное UTC время отдельных прогонов не записано, длительность не выдумывается
- **Направление / статус:** kinetics; ограниченная локальная независимая математическая и численная проверка, CI нового commit ещё требуется подтвердить
- **Ветка / стартовый SHA:** research/kinetic-equilibrium-audit-20261121, [8bc98b1](https://github.com/IskanderAbdullin/Dark-Matter/commit/8bc98b113bbef8e55d8f17640b771df08a744010); итоговый commit содержит эту запись
- **Цель → план:** [GOAL](../../../../GOAL.md), [PROJECT_MAP](../../../../PROJECT_MAP.md), [PLAN](../../../../PLAN.md), [подробная научная программа](../../../../plans/RESEARCH_ROADMAP_DETAIL.md), [специальное кинетическое ТЗ](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md), [научные инварианты](../../../../kinetics/RESEARCH_GUARDRAILS.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json)

## Предварительная гипотеза и измеримые ворота

Энтропийно-согласованная χℓ 2→2 проекция [прошлого этапа](2026-11-22-chi-lepton-logit-entropy-repair.md) доказала H на мгновенном источнике, но не гарантировала ни физических Pauli bounds после dt, ни сохранения N/E при реальном time evolution. Цель: пересчитывать полный on-shell finite-vector 5D Pauli source от **текущего** fχ,fℓ, использовать истинный численный density-of-states вес и принимать только корректные по 0<f<1, полной энтропии и спектральной локальной ошибке временные шаги. Альтернативы: implicit event graph, entropy-stable IMEX, exact reaction quadrature. Критерий: настоящая адаптивная rejection, не имитация; численная Cχ,Cℓ/E и entropy до проверяемых ошибок, отдельная отрицательная граница h/UV.

## Математическая физика, масса и фазы

Каждый E-bin веса w_i= g/(6π²)[(max(0,U²−m²))^(3/2)−(max(0,L²−m²))^(3/2)] с g=4, L=max(m,E_i−h/2), U=min(Emax,E_i+h/2). Массы χ=.1 и лептона=0 в условных единицах; gχgℓ=4×4 пока только прежняя numeric нормировка, не correct ν flavor. Текущий непрерывный collisional F/R получает f_event=expit(Σα logit(f_i)) на точных relativistic on-shell outgoing, векторный propagator finite mV; такое же α для deposition χ и ℓ dN. Динамика df_i=dt dN_i/w_i, текущие distributions меняются. Доказательство мгновенной H в [логит-этапе](2026-11-22-chi-lepton-logit-entropy-repair.md); finite step проверяется реальной entropy из xlogy. Отказ от невалидного шага, снижение dt, без clipping, ресета FD или posthoc conservation repair. Два half-Euler vs один Euler на all species max-spectrum error, прошедший local tolerance даёт два half-Euler; итог остаётся **first-order method** и не гарантирует global convergence.

## Воспроизводимые источники, входы и команды

[Исходный код](../../../../kinetics/paired_chi_lepton_bounded_time.py), [проверки](../../../../tests/test_paired_chi_lepton_bounded_time.py), [подробный аудит](../../../../kinetics/PAIRED-CHI-LEPTON-BOUNDED-TIME-AUDIT.md), [машинный JSON](../../../../kinetics/results/paired-chi-lepton-bounded-time-2026-11-22.json).

Команды: python -m pytest -q tests/test_paired_chi_lepton_bounded_time.py; python -m kinetics.paired_chi_lepton_bounded_time. Числовые параметры: Mχ=.1,mℓ=0,G=.001,mV=.1,Tref=1,nrad4,nmu4,ncos12,nφ12,cut12,N=10 energy intervals, numerical time1e10, initial step1e10, atol=1e−6,1e−7,1e−8. NumPy/SciPy, детерминированная схема без seed/MC. Временной масштаб никакого отношения к physical cosmological seconds/Γel/H не имеет. Локальная независимая sandbox reconstruction выполнена; данные и официальный скрипт в Git, внешний scratch не архивирован (not_archived: временный дубликат, не является существенным первичным набором).

## Факты и ошибки

Oscillatory nonthermal fχ=expit(-E+.3sin4E), fℓ=expit(-ω+.3cos4ω): atol1e−6/1e−7/1e−8 → accepted4/11/33, rejected2/7/15, maximum local spectral error 9.70e−7/9.77e−8/9.79e−9, entropy gain +4.21465e−4/+4.21146e−4/+4.21049e−4. Smooth nonthermal fχ=expit(-E/.9+.15exp(-E/.9)), fℓ=expit(-ω+.05exp(-ω)): accepted2/7/24, rejected1/4/12, S gain +9.94502e−5/+9.92979e−5/+9.92540e−5. Для всех проверенных прогонов сохранены weighted Nχ,Nℓ и energy с abs≤4.5e−16, entropy монотонна, все f остаются строго внутри (0,1). Значения, units, границы и следующие exact commands — в JSON/аудите. Нейтринный и χ2 спектры реальной Вселенной не заданы.

## Отрицательный и физический вывод

Это действительно **двухвидовая временная Pauli-схема для условного χ+ℓ**, а не terminal χ1/χ2 Yukawa+SM FRW. Важнейшая численная ошибка уже доказана предыдущим stage: при 10 интервалах осциллирующее nonthermal f даёт **−93.4%** ошибку энтропийного production относительно непрерывного physical event integral; 320 интервалов ещё −1.38%, и входной 5D nr4 не сходится. Новый safety/time-step controller не исправляет ни одну из этих сеточных ошибок. Реальные SM flavors, UV-valid g,mV, produced fχ2, on-shell Z′, χ1↔χ2+86eV φ Bose/Pauli, ν/e±γ reservoir, H, Γelastic, ΩDM, BEC и CLASS отсутствуют. CLASS-ready=false. Не называть accepted atol глобальной time error guarantee.

## Следующий исполнимый шаг

Проверить настоящий GitHub Scientific CI, сохранить traceback, если упадёт. P0: независимая 3+ energy bins/radial/angle/Emax and at least two smooth/oscillatory Fermi shapes, full weak/strong N/E/S and global time spectrum error, first order vs suitable higher-order entropy-safe method; затем физическое χ1+χ2+real φ Bose и full SM/Z′ FRW. Не переписывать отрицательные выводы; [EVIDENCE](../../../../kinetics/research-control/EVIDENCE.md) и [BACKLOG](../../../../kinetics/research-control/BACKLOG.md) являются сохранённым источником истории.

## Поправки

Продолжает, не отменяет [naive interpolation H-theorem counterexample](2026-11-22-chi-lepton-conservative-entropy-negative.md) и [logit interpolation partial repair](2026-11-22-chi-lepton-logit-entropy-repair.md). Конкретная новая проверка доказывает только Pauli-bounds, discrete conservation and entropy в приведённых dt/grid режимах.
