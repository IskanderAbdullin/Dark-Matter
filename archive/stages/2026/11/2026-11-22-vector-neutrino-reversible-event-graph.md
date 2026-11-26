# Этап: 2026-11-22 — неравновесная microreversible Z′–ν–antiν схема

- ID: DM-KINETICS-20261122-REVERSIBLE-EVENT-GRAPH. Точная UTC длительность не фиксировалась
- Направление: kinetics; статус: независимый локальный численный пилот; новый GitHub Scientific CI проверить после коммита
- Начальный SHA: 64c39424b773c2f1e93127c1c0060fa44d58d0a7, ветка research/kinetic-equilibrium-audit-20261121
- Навигация: [GOAL](../../../../GOAL.md), [PLAN](../../../../PLAN.md), [специальный GOAL](../../../../kinetics/research-control/GOAL.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json)

## Предварительный научный вопрос и альтернативы

Можно ли эволюционировать разные произвольные fν(k), fantiν(k), bV(p), без скрытого FD thermal reset, с физическим on-shell Z′↔ν+antiν и точными дискретными number/energy/entropy constraints? Критерии: Pauli 0≤f≤1, Bose b≥0; точные C=Nν+Nantiν+2NV, L=Nν−Nantiν и энергия, монотонная энтропия; thermal fixed point μV=μν+μantiν; честные h/dt/cutoff и независимый непрерывный weak-контроль. Альтернативы: прямое непрерывное strong collocation, IMEX, conservative projection и higher-order implicit method. Отдельно проверить, что положительность и H-теорема НЕ гарантируют континуальной точности.

## Физическая схема и методы

Массивный spin-1 Z′ gV=3, два активных ν и два antiν состояния (gν=gantiν=2). Узкая модель Z′↔ν+antiν, нет других SM/χ каналов и расширения. Энергетическая lattice E_V,j=jh, kν,i=ih, kantiν,l=lh, j=i+l; при p_j=sqrt(E_j²-mV²) угловое условие abs(k_i-k_l)≤p_j. Положительные веса W_jil=wV,j Γ0(mV/E_j)h/p_j, wV=3pE h/(2π²), wnu=2k²h/(2π²). F=b(1-fν)(1-fantiν), R=(1+b)fνfantiν. Каждый event extent ξ удовлетворяет имплицитному ξ=dt W[F(ξ)-R(ξ)] и находится Brent-bracket в точной области допустимых occupations без posthoc clipping. Монотонность F↓,R↑ обеспечивает уникальный root; dS/dξ=ln(F/R), так что каждый элементарный шаг не уменьшает энтропию. Тождество E_j=k_i+k_l обеспечивает энергию без проекционных поправок.

Источники: [Froustey, Pitrou & Volpe (2020)](https://arxiv.org/abs/2008.01074), [Discrete Quantum Kinetic Equation (2023)](https://doi.org/10.1007/s44007-023-00070-1). Исходники статей не копировали, применяли собственную реализацию. Γ0=1 — только условная единица времени, не физическая Γelastic/H.

## Данные и исполнение

[Скрипт](../../../../kinetics/vector_neutrino_reversible_event_graph.py), [тесты](../../../../tests/test_vector_neutrino_reversible_event_graph.py), [отчёт](../../../../kinetics/VECTOR-NEUTRINO-REVERSIBLE-EVENT-GRAPH-AUDIT.md), [JSON](../../../../kinetics/results/vector-neutrino-reversible-graph-2026-11-22.json).

Команды: python -m pytest -q tests/test_vector_neutrino_reversible_event_graph.py; python -m kinetics.vector_neutrino_reversible_event_graph.

NumPy, SciPy, pytest; случайность отсутствует, seed не нужен. mV=.4, Emax=4, h=.4/.2/.1/.05, Γ0=1, T=1. Синтетические начальные fν=expit(−k+.15exp(−k)), fantiν=expit(−k−.1exp(−k/2)), bV=.75 BE(E). Непрерывный независимый reference 144×144 Gauss incoming p/CM angle. Time to Γ0t=4; dt=1,.5,.25,.125. Дополнительный cutoff-контроль Emax=4/6/8.

## Фактические результаты и численные ошибки

При h=.2, 189 микро-событий: начальные C=0.677688967485, L=0.00629871876221, E=1.087744100397, S=1.617320525670. В локальном Python-прототипе все C/L/E на конце сохраняются до abs O(1e−16), entropy монотонна, Bose/Pauli bounds соблюдены; chemical thermal fixed point и empty-Z′ inverse creation прошли. **Новый репозиторный код ещё требует реального GitHub CI — локальную проверку не выдавать за зелёный CI.**

Результат независимого continuous weak при Emax=4: dNV=.0058088005769424, dEV=.0115779873168571, dS=.00144991988902934. Относительные ошибки дискретного dNV на h=.4/.2/.1/.05: −21.0550%, −8.0845%, −0.8873%, +0.3778%. При h=.05 dEV ошибка +0.4095%, dS +0.3305%. Грубая h=.4 энергия ошибается −12.66% даже при точной дискретной энергии. Emax4 vs Emax8: источник N изменяется на ~14%, энергия ~29%, нет хвостовой сходимости.

При dt=1/.5/.25/.125 low-E bV после τ=4: .8887712453/.8722947720/.8615340320/.8578753698; существенный timestep bias. Внутренняя энтропия выросла +.000693635/+.000748727/+.000788669/+.000812388 соответственно. Это **не** физическая временная сходимость.

## Отрицательные результаты, ограничения

Exact lattice energy, positivity, quantum H theorem НЕ доказывают точность непрерывного collision: >12–21% ошибка на грубой сетке, ~.3–.4% даже на h=.05, nonmonotone crossing и нет строгого continuum bound. Сильный dt bias последовательного имплицитного расщепления; требуются адаптивные q/dt/cutoff и другие сетки. Полная SM ν–e±γ thermal chemistry, Z′ dark channels, χ1/χ2/φ Bose, Friedmann/entropy, Γelastic, 86-eV f(q), ΩDM/BEC/CLASS отсутствуют. CLASS-ready=false.

## Вывод и следующий исполнимый шаг

Численно построена локально воспроизводимая **консервативная дискретная схема**, но не validated continuum time evolution. Главные ворота PLAN A/B не меняются; добавлены отрицательные численные тесты. Сначала получить реальный Scientific CI; если failure, сохранить причину и исправить без ослабления критериев. Затем — независимая 3+ сеточная оценка и angular projection, higher-order coupled implicit/IMEX с bounded populations, energy/charge/entropy; проверить реальную динамику νbar/ν и восстановить weak источники. Только после этого физическая ν–SM и χ1/χ2/φ+FRW. Финальный SHA — коммит Git с этим файлом.

## Поправки

Не заменяет [прошлый nonthermal strong stage](2026-11-22-vector-neutrino-nonthermal-strong.md) и его заархивированный отрицательный CI/syntax fix; нет нового terminal dark matter вывода.
