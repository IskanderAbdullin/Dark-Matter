# Этап: 2026-11-22 — искусственные нейтринные conserved modes в Z′ collision lattice

- **ID / UTC:** DM-KINETICS-20261122-NU-GRAPH-INVARIANTS; точное UTC время расчёта не измерялось
- **Область/статус:** kinetics, независимый математический анализ + локально исполненные численные и exact-rank контроли; новый GH CI требуется проверить
- **Исходный HEAD:** research/kinetic-equilibrium-audit-20261121, [b37b7dd](https://github.com/IskanderAbdullin/Dark-Matter/commit/b37b7dd1715340680962db6061e78ffda1e30ff2); итоговая ссылка — обычный коммит этого этапа
- **GOAL/PLAN/NEXT_ACTION:** [общая цель](../../../../GOAL.md), [живой план](../../../../PLAN.md), [кинетическая цель](../../../../kinetics/research-control/GOAL.md), [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json)

## Предварительная гипотеза, независимый эталон и альтернативы

Физический constrained equilibrium не должен быть признан достижимым только потому, что reaction solver сохраняет C=Nnu+Nantinu+2 NV, L=Nnu−Nantinu, E и растит энтропию: конечная lattice может иметь дополнительные запрещённые физически conserved modes. Проверяемая гипотеза: исходная целочисленная k-grid с kmax=Emax замораживает последние отдельные ν и antiν modes, а смещённая half-offset-grid устраняет именно эти два null направления. Критерий: разложить точную целочисленную incidence матрицу S и определить полную dim ker(S^T) с независимыми методами, сравнить с аналитическими C,L,E, найти изолированные строки. Альтернативы: численный SVD threshold artifact, missing angular support, запрет decay at different mV/Emax, реально сохранённая физическая charge.

## Модель/доказательство

Z′↔ν_i+antiν_l на on-shell E_j=k_i+k_l, столбец S=(-e_V,+e_ν,+e_anti). C=(2,1,1), L=(0,1,-1), E=(E_j,k_i,k_l) являются точными left-null vectors. Любой дополнительный y с S^T y=0 сохраняется независимо от kinematical distribution, dt и энтропии. Для integer ν-grid k=ih (i≥1), последний ν=k=Emax не может быть дочерним при E_V≤Emax и другой дочери с k>0, а также последний antiν; обе строки в S тождественно нулевые. Для half-shifted k=(i+1/2)h верхняя мода k=Emax−h/2 может парно образоваться при E_V=Emax с k_other=h/2, если это разрешает |k1−k2|≤p_V. В заявленном mV=.4,Emax4 режиме разрешает. Равновесные/позитивные свойства ранее изученного discrete quantum event solver не меняют эти линейные факты.

## Источники, код и воспроизводимость

[Исполняемый стоихиометрический анализ](../../../../kinetics/vector_neutrino_event_invariants.py), [тесты](../../../../tests/test_vector_neutrino_event_invariants.py), [подробный отчёт](../../../../kinetics/VECTOR-NEUTRINO-EVENT-INVARIANTS-AUDIT.md), [JSON](../../../../kinetics/results/vector-neutrino-event-invariants-2026-11-22.json). Команды: python -m pytest -q tests/test_vector_neutrino_event_invariants.py; python -m kinetics.vector_neutrino_event_invariants.

Только NumPy/SciPy/test, seed не нужен. Два точных modular prime ранга 1000003 и 1000033 + floating SVD, дополнительный локальный independent SymPy exact rank at h=.4/.2 (SymPy не обязательная runtime зависимость). Входы mV=.4,Emax=4,h=.4,.2,.1, обе сетки с identical числами дискретных состояний.

## Фактические результаты и отрицательные проверки

h=.4: 29 states, integer edges45/rank24/nullity5, shifted edges54/rank26/nullity3. h=.2: 58 states, integer 189/rank53/nullity5, shifted 207/rank55/nullity3. h=.1: 116 states, integer 774/rank111/nullity5, shifted 808/rank113/nullity3. Оба prime ранга и SVD точно согласованы, локальный SymPy at h=.4/.2 подтверждает integer 24/53, shifted26/55. C/L/E имеют остатки S^T v≈0. У integer-grid две верхние отдельные ν и antiν строки S полностью нулевые: это **два дополнительно искусственных conserved occupation**. У shifted-grid в этих тестах нет isolated ν rows и nullspace dim=3.

**Негативная физическая трактовка:** первоначальная integer-grid может монотонно увеличивать H entropy, строго сохранять настоящие C/L/E и ВСЁ ЖЕ никогда не забыть произвольные значения двух конечных верхних occupation — они участвуют в ложных дополнительных химических ограничениях. Это не настоящие новые симметрии нейтринной физики. Никакое уменьшение dt или «более точный» ODE solver не исправит топологию reaction network.

## Вывод и последующая задача

Для этих параметров shifted half-offset angular-cell lattice является топологически более физически пригодной: она устраняет две доказанные frozen boundary invariants. Но nullity=3 для трёх тестовых h **не доказывает** continuum ergodicity, отсутствие малых собственных значений, динамической ν thermalization в Standard Model и физическую модель eV bosonic DM. Схема всё ещё требует full cutoff/dt/energy/spectral test и реального gauge/SM/χ1/χ2/φ+FRW. Следующий исполнимый этап: подтвердить новый CI, затем решить independent constrained entropy maximum FD/BE с C/L/E на shifted graph, интегрировать к нему реальное nonthermal time distribution при dt/q refinements и проверить отсутствие дополнительных almost-invariants. CLASS-ready=false.

## Поправки/связи

Исходный [reversible event graph](2026-11-22-vector-neutrino-reversible-event-graph.md) математически сохраняет C/L/E, но теперь документированы **две дополнительные искусственные invariant modes** в его integer сетке; не переписывать архив раннего этапа без append-only связи. Позитивный angular overlap [stage](2026-11-22-vector-neutrino-shell-overlap.md) получает новое независимое топологическое обоснование.
