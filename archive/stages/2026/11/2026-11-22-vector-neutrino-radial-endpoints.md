# Этап: 2026-11-22 — ошибка конечного Emax и положительная Z′-радиальная квадратура

- **ID / UTC:** DM-KINETICS-20261122-ZPRIME-ENDPOINT-MEASURE; точное UTC начало/конец не записывались
- **Область, статус:** kinetics; выполнен независимый локальный численный контроль, новый репозиторный CI ещё проверяется
- **Ветка / стартовый SHA:** research/kinetic-equilibrium-audit-20261121, [174e9b4e](https://github.com/IskanderAbdullin/Dark-Matter/commit/174e9b4ed1fc2cc8b11f02caab7635cfda8ad65d); итоговый SHA — fast-forward коммит этого stage
- **Цель/план:** [GOAL](../../../../GOAL.md) → [PLAN](../../../../PLAN.md) → [kinetics GOAL](../../../../kinetics/research-control/GOAL.md) → [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json)

## Гипотеза и независимый критерий

После точной угловой нормы h=.05 ошибка E-момента остаётся +0.592%. До расчёта: проверить, не является ли главным источником систематического сдвига полная ширина h, присвоенная последнему векторному узлу E_j=Emax при сравнении с эталонным непрерывным интегралом только по E≤Emax. Альтернативы: оставшаяся погрешность k интерполяции, p/E threshold, угловая неравномерность, finite cutoff, спектральный dt-bias. Критерий: на той же ν/antiν energy lattice и том же F−R отдельно сравнить оригинальный полный endpoint, верхнюю половинную трапецию и физически точный интеграл vector 0th density over clipped E-bin. Проверить три и более h и оставить даже случаи ухудшения N/entropy.

## Физика, нормировка и независимые ссылки

Вектор gV=3, два активных ν и antiν (каждой g=2); только Z′↔ν+antiν, Γ0=1 условная единица. Физическая энергия E=sqrt(p²+mV²) даёт d³p/(2π)³=E sqrt(E²−mV²)dE/(2π²). Положительная exact cell 0th weight:
wV,j=3/(6π²)*[(U_j²−mV²)^(3/2)−(L_j²−mV²)^(3/2)], где U=min(Emax,E_j+h/2), L=max(mV,E_j−h/2).
Все event W масштабируются тем же wV,jnew/wV,jold, сохраняя nonnegative rates, на-shell E_j=k_i+k_l, C=Nν+Nantiν+2NV, L и E. Изолированно корректируется **нулевой момент E-bin**, не полный E-weighted integral. Источник контекста quantum kinetic: [Froustey/Pitrou/Volpe](https://arxiv.org/abs/2008.01074). Чужого кода не копировали.

## Машинные результаты, исходники и воспроизведение

- [Скрипт](../../../../kinetics/vector_neutrino_shell_radial_endpoints.py)
- [Три регрессионных теста](../../../../tests/test_vector_neutrino_shell_radial_endpoints.py)
- [Подробный аудиторский отчёт](../../../../kinetics/VECTOR-NEUTRINO-SHELL-RADIAL-ENDPOINTS-AUDIT.md)
- [Компактный JSON с полными h-N/E/S числами](../../../../kinetics/results/vector-neutrino-shell-radial-endpoints-2026-11-22.json)

Команды: python -m pytest -q tests/test_vector_neutrino_shell_radial_endpoints.py; python -m kinetics.vector_neutrino_shell_radial_endpoints.

Синтетические начальные fν=expit(−k+.15exp(−k)), fantiν=expit(−k−.1exp(−k/2)), b=.75 BE(E), mV=.4,Emax=4,h=.4,.2,.1,.05,T=Γ0=1, continuum independent p×cosθ 144×144. Seed не нужен. Start SHA закреплён, вторичные результаты в постоянном Git, Actions artifacts не единственное хранение.

## Фактические результаты и границы ошибок

На h=.2 исходная относительная энергия к независимому weak +2.2292%, верхняя трапеция −0.2765%, точно интегрированная физическая радиальная мера −0.2288%. На h=.1 исходная +1.1547%, трапеция −0.0980%, exact cell −0.0866%. На h=.05 исходная +0.5922%, трапеция −0.0342%, exact cell −0.0320% (улучшение в ~18.5 раза по абсолютной ошибке). Но число на h=.05 остаётся −0.1417% вместо исходных +0.1640%, а на h=.4 **ухудшается** с −0.712% до −2.694%. То есть улучшение одной меры не доказывает универсальную N/E/S точность. Угловое ∑W/(wV Γ0mV/E)−1 по-прежнему <1e−12. Conservation/entropy на bounded event step не изменены.

## Отрицательные случаи и неизвестная физика

Не представлена нулевая векторная мода p=0 и первый пороговый интервал; векторная энергия каждого целого cell по-прежнему коллоцирована E_j, ν k² measure и vector angular occupation piecewise constant. Отдельно остаются ошибки Emax/cutoff и dt, взаимодействия с конечной ν–e±γ плазмой, реальные гейдж-переносчик g/m(T), обе χ1/χ2, real φ Bose, H/entropy, freezeout и Ω/BEC/CLASS. Ни один процент **не является** ошибкой реликтовой плотности. CLASS-ready=false.

## Вывод, изменение плана, следующий checkpoint

Продемонстрировано, что значительная часть оставшегося weak E-bias была ошибкой boundary quadrature, а не нарушением Паули или энергии. Живой PLAN не меняет конечную цель, а уточняет численный приоритет: после правильного k-angular overlap брать finite E-bin clipped measure, затем независимую cutoff/h/dt экстраполяцию. Следующий шаг: проверить реальный CI нового коммита и двух предшественников, расширить физические форму f и mV/T, включить недостающий threshold shell и истинные energy-weighted cell moments без потери exact C/L/E и H; затем nonthermal ν–SM, χ1/χ2/φ/FRW. Все отрицательные данные сохранить.

## Поправки

Не отменяет [первый microreversible event graph](2026-11-22-vector-neutrino-reversible-event-graph.md) или [angular overlap stage](2026-11-22-vector-neutrino-shell-overlap.md), а уточняет их область точности; предыдущие CI-промахи остаются в архиве.
