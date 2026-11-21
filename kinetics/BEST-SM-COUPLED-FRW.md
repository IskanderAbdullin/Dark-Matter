# Первая потенциально замкнутая BEST χ1↔χ2+massive φ + SM radiation модель с явным Q и T_SM

**Цель:** перестать выбирать между (a) правильной, но почти collisionless матерью в предыдущих реальных BEST RD runs и (b) искусственно поддерживаемой термальной FD матерью с бесконечным неучтённым SM резервуаром. Новый [best_sm_coupled_frw.py](best_sm_coupled_frw.py) **в одном time evolution** соединяет физическую квантовую χ1↔χ2+massive φ реакцию независимого pinned BEST и конечную SM radiation heat-capacity bath, которая получает ровно противоположную энергию от number-conserving elastic χ1 thermal contact.

## Что действительно будет интегрироваться

- BEST v1.3.5 pinned Git SHA c70875be714e7a6bc6d38578c34bb0dc9f5aefa9, только две MPI задачи. RD фон \(a(t)=\sqrt{t/10^{-4}}\), tΓ=0→4×10⁻⁵, a_end=sqrt(1.4).
- χ1 массивна M/T0=.1, FD g=4; χ2 massless FD g=4; φ массивна mφ/T0=.001 BE g=1, имеет отдельный BEST zero-mode; химия χ1↔χ2+φ обратима и содержит Bose stimulation, Pauli blocking, true on-shell 1→2 kinematics.
- Реальное up to date \(f_{\chi1}(q)\) не сбрасывается на FD(T) насильно. Прежде каждого elastic substep подбирается \(\mu_{\chi1}\) так, чтобы FD(T_SM, µχ1) имела **то же конечносеточное комовинговое Nχ1**, что текущая мать.
- Эффективная elastic γ_el≥0 задаёт exact convex exponential step между fχ1 и number-projected FD. Он сохраняет Nχ1, 0≤f≤1 и уменьшает Fermi KL.
- Изменение comoving dark energy непосредственно передаётся \(\Delta\mathcal E_{SM}=-\Delta\mathcal E_{\chi1}\); \(T_{SM}\) затем вычисляется из \(\mathcal E_{SM}=g_{\star}(\pi²/30)(aT_{SM})^4\). Это уже не «бесплатный» thermal reset; gstar=55 принято постоянным и одинаковым для entropy и energy.
- Разрешены Γelastic=0 (контроль: полностью collisionless мать) и большая γelastic/H (быстрое упругое перераспределение). Глобальный integrator использует **проверенный реальным BEST метод второго порядка exprb2** внутри трёхчастичной химии.

## Как разделён шаг

На каждом временном интервале Δt:
1. При текущем a, TSM выполнить **половину** elastic χ1⇄SM BGK step и точно отразить heat Q в конечной радиационной теплоёмкости.
2. Продвинуть через настоящий BEST полный χ1↔χ2+φ интервал Δt, одновременно меняя a и перемещая все три species/массовую zero-mode в комовинговом пространстве.
3. При новом a и обновлённом из comoving SM energy значении TSM выполнить вторую половину elastic step и обновить BEST interpolatorχ1 **на каждом MPI ранге**.

Это Strang-подобный operator splitting; само условие симметрии не делает автоматически полную схему строго второго порядка при каждом Monte Carlo source. При вычислении после каждого шага проверяется отдельный physical FRW ledger:
\[
\boxed{\frac{d}{d\ln a}\left[a^4(\rho_{\rm dark}+\rho_{SM})\right]
=a^4\left[(\rho_{\rm dark}+\rho_{SM})-3(P_{\rm dark}+P_{SM})\right].}
\]
SM radiation имеет \(P_{SM}=\rho_{SM}/3\) и нулевой pressure-work contribution. Каждая упругая часть имеет **нулевую** общую переданную энергию в том же a; сама χ1↔χ2+φ реакция сохраняет суммарную физическую энергию до разрешения BEST, но не подгоняется к точному нулю.

Два числа \(N_{\chi1}+N_\phi\) и \(N_{\chi1}+N_{\chi2}\) не меняются даже от elastic contact, потому что он **упругий** и сам Nχ1 сохраняет. При химическом номер-изменяющем SM канале потребуются новые источники в этих уравнениях — недопустимо объявлять их консервативными тогда.

## Обязательные тестовые сценарии и критерии

- **Gamma_el=0**: реальный BEST+SM должен воспроизвести уже проведённый closed FRW experiment при той же сетке, времени и source. При этом heat transfer Q≈0 и comoving SM energy постоянна; физически TSM∝1/a.
- **Gamma_el=100000, H0=5000**: оценить Γelastic/H≈20 в начале, сравнить изменение x/Tχ, \(\mu_{\chi1}\), числа/энергии/soft fφ/zero-mode с нулевой γelastic.
- **Nt=4 vs8 or16**: честная временная чувствительность, проверка отдельных split heat ledgers и Fermi KL/total entropy, двух numbers, полной физической FRW pressure-work residual. Даже если изменение soft fφ чрезвычайно мало, это должно быть численно измерено.
- Отсутствие отрицательных populations, нереального нагрева SM, нет скрытого number creation или бесплатной rethermalization.

## Что эта версия не может подтверждать

Это **впервые честно связанный феноменологический** BEST + finite-radiation temperature feedback candidate, но всё ещё **не физически предсказанный thermal-plasma механизм Yin (2023)**: Γelastic не выведена из SM-χ1 scattering matrix element, number-changing χ1 replenishment отсутствует, gstarρ(T) и gstarS(T) постоянны, \(\Gamma_{1↔2}(T)\) условна, Hubble задан, а не выведен из реальных полных densities. Нет длительного freeze-out, доказанной BEC, mass model 86 eV, CLASS, observational likelihood. Физические выводы должны оставаться условными до всех этих проверок.

Код: [best_sm_coupled_frw.py](best_sm_coupled_frw.py), компоненты [number_conserving_sm_bath.py](number_conserving_sm_bath.py) и [sm_radiation_heat_feedback.py](sm_radiation_heat_feedback.py). Обычная regression suite ничего не подделывает: без --execute она строит только protocol, а при отсутствующем pinned BEST отказывает.
