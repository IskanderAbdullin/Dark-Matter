# Учёт реальной теплоёмкости SM-плазмы при упругой ретермализации материнской χ1

**После [числосохраняющего химического projection](SM-NUMBER-CONSERVING-BATH.md)** следующий узел модели — не считать окружающую SM plasma бесконечным неподвижным резервуаром температуры. При обмене кинетической энергией χ1↔SM мать может остывать/нагреваться, а её физическая энергия должна возвращаться в динамическую переменную \(\rho_\mathrm{SM}\), влияя на \(T_{SM}\).

Разработана независимая [sm_radiation_heat_feedback.py](sm_radiation_heat_feedback.py): **один фиксированный-a collision substep** с конечной, а не бесконечной теплоёмкостью *идеализированной радиационной* SM plasma. Этот этап не является полным FRW BEST+SM решателем: переменные тёмной реакции 1↔2, расширения и Hubble пока не интегрируются одновременно.

## 1. Уравнение энергии с постоянным gstar

Для ideal radiation с одним постоянным \(g_{\star\rho}=g_{\star S}=g_\star\), \(T_b>0\):
\[
\mathcal E_{SM}=a^4\rho_{SM}=\frac{\pi^2}{30}g_\star(aT_b)^4.
\]
Масса χ1 постоянна; упругий step [SM-NUMBER-CONSERVING-BATH.md](SM-NUMBER-CONSERVING-BATH.md) даёт **фактически вычисленное** изменение комовинговой энергии
\[
\Delta\mathcal E_{\chi_1}=
\frac{g_1}{2\pi^2}\sum_i w_iq_i^2
\sqrt{q_i^2+(aM_1)^2}\,[f_{\chi_1,i}^{new}-f_{\chi_1,i}^{old}].
\]
Тогда единственный допустимый тепловой баланс collision substep:
\[
\boxed{
\mathcal E_{SM}^{new}
=\mathcal E_{SM}^{old}-\Delta\mathcal E_{\chi_1}},
\qquad
\boxed{
T_b^{new}=
\left[
\frac{\mathcal E_{SM}^{new}}
{(\pi^2/30)g_\star a^4}
\right]^{1/4}}.
\]

То есть **тёплая** χ1 при cooling отдаёт энергию SM и повышает Tb, а слишком **холодная** χ1 при thermalization отбирает энергию из SM. Скорость обмена определяет пока **феноменологическая** \(\gamma_{elastic}\), а не микрофизический матричный элемент. Степень свободы gstar=55 здесь пример constant relativistic bath, а не реальная параметризация всех SM thresholds.

## 2. Невозрастающая свободная энергия / энтропия

Новая физическая переменная позволяет независимо оценить *полный* комовинговый энтропийный ledger:
\[
S_{SM}=a^3\frac{2\pi^2}{45}g_\star T_b^3,
\]
\[
S_{\chi_1}=-\frac{g_1}{2\pi^2}\int dq\,q^2
\{f\ln f+(1-f)\ln(1-f)\}.
\]
Удобный test case: M=.1, a=1.18, initial Tb=.82, gstar=55, Γelastic=1.2, dt=.4 и слегка неравновесный FD parent. Реальный conservative step по проектируемой FD temperature/chemical µ демонстрирует:
- **число χ1 не меняется** по той же сетке даже при конечной температурной перенастройке;
- dark transfer \(\Delta\mathcal E_{\chi1}\approx-0.01308134\), энергия \(\mathcal E_{SM}\) повышается ровно на противоположное значение;
- temperature bath повышается с .82 до **~.8201690**;
- total entropy \(S_{SM}+S_{\chi_1}\) возрастает на **~0.00036333**.

Метод **явно отклоняет** too-large nonphysical steps, если суммарная энтропия падает из-за нарушения допущения о постоянной Tb на протяжении субшага, вместо объявления такой динамики корректной.

## 3. Что уже завершено и что нужно затем

Прототип построен поверх **точного positive/Pauli-safe projected elastic substep**, поэтому для любого конечного γΔt сохраняет \(0\le f_{\chi_1}\le1\) и χ1 number. При фиксированном a SM heat capacity / temperature обновляются из измеренного Δenergy, без скрытого source. Специальные regression tests проверяют cooling χ1→SM, heating SM→χ1, полный energy closure, total entropy, stiff Γ, invalid inputs.

До полноценной космологической модели нужны:
- отдельный истинный \(C_{1\leftrightarrow2}\) BEST + finite SM thermal contact с правильным split/implicit energy feedback и QS/QH balance;
- varying \(g_{\star \rho}(T)\), \(g_{\star S}(T)\), быстрые SM scattering и physically derived \(\gamma_{elastic}(T)\), включая thermal masses;
- давление при расширении \(d(a^4\rho)/d\ln a=a^4(\rho-3P)\), а не ложное постоянство энергии массивной χ1;
- Hubble из настоящих \(\rho_{\rm tot}\), либо количественно обоснованное neglected dark backreaction;
- только после этого настоящий 86-eV DM production, Bose burst, freezeout f(q), hot tail и CLASS P(k).

**Итог текущего строго ограниченного шага:** полноценный **физически учитываемый энергообмен матери↔идеализированной плазмы на отдельном collision substep** создан; **космологического SM+dark coupled решения пока нет**. Нельзя смешивать этот факт с готовностью relic.

Код: [sm_radiation_heat_feedback.py](sm_radiation_heat_feedback.py). Тесты: [test_sm_radiation_heat_feedback.py](../tests/test_sm_radiation_heat_feedback.py). Выполнение примера: python kinetics/sm_radiation_heat_feedback.py.
