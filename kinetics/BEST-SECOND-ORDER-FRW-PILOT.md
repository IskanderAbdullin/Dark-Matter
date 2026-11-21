# Проверка точности второго порядка для настоящей ранней FRW Bose эволюции: BEST exprb2

До сих пор независимый BEST v1.3.5 запускал метод exprb в предписанной RD Вселенной. [Реальная последовательность 4→8→16→32→64→128](BEST-RADIATION-FRW-PILOT.md) обнаружила **первый порядок по Δt**: даже при 128 шагах fφ(qmin) отличается от решения на 64 шагах на 0.508%. Это не означает, что кинетика неверная: first-order integrator требует слишком мелкий time grid для холодного IR.

Сам pinned BEST v1.3.5 уже содержит **собственный двухстадийный экспоненциальный метод exprb2** — это НЕ наша новая численная схема. По опубликованным исходникам и [описанию BEST](https://github.com/best-hep/best), он строит predictor exprb, повторно измеряет collision integral на end-of-step predicted state и добавляет φ₂-коррекцию к экспоненциальному шагу. Он требует два Vegas rate passes/step, но заявлен метод второго порядка при гладких rate и конечном dt. В отличие от diagonal exprb, его φ₂-коррекция не автоматически гарантирует положительные occupations; BEST имеет clamps, и каждый результат требует контроля born number, Pauli bounds и energy residual.

## План строго независимой диагностики

Физика совершенно идентична существующему реальному RD experiment:

- χ1(M/T0=.1, FD g=4) ↔ χ2(massless, FD g=4) + massive φ(m/T0=.001, BE g=1);
- одна закрытая обратимая 1↔2 реакция, Γrest=1, начальная mother FD и маленькая daughter occupation 1e−12 e−q;
- a(t)=sqrt(t/1e−4), t0=1e−4, tfinal=1.4e−4, afinal=sqrt(1.4);
- 56 log momentum bins, qmin=1e−4, qmax=30, Vegas neval=6000, 2 MPI ranks;
- **exprb2** при Nt=4,8,16,32, так что ΔtΓ=1e−5, 5e−6, 2.5e−6, 1.25e−6.

[best_radiation_frw_pilot.py](best_radiation_frw_pilot.py) теперь явно разрешает метод BEST exprb2 в протоколе и CLI, но не изменяет пакет BEST. [GitHub Actions](../.github/workflows/best-frw-massive.yml) исполнит все четыре случая и сохранит genuine JSON/negative logs. Сверим error orders для полного nφ, IR fφ(qmin), подсеточной massive zero-mode N0, стехиометрических зарядов и независимой physical FRW pressure-work quadrature.

Основные численные контрольные пределы от уже **выполненных** независимых exprb 64/128 и прямого аналитического open-bath FRW эталона:

- exprb Nt128: fφ(qmin)=1.02467504042, nφ=7.212682366e−7, N0=1.698379206e−14;
- Richardson первого порядка 64/128: fφ(qmin)≈1.019468120;
- точный континуальный *приближённый* redshifted collisionless mother source: fφ(qmin)=1.019683988, regular IR volume N≈1.709258419e−14.

Важно: exprb2 второй порядок **для численного интегратора**, а не для природы Monte Carlo source; если 4→8→16→32 не показывает разницы, убывающей примерно ×4, либо ограничения сетки/Vegas/finite qmin доминируют, либо неправильно обработана FRW time dependence. Не подгонять пороги под результат: honest negative evidence полезнее искусственно показанной сходимости.

Даже отличная temporal accuracy не заменяет реальных gstarS(T), SM bath Q, Γ(T) из лагранжиана, 86-эВ масс и начальной популяции, долгого физического freeze-out, поздней конденсации и CLASS P(k). Эти пункты по-прежнему НЕ выполнены.
