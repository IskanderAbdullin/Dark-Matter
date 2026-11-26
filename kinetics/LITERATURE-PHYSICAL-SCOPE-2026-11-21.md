# Литературная проверка физической области применимости: Yin (2023), Sakurai–Yin (2025), BEST (2026)

**Дата систематизации: 2026-11-21.** Цель — не перепутать реально проверенный quantum collision kernel и динамику игрушечной closed χ1↔χ2+massive φ сети с **предсказанием холодной еВ тёмной материи**. Первые опубликованные работы мотивируют механизм и накладывают дополнительные физические обязательства, которых в текущем BEST benchmark ещё нет.

## 1. Первоисточники

1. **Wen Yin, Thermal production of cold “hot dark matter” around eV**, *Journal of High Energy Physics* 05 (2023) 180; DOI [10.1007/JHEP05(2023)180](https://doi.org/10.1007/JHEP05(2023)180); preprint [arXiv:2301.08735](https://arxiv.org/abs/2301.08735). Формулирует cold production of bosonic eV-range DM вследствие Bose enhancement at very low momentum в распадах/обратных распадах 1↔2 **из горячей thermal plasma**. Явно обсуждает: ignition, burst, saturation/quasiequilibrium, thermalization after burst; radiation-dominated background, relic abundance, structure-formation/free-streaming and hot subcomponents.
2. **Kodai Sakurai, Wen Yin, Stimulated emission of dark matter via thermal scattering: novel limits for freeze-in and eV cold dark matter**, *JHEP* 03 (2025) 202; DOI [10.1007/JHEP03(2025)202](https://doi.org/10.1007/JHEP03(2025)202); [arXiv:2410.18968](https://arxiv.org/abs/2410.18968). Показывает, что small-angle 2→2 scattering с Bose final state тоже вызывает stimulated emission, т.е. механизм шире простого 1↔2 распада, а abundance и phase-space distribution стандартного freeze-in могут существенно измениться. Обсуждает влияние массы DM на стимуляцию.
3. **Jong-Hyun Yoon, Boltzmann Equation Solver for Thermalization**, *Computer Physics Communications* **327** (2026) 110295; DOI [10.1016/j.cpc.2026.110295](https://doi.org/10.1016/j.cpc.2026.110295); [arXiv:2603.28848](https://arxiv.org/abs/2603.28848); [BEST GitHub](https://github.com/best-hep/best). Внешняя проверяемая реализация momentum-resolved quantum Boltzmann for n↔m scattering with BE/FD statistics, massive bosons and zero-mode, comoving momenta, RD cosmology and time integrators. В данной ветке использован строго Git SHA **c70875be714e7a6bc6d38578c34bb0dc9f5aefa9** и MPI из двух процессов.

## 2. Что реально проверено нашей моделью

- **Квантовая микрофизика 1↔2:** независимый BEST против аналитического on-shell integral воспроизвёл number/energy source, разные chemical affinities, Bose gain, Pauli blocking, inverse χ2+φ→χ1, finite mφ и реальный q→0 occupation. Измерены не только общий number, но и дифференциальный спектр и энергетические невязки.
- **Три динамических particle species:** χ1, χ2 и massive φ действительно меняют occupations во времени. Отдельная massive Bose zero-mode N0 обновляется, сохраняются две реакции-стехиометрические комбинации числа, проверяется physical energy.
- **Предписанное RD расширение:** независимый BEST использует a(t)=sqrt(t/t0), не подменяя comoving mother FD новой thermal FD на каждом шаге. Проверен pressure-work закон для комовинговой энергии; t-grid=4→128 показал first-order ошибки мягкого спектра, а не ошибочно объявил холодный DM готовым.
- **Отрицательный результат о конденсате:** при tΓ=4×10⁻⁵, H=0 получено N0∝qmin³ на трёх независимых qmin. Это regular low-q phase-space particle number, не макроскопический δ³(p) condensate. В предписанном RD отдельный такой тест готовится и верифицируется.
- **Независимый континуальный FRW ранний эталон:** frozen collisionless mother в комовинговом импульсе и source \(A(q,t)(1+f)\) дают fφ(qmin)≈1.019684; extrapolation реального внешнего BEST 64/128≈1.019468. Это свидетельство правильного красного смещения, но не физически полной плазмы.

## 3. Главное физическое расхождение с моделью thermal plasma Yin (2023)

В текущем FRW эксперименте мать χ1 стартует термальной, но затем **свободно красносмещается**, за исключением самого 1↔2. Такая постановка означает, что её температура/химпотенциал не поддерживаются внешней SM plasma. В модели cold-thermal production мать может эффективно **перетермализовываться** за счёт других, быстрых процессов: это уже физически другая система с собственными \(C_{\rm SM\leftrightarrow\chi_1}\), \(\dot{\rho}_{SM}=-Q\), температурой \(T_{SM}(t)\), gstarS(T), а при выраженной backreaction и Friedmann H(t) из полной энергии.

**Не разрешается** просто перезаписывать fχ1 на thermal FD после каждого шага и продолжать заявлять, что dark sector закрыт: это создаёт/уничтожает энергию и число из необъявленного резервуара. Важно строить два контролируемых физических сценария: (a) collisionless-mother исходное приближение и (b) реально thermalized mother с explicit SM exchange, затем исследовать, в каких параметрах применяется каждый.

Поэтому отличная сходимость недавнего thermal SNAPSHOT и collisionless FRW pilot не означает уже подтверждённого burst regime Yin. Для него важна **относительная скорость ретермализации матери к H и Γ**, физическая thermal mass, phase-space Bose growth и обратные распады вплоть до насыщения.

## 4. Дополнительные механизмы и строгая программа валидации

В работе Sakurai–Yin 2025 механизм soft stimulated emission также возникает в **2→2 small-angle scattering**. Для конкретного лагранжиана этот канал может конкурировать с 1↔2 и радикально изменить как количество, так и форму fφ(q) и сопутствующий hot component. Его нельзя автоматически исключить, не оценив thermal rates/couplings и реальные пороги.

Путь к физической 86-эВ гипотезе требует:

1. Выбрать конкретные \(M_{\chi_1},m_{\chi_2},m_\phi=86\,\mathrm{eV}\), dof, spin, запрещённые/разрешённые каналы, лагранжиан и настоящие \(|\mathcal M|^2\), \(\Gamma(T)\) и thermal effective masses; здесь mφ/T0=.001 и M1/T0=.1 пока просто экспериментальные.
2. Включить реальные plasma interactions или явную физическую аргументацию collisionless χ1 (rates vs H/Γ), plus SM entropy/energy ledger с переменными gstarS и gstarρ; не вводить скрытых внешних sources.
3. Исследовать нелинейную эволюцию **от ignition до saturation** и long-time freeze-out при достаточно тонких dt, Nq, Vegas, qmin, qmax; проверить physical zero-mode n0 и отдельно hot-tail.
4. По физическому лагранжиану проверить значимость и/или добавить стимулированные 2→2 каналы; их влияние на soft и thermal populations нельзя переносить из 1↔2 автоматически.
5. Получить абсолютное число bosonic relic с правильной нормировкой и реалистичный замороженный fφ(q), а также нулевую моду, если она действительно survives qmin→0; **без произвольной перенормировки на ΩDM**.
6. Только затем давать этот достоверный primordial phase-space distribution в CLASS или аналог, исследовать P(k), free streaming, Lyman-α/галактики/CMB и сравнивать с наблюдениями. Старый 86-эВ thermal-only benchmark является горячим и не доказывает cold mechanism.

## 5. Приоритет следующей физической задачи

Сразу после закрепления FRW qmin и exprb2 precision **самый важный нечисленный пробел** — конкретный SM bath rethermalization collision operator и строгий контроль передачи энергии/энтропии. Не стоит объявлять улучшение точности временной сетки доказательством физического thermal plasma source без этого шага.

Ссылки на собственные исследования: [BEST-CHEMICAL-AFFINITY.md](BEST-CHEMICAL-AFFINITY.md), [BEST-MASSIVE-DAUGHTER-REFERENCE.md](BEST-MASSIVE-DAUGHTER-REFERENCE.md), [BEST-CLOSED-MASSIVE-DYNAMIC.md](BEST-CLOSED-MASSIVE-DYNAMIC.md), [BEST-RADIATION-FRW-PILOT.md](BEST-RADIATION-FRW-PILOT.md), [RD-COLLISIONLESS-BOSE-REFERENCE.md](RD-COLLISIONLESS-BOSE-REFERENCE.md).
