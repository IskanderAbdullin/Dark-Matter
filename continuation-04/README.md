# Продолжение IV: кинетическая устойчивость, примордиальная фазовая функция, инварианты Власова и конечные Bose halos

**Дата:** 2026-11-21. **Исходный источник:** I. G. Abdullin, V. A. Popov, *Boson dark matter halos with a dominant noncondensed component*, JCAP 11 (2021) 055, https://doi.org/10.1088/1475-7516/2021/11/055. Использован повторно открытый исходный файл \`rv2_JCAP_041P_0721.tex\` (уравнения distribBE, dispV0, Bose-gas, inhR, DMabundance, mLowerBEC).

**Основной выбор:** вместо новых подгонок продолжена математическая проверка исходной теории и её происхождения из микрофизики. Нет никаких GitHub Actions, пайплайнов или real-money вычислений. Этот пакет в отдельной ветке \`research/bose-halos-theory-2026-11-21\`, без изменений \`main\`.

## Новые независимые вычисления
1. [01 — Кинетический критерий Джинса](papers/01-kinetic-jeans-of-near-critical-bose-gas.md). Аналитически линеаризованы Vlasov–Poisson и repulsive Hartree–Fock; получен точный *локальный, однородный, при \(\omega=0\)* предел \(k_J^2 R^2=A/(1+bA)\), \(A=\operatorname{Li}_{1/2}(z)/\operatorname{Li}_{3/2}(z)\). Для печатных F583-4 параметров \(\nu=2.61,R=0.6\) кпк: \(\lambda_J\approx118.47\) пк; это не подтверждение нестабильности всей сферической галактики.
2. [02 — Замороженный релятивистский реликт ≠ поздний нерелятивистский Bose-газ](papers/02-adiabatic-relic-and-abundance.md). Строго сравниваются функции \(f(p)\), состояния после декуплинга и теорема об отсутствии спонтанного BEC при *только* однородном адиабатическом расширении уже нерелятивистского свободного газа. Проверена зависимость коэффициента космологического изобилия от определения \(g_{*S,0}\) и числа степеней свободы.
3. [03 — Шестимерные Casimir-уровни](papers/03-six-dimensional-casimir-tail.md). Высокозаполненные области у галактического идеального критического центра имеют объем фазового пространства \(\propto F^{-3}\), **как** у однородного релятивистски замороженного газа, хотя распределения по \(p\) отличаются. Выведен *ненулевой* подведущий член \(\propto F^{-7/2}\), отличающий полные спектры уровня заполнения в закрытой collisionless постановке.
4. [04 — Конденсат с конечной массой при *точно критическом* центре](papers/04-critical-lowered-bose-casimir.md). Ранее при \(\psi_0<\eta\) получены четыре профиля. Теперь при \(\eta=\psi_0=4\) найдено конечномассовое решение \(\xi_t=40.299, Q_t=26.357\), получена новая центральная линейная асимптотика \(x=1-0.589810\xi+O(\xi^2)\). Построена обратная выпуклая Casimir-функция для lowered-Bose-анзаца, **без заявления о доказательстве полной устойчивости**.
5. [05 — Тидальная выживаемость Bose stars](papers/05-bosestar-tidal-survival.md). Радиус из уравнений GP/TF связан с требованием на число частиц, но даже TF допустимый объект может быть разрушен приливным полем. Для условного примера \(R_c\simeq133.45\) a.u., \(V_c=44.1\) км/с, \(r=1\) кпк, \(M_{\rm tidal}\gtrsim4.87\times10^{20}\) кг (только для сферического приближения flat-rotation и круговой орбиты).
6. [06 — Карта источников, ошибок модели, проверок](literature/prior-art-and-next-tests.md) + [журнал исследований](RESEARCH-LOG.md).

## Степени достоверности
**Доказано в названных упрощенных уравнениях:** идентичность полилогарифмического отклика, стационарное BE vs frozen-rel спектральное различие, Liouville phase-volume асимптотика, \(df/dE<0\) и выпуклость обратной Casimir, геометрическая граница по приливному ускорению.

**Выполнено численно:** momentum quadrature независима от polylog (относительная ошибка \(\sim2\times10^{-16}\)); imaginary-frequency Vlasov dispersion, finite-radius critical lowered Bose ODE, 6-D Casimir integral, физические размерности. JSON-результаты находятся в \`results/\`.

**Не выполнено:** глобальные собственные моды сферического гало, нестационарная кинетическая симуляция GP–Boltzmann, CLASS/CAMB, проверка на независимых галактиках с nuisance, реальное образование и функция масс Bose stars. Ни один вычисленный масштаб здесь не является наблюдательным исключением массы 86 эВ.

**Строгие ограничения:** результаты Vlasov–Jeans — **однородная локальная** модель (не глобальная); Casimir-объем — **fine-grained и замкнутое множество частиц** (не coarse-grained galaxy fit); physical cross section — **иллюстративный верхний ориентир**, не измерение; исходные \(R,V_0,\nu,m\) приведены в статье с округлением и не все значения строго алгебраически независимы.

## Воспроизведение без облачного CI
\`python -m pip install -r requirements.txt\`; \`python code/kinetic_jeans.py\`, \`python code/kinetic_growth.py\`, \`python code/frozen_relic.py\`, \`python code/casimir_ir_tail.py\`, \`python code/casimir_mass_anchor.py\`, \`python code/critical_lowered.py\`, \`python code/bosestar_tides.py\`.
Численные результаты локально проверены до помещения в GitHub; см. README, записи точности/модели и требования к пакетам.

### Дополнительные результаты VII–VIII
- [06](papers/06-jeans-wkb-validity.md): строгая проверка применимости однородного критерия Джинса; центральный высокооткликающий слой существенно меньше соответствующей длины волны.
- [07](papers/07-lowered-bose-virial.md): точное тождество для кинетического давления, независимая численная проверка теоремы вириала с относительной невязкой 4.13×10⁻¹¹.
- [08](papers/08-universal-wkb-asymptotic.md): аналитическая универсальная асимптотика λ_J/r_half∝(ν_crit−ν)^{-1/2}; подтверждена точными корнями полилогарифмов.
- [literature/prior-art-and-next-tests.md](literature/prior-art-and-next-tests.md), [NEXT-DECISIVE-QUESTIONS.md](NEXT-DECISIVE-QUESTIONS.md), [RESEARCH-LOG.md](RESEARCH-LOG.md): источники, план фальсифицирующих проверок и исследовательский журнал.

- [09](papers/09-critical-cutoff-family.md): строгий точный предел от конечномассового критического усечённого Bose-гало к бесконечной изотермической оболочке; 7 независимых численных профилей. Дополнительно **доказана строгая монотонность** центральной квантовой нормировки по параметру отсечения, но не общая монотонность внешнего радиуса и массы.

## Additional results IX–XII

- [09: critical cutoff family](papers/09-critical-cutoff-family.md): rigorous monotonicity of central phase normalization, seven numerical critical equilibria.
- [10: Casimir mixing constraints](papers/10-casimir-mixing-constraints.md): Dehnen excess-mass necessary condition for a very specific closed thermal relic patch. Stronger finite-occupation mismatch and an exact-support no-go; NOT a physical particle exclusion.
- [11: universal shallow-cutoff scaling](papers/11-shallow-cutoff-universal-scaling.md): a parameter-free limiting ODE with r_t/R approaching 11.860598 sqrt(eta), and Q_t approaching 2.895576 eta^(3/2).
- [12: smoothly lowered Bose](papers/12-smooth-quadratic-lowered-bose.md): convex tangent-subtracted positive distribution, continuously differentiable at the escape boundary; a calculated edge exponent 7/2, self-consistent finite mass, and virial test.
- [Consolidated Russian report](CONSOLIDATED-RESULTS-RU.md) and [recent mixing/stability literature](literature/new-mixing-and-stability-sources.md).

Reproduce using the Python scripts and JSON under code/ and results/. No GitHub Actions were used or requested.

## Two additional infrared mathematical results

- [13: finite-cutoff critical occupation Lp classes](papers/13-critical-ir-lp-integrability.md): exact critical six-dimensional Bose occupation is Lp for p<3 but not for p>=3 (finite outer cutoff). Central local momentum integral of f² is singular, yet the full local six-dimensional f² integral is finite. Distinguish this from actual collision rates.
- [14: weak-repulsion Hartree–Fock tail](papers/14-hartree-fock-ir-phase-volume.md): fixed positive coupling and exact critical center yield an IR phase-volume F^(-9/4) rather than ideal F^(-3), with different Lp threshold 9/4. This is a **nonuniform formal limit** reached only at enormous occupations for the illustrative physical coupling and cannot be treated as a real observable prediction without better quantum critical theory.
- [Reproducibility script](code/ir_phase_volume_exponents.py) and [stored numerical check](results/ir_phase_volume_exponents.json): high-precision implicit EOS verification and phase-volume scaling.

These two appendices bring the mathematical notes in this continuation to **14**. All results and caveats are documented; no GitHub Actions were run.
