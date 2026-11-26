# Этап: 2026-11-22 — консервативный χ–ℓ 2→2 binned source нарушает H-теорему

- **ID/UTC:** DM-KINETICS-20261122-CHI-LEPTON-ENTROPY-PROJECTION-NEGATIVE. Дата 2026-11-22, начало/конец UTC не измерялось.
- **Направление/статус:** kinetics, scoped analytic proof and locally reproducible numerical counterexample; full Scientific CI нового HEAD требуется проверить.
- **Ветка/исходный SHA:** research/kinetic-equilibrium-audit-20261121, [b4d6549](https://github.com/IskanderAbdullin/Dark-Matter/commit/b4d6549504f60f8dbf87658db2f3df0e7aaa4ca7). Итоговый commit будет содержать этот stage.
- **Контекст:** [общий GOAL](../../../../GOAL.md) → [PROJECT_MAP](../../../../PROJECT_MAP.md) → [PLAN](../../../../PLAN.md) → [kinetics научное ТЗ](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md) → [RESEARCH_GUARDRAILS](../../../../kinetics/RESEARCH_GUARDRAILS.md) → [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Научный вопрос и критерий до испытания

Построить первый paired weak-to-binned источник из полного 5D nonthermal finite-mV χℓ→χℓ Pauli event, исходные/конечные энергии χ и лептона размещать на энергетической сетке положительной barycentric линейной интерполяцией. До расчёта: обязана сохраняться масса/число каждого вида `sum dNχ=0`, `sum dNℓ=0` и энергия `ΣEχ dNχ=−ΣEℓ dNℓ=Qχ_weak`; отдельно проверить знак дискретной `Sdot_bins` против **истинного** `Sdot_event≥0`. Конкурирующие предположения: linear moment conservation может быть недостаточным для entropy; cell-interpolation log-odds может разрушить H theorem. Фиксировать отрицательный результат, не добавлять clipping.

## Микроскопическая модель и независимый метод

[Существующий локально подготовленный χℓ weak kernel](../../../../kinetics/nonthermal_chi_lepton_5d_weak.py) рассчитывает `J=F−R` с настоящей кинематикой `E+ω=E′+ω′`, spin-averaged finite-vector t-channel матрицей и произвольными строго внутренними fermion fχ(E),fℓ(ω). Для каждого события с весом W переносится `WJ` из входных в выходные χ и ℓ узлы. Приняты неотрицательные линейные коэффициенты `α_i=(E_{i+1}−E)/h`, `α_{i+1}=(E−E_i)/h`. Тождества `Σα_i=1`, `ΣE_iα_i=E` гарантируют точные N и E **в дискретной проекции** независимо от формы fχ,fℓ.

Истинная микроскопическая энтропия `Σ W(F−R)ln(F/R)≥0`. Для предполагаемого следующего сильного оператора naïve bin энтропия `Sdot_bins=−Σ dNbin ln(f(Ebin)/(1−f(Ebin)))` — это интерполяция log-odds. Для термального FD `ln[f/(1−f)]=(μ−E)/T` является линейной функцией энергии, поэтому сравнение совпадает по моментам. Для nonlinear nonthermal log-odds **математической гарантии положительности нет**. Это конкретный контрпример неполной физической strong reconstruction, не отказ микроскопической H theorem.

Физический model scope: χ Dirac Mχ=.1 (условные единицы), massless lepton mℓ=0, finite vector `G=.001,mV=.1`; angular/azimuth 12×12, incoming p,k,μ Gauss 4×4×4, kinetic E cutoff12, gχgℓ=4×4 историческая численная нормировка, no verified UV/lab couplings. Литература контекста: [численная Boltzmann positivity и entropy](https://www.sciencedirect.com/science/article/pii/S0021999103000652). Код не копировали из сторонних публикаций.

## Данные, команды, собственные измерения и ошибки

- [Полный executable](../../../../kinetics/paired_chi_lepton_barycentric_5d.py), [регрессии](../../../../tests/test_paired_chi_lepton_barycentric_5d.py).
- [Подробная независимая научная заметка](../../../../kinetics/PAIRED-CHI-LEPTON-BARYCENTRIC-5D-AUDIT.md).
- [Машинный JSON: пять/шесть сеток и отрицательный случай](../../../../kinetics/results/paired-chi-lepton-barycentric-5d-2026-11-22.json).
- Команды: `python -m pytest -q tests/test_paired_chi_lepton_barycentric_5d.py` и `python -m kinetics.paired_chi_lepton_barycentric_5d`.
- Python NumPy/SciPy/pytest, seed нет. Синтетические fermion occupations, ни наблюдательные, ни физически эволюционированные данные.

**Контрпример, заданный до проверки H:** `fχ(E)=expit[-E+0.3 sin(4E)]`, `fℓ(ω)=expit[-ω+0.3 cos(4ω)]`, всегда 0<f<1. Микроскопическая `Sdot_event=+6.43036644×10^-13`. Число и энергия χ/ℓ сохраняются на ВСЕХ выбранных bins до абсолютной ошибки ~1e−26 и ~1e−28. При N=10 интервалов `Sdot_bins=−3.60494354×10^-14`, при N=20 `−2.76287023×10^-15` — **нефизический отрицательный энтропийный поток**. При N40/80/160/320 `Sdot_bins=3.01614948e−13/5.82988043e−13/6.23071955e−13/6.38529985e−13`; при N320 относительная ошибка всё ещё −0.70%.

Для гладкого `fχ=expit(-E/.9+.15e^-E/.9)`, `fℓ=expit(-ω+.05e^-ω)`, микроскопический `Sdot_event=+9.7847996053e−15`, bin-S ошибки N10,N20,N40,N80,N160 = +0.861%, +3.260%, +0.0744%, +0.1559%, +0.00449% — **немонотонная** сеточная картина. При FD с affine log-odds bin-entropy совпадает с event энтропией до округления, отдельный контроль.

## Отрицательный вывод, принятый пересмотр плана и исполнимый шаг

Наивная энерго-сохраняющая barycentric проекция физического 2→2 quantum weak source **не является** H-theorem-совместимой нелинейной strong collision scheme. Её нельзя непосредственно интегрировать по времени с полноценными Pauli bounds и объявлять физической χ gamma. Ни отрицательный bin entropy, ни точная energy conservation сами по себе не опровергают 86-eV boson mechanism — это **дискретизационная ошибка конкретного метода**. Сохраняется отсутствие actual χ2 daughter f, ν/e±γ/Z′ bath, χ1↔χ2+real φ Bose, FRW, cutoff convergence, BEC/Ω/CLASS.

Следующий **P0**: вывести microreversible entropy-consistent paired `Cχ,Cℓ` с истинными фазовыми мерами, occupation-to-event energy interpolation и совместными законами N/E/S; проверить H theorem при осциллирующих nonthermal f и переход в FD, выполнить независимый comparison с continuous weak 5D на ≥3 energy/angular/cutoff сетках, а затем ввести обе χ species/real φ Bose с настоящим SM+FRW. [STATE](../../../../kinetics/research-control/STATE.json) и [EVIDENCE](../../../../kinetics/research-control/EVIDENCE.md) обновлены без удаления отрицательных результатов. `CLASS-ready=false`.

## Поправки и связь с предыдущим этапом

Это необходимое **уточнение/ограничение** предыдущего [χℓ nonthermal weak stage](2026-11-22-chi-lepton-nonthermal-weak.md). Он верен как интегральный функционал, но его моментная проекция на произвольный Pauli-неравновесный спектр не сохраняет H-теорему автоматически. Активный Git SHA хранит оба исходника и оба машинных JSON.
