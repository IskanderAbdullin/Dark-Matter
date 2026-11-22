# Этап: 2026-11-22 — χ–ℓ Pauli logit-согласованная реконструкция возвращает H-теорему

- **ID/UTC:** DM-KINETICS-20261122-CHI-LEPTON-LOGIT-H-REPAIR. Дата 2026-11-22; точные UTC время начала/окончания не измерялись.
- **Направление / статус:** kinetics; аналитический метод + независимая локальная numerical проверка; новый Scientific CI репозитория требуется проверить.
- **Исходный SHA/ветка:** [05cac6ae](https://github.com/IskanderAbdullin/Dark-Matter/commit/05cac6ae4412f531188d991127bae921034ceb8f), research/kinetic-equilibrium-audit-20261121. Итоговый SHA — commit с этим stage.
- **Карты проекта:** [GOAL](../../../../GOAL.md) → [PROJECT_MAP](../../../../PROJECT_MAP.md) → [PLAN](../../../../PLAN.md) → [научные ограничения](../../../../kinetics/RESEARCH_GUARDRAILS.md) → [kinetics спецификация](../../../../kinetics/CONTINUOUS-RESEARCH-SPEC.md) → [NEXT_ACTION](../../../../kinetics/research-control/NEXT_ACTION.json).

## Заданный до опыта вопрос и критические альтернативы

[Предыдущий контрпример](2026-11-22-chi-lepton-conservative-entropy-negative.md) доказал, что положительная barycentric deposition, exact Nχ,Nℓ и Qχ+Qℓ **не гарантируют** H-теорему при arbitrary nonlinear Fermi distribution; при N10 bin Sdot отрицательна при истинной event Sdot>0. Исследовать строго согласованный альтернативный путь: реконструировать в каждой физической точке 2→2 event не f по прямой, а **logit(f_i)** теми же barycentric коэффициентами, которыми потом проектируется source. Критерии: точная положительность Sdot_bins=Sdot_event, exact number/energy, FD affine-logit test; при этом не маскировать физическую ошибку восстановления oscillatory shape.

## Физика, доказательство и применимость

Векторный portal χℓ↔χℓ с real φ U1'-compatible qχ1=qχ2, arbitrary nonthermal fχ(E), fℓ(ω), полный finite-mV, CM boost, Pauli F/R. Для каждой энергии x между grid nodes x_i,x_{i+1}: α_i=(x_{i+1}−x)/h, α_{i+1}=1−α_i; `ℓ(x)=Σα_i logit(f_i)`, `f_event(x)=expit[ℓ(x)]`. Тот же α применяется к forward/reverse deposited `dN_i` на двух видов χ и ℓ. Тогда каждое событие отдельно даёт `Sdot_bins=J[ℓχ(E)+ℓℓ(ω)−ℓχ(E′)−ℓℓ(ω′)]=Jln(F/R)≥0`. Никакой posterior entropy clipping нет; Pauli у reconstructed f гарантирована конструкцией. Входы f_i∈(0,1), границы f=0/1 пока `fail closed`.

Важно: это **instantaneous** bin rate projection при заданных node occupations, не доказательство того, что произвольный конечный timestep сохраняет bounds 0≤f_i≤1 при использовании этих `df_i/dt`. Для настоящего physical strong 2→2 operator ещё нужны phase-space cell measures, Fermi-bounded time integrator, matched weak/strong entropy/number/energy на сходящихся сетках и full SM/Z′/χ1χ2φ FRW.

## Файлы и воспроизводимый эксперимент

- [Исходный код](../../../../kinetics/paired_chi_lepton_logit_galerkin_5d.py)
- [Тесты](../../../../tests/test_paired_chi_lepton_logit_galerkin_5d.py)
- [Научный отчёт](../../../../kinetics/PAIRED-CHI-LEPTON-LOGIT-GALERKIN-5D-AUDIT.md)
- [Компактный первичный JSON](../../../../kinetics/results/paired-chi-lepton-logit-galerkin-2026-11-22.json)
- Команды: `python -m pytest -q tests/test_paired_chi_lepton_logit_galerkin_5d.py`; `python -m kinetics.paired_chi_lepton_logit_galerkin_5d`

NumPy/SciPy/pytest, нет RNG/seed. `Mχ=.1,mℓ=0,G=.001,mV=.1,cut=12,Te=1,nrad=4,nmu=4,ncos=12,nphi=12`, energy grids bins=10/20/40/80/160/320. Исходные физически допустимые synthetic `fχ=expit(-E+.3sin(4E))`, `fℓ=expit(-ω+.3cos(4ω))`. Независимый эталон — те же analytically prescribed f(E) **без** grid interpolation в непрерывном 5D event integrator; эксперимент фиксирует именно discretization shape bias, не физическую сеточную непрерывную сходимость исходного 5D источника.

## Измеренный положительный/отрицательный результат

При всех шести bins `Sdot_bin/Sdot_event−1` у *logit-consistent реконструированного* источника находится на уровне округления (обычно ≤4×10⁻¹⁴), Sdot>0, Nχ,Nℓ и Qχ+Qℓ сохраняются до машинной ошибки. FD с отдельными химическими потенциалами (аффинные logit) реконструируется между узлами точно.

**Не закрытый, но измеренный отрицательный физико-численный результат:** для original analytic oscillatory f continuous 5D `Sdot=+6.430366438e−13`, а новая logit grid Sdot при N10 = `+4.269218391e−14` (**−93.36%**), N20 `+8.019647525e−14` (−87.53%), N40 `+1.938153454e−13` (−69.86%), N80 `+5.285160749e−13` (−17.81%), N160 `+6.037665127e−13` (−6.11%), N320 `+6.341565574e−13` (−1.38%). Исправление H-теоремы **не** доказывает физическую точность coarse energy grids. При N10 Qχ также значительно занижен (1.89814e−14 против 2.55857e−13 оригинального 5D аналитического f). Нет 3+ входящих radial/CM/cutoff проверок.

## Научный вывод, изменение плана и точный следующий шаг

**Поддержано с указанной областью:** одновременное использование logit-interpolation и тех же conservation-preserving event deposition возвращает на дискретном уровне истинную микроскопическую quantum H-теорему без скрытых коррекций. **Остаётся отрицательным:** крупная shape alias при N10–80 и отсутствие time/Pauli bound, реального couplings и full physical χ1/χ2/φ+SM+FRW. Этот метод — кандидат для следующего solver, но не научная валидация физической 86-eV Bose DM.

**NEXT_ACTION:** проверить реальный GH CI, построить quantum-entropy-compatible finite-volume **временную** схему на node Fermi occupations с hard bounds и истинной плотностью состояний; сравнить с independent weak 5D for ≥3 energy/radial/cutoff and oscillatory shapes, затем coupled χ1+χ2+φ Bose/e±γν/Z′+FRW. Сохранять отсутствие подтверждённой Γel/ΩDM/BEC/CLASS. Все отрицательные результаты архивированы [EVIDENCE](../../../../kinetics/research-control/EVIDENCE.md) и [BACKLOG](../../../../kinetics/research-control/BACKLOG.md).

## Поправки

Продолжает и не стирает [более ранний falsifying barycentric stage](2026-11-22-chi-lepton-conservative-entropy-negative.md); новый entropy-safe метод устраняет **структурную** ошибку H в том же limited model, но оставляет большую ошибку спектра и абсолютных rates.
