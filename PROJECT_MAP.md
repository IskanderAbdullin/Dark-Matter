# Dark-Matter — карта научного контекста

**Это маршрутизатор:** при конкретной задаче открывать только нужное. Главная научная цель находится в [GOAL.md](GOAL.md); актуальный маршрут — в [PLAN.md](PLAN.md).

| Работа | Что читать после карты |
|---|---|
| **Общая цель и критерии научного завершения** | [GOAL.md](GOAL.md); подробная программа в [plans/RESEARCH_ROADMAP_DETAIL.md](plans/RESEARCH_ROADMAP_DETAIL.md) |
| **Космологическая кинетика / eV-бозон** | [kinetics/GOAL](kinetics/research-control/GOAL.md) → [kinetics/CONTINUOUS-RESEARCH-SPEC](kinetics/CONTINUOUS-RESEARCH-SPEC.md); при оценке результатов [научные инварианты](kinetics/RESEARCH_GUARDRAILS.md) |
| **Продолжение уже начатого кинетического исследования** | [STATE](kinetics/research-control/STATE.json) → [NEXT_ACTION](kinetics/research-control/NEXT_ACTION.json); только нужные записи [EVIDENCE](kinetics/research-control/EVIDENCE.md)/[BACKLOG](kinetics/research-control/BACKLOG.md), Git-HEAD и CI |
| **Галактические гало, кривые вращения, профили** | [code/README](code/README.md), [начальная критика](00-critique-and-roadmap.md), [предсказания](01-new-numerical-predictions.md), результаты по задаче в [continuation-02](continuation-02/README.md) и последующих исследованиях |
| **Космологические наблюдательные ограничения / CLASS** | [cosmology/README.md](cosmology/README.md) → соответствующая проверка/результат в `cosmology/` |
| **Воспроизвести вычисление / проверить тесты** | [TESTING.md](TESTING.md), [kinetics/README.md](kinetics/README.md), версия кода и конкретный скрипт |
| **Полная научная история / отрицательные результаты** | [archive/README.md](archive/README.md) → [CHRONOLOGY](archive/CHRONOLOGY.md) → конкретный stage, commit SHA, JSON/CSV |
| **Выбор нового маршрута** | [PLAN.md](PLAN.md) → [подробный исследовательский план](plans/RESEARCH_ROADMAP_DETAIL.md) |

Не читать всю хронологию или 40-КБ STATE при любой мелкой правке. **Но** перед физическими выводами проверять соответствующие границы применимости, актуальные входные данные и версию кода. GitHub Actions разрешены для публичного репозитория, но временные artifacts не заменяют постоянный научный Git-архив.
