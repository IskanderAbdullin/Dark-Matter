# PLAN — карта исследований Dark-Matter

**Живой маршрут, не замена [GOAL.md](GOAL.md).** Научная цель охватывает космологическое происхождение бозонной тёмной материи и галактические распределения/гало. Нельзя считать всё исследование законченным из-за успеха отдельной модели.

| Если задача | Открыть |
|---|---|
| Определить новый научный этап, критерии и очередность | [Подробный изменяемый план](plans/RESEARCH_ROADMAP_DETAIL.md) |
| Перейти к кинетике, микрофизике и freeze-out | [Кинетическое ТЗ](kinetics/CONTINUOUS-RESEARCH-SPEC.md), [специальная цель](kinetics/research-control/GOAL.md) |
| Продолжить фактический расчёт | [STATE](kinetics/research-control/STATE.json) → [NEXT_ACTION](kinetics/research-control/NEXT_ACTION.json); открыть нужную часть [EVIDENCE](kinetics/research-control/EVIDENCE.md) |
| Проверить FRW/космологическую наблюдаемость | [cosmology/README](cosmology/README.md) |
| Исследовать гало, профили, кривые вращения | [PROJECT_MAP](PROJECT_MAP.md) → ветка галактических исследований |
| Проверить историю и воспроизводимость | [archive/CHRONOLOGY](archive/CHRONOLOGY.md), [archive/README](archive/README.md) |

**Правило:** PLAN/подробные научные планы можно менять при появлении доказательств, сохраняя причину, альтернативы, научные ворота и следы Git. GOAL и разрешения не менять без пользователя. Исследование не прекращается на первом CI или коммите: следующий релевантный шаг выбирать самостоятельно в пределах текущего запуска.

**Научный приоритет P0 после 2026-11-24:** минимальный real-φ Yukawa + одинаковый векторный U1 заряд χ1/χ2 сохраняет стабильного тяжёлого χ2 (20 МэВ). Его симметричный тепловой toy-реликт при y=.001 сильно перепроизводится даже без φ; см. [постоянный отрицательный физический аудит](kinetics/PHI-BOSE-STABLE-CHI2-THERMAL-RELIC-GATE-2026-11-24.md), [условный векторный SM portal](kinetics/PHI-BOSE-VECTOR-PORTAL-PARTIAL-SM-GATE-2026-11-24.md). **Проверять суммарную Ωχ2+Ωφ и реальные начальные условия в полной FD/BE+SM+FRW модели**, не игнорировать стабильного партнёра. Правило остановки только при объективно полной блокировке — в [AGENTS](AGENTS.md).

**Новый gate P0 (2026-11-24):** реальные CLEO R_had(9.996GeV) 3.497 ± статистическая/систематические ошибки измерены (но отличаются на 4MeV от mV10); точный Noether dark U1 ограничивает тяжёлую χ2 асимметрию вне MB; электронный A'10GeV νtransport/H≈1.56e−7 при χ2 химическом x8 тогда как точный парный Yukawa Γ/H≈4.97. **Нельзя** называть прежний toy Ωχ2pair h²≈880 настоящим relic указанного vector SM модели, поскольку общий Tχ2,φ,SM не выведен. Следующий шаг: отдельная dark/SM entropy+FRW, сумма Ωheavy+Ωφ и наблюдения. [Первичный исследовательский отчёт](kinetics/PHI-BOSE-CLEO-HADRONIC-U1-ASYMMETRY-AND-THERMAL-GATES-2026-11-24.md).

**P0 update (2026-11-24, 21:34 UTC):** [Variational entropy-energy/count gate](kinetics/PHI-BOSE-ENTROPY-ENERGY-GATE-2026-11-24.md) shows initial thermal + isolated benchmark requires at least ~0.89 neutrino energy units during BBN independently of later dark thermalization; scalar thermal freezeout Ω≈0.848; exact entropy-only today bound Ω≥0.0837 does not alone exclude. Prioritize actual initial entropy/heat export, lower-entropy nonthermal histories, and time-dependent BBN; count gate conditional on number nondecrease, not proven. No universal exclusion.
