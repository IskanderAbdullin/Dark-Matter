# 21. Проверка известной теоремы о конечном радиусе и массе: случай Bose F_k

21 ноября 2026. Это теоретический аудит предыдущих результатов. Общий порог Lane–Emden n=5 НЕ заявляется новым.

Для классического изотропного политропного распределения f(E)=(E_t−E)_+^k при −1/2<k<7/2 сферические решения имеют конечные массу и радиус; k=7/2 отвечает Lane–Emden n=5 с бесконечным радиусом и конечной массой. При k>7/2 классическое политропное решение не имеет одновременно конечной массы и радиуса. Более общая теорема Rein–Rendall (2000) даёт условия на фазовое распределение, локальные у энергетической границы, обеспечивающие конечную массу и носитель; применять её к другой формуле f(E) можно лишь после проверки точных гипотез.

Положительный Bose Taylor-remainder F_k из заметки 15 удовлетворяет F_k(E)~C_k(E_t-E)^k при E→E_t^−, C_k>0, и поэтому даёт внешний показатель плотности n=k+3/2. Но полный профиль НЕ является политропой: внутри он бозонный, а строго критический центр обладает неограниченным occupation. Поэтому n=5 threshold — мощная аналогия и математическая гипотеза для Bose, а не автоматически готовая классификация для всех k и глубин.

Численный статус: k1,k2 дают конечный радиус; k3 может требовать ~10^9 R для первого нуля; k4–6 на испытанных центрах не пересекли ноль до 10^10 R, а их масса продолжила расти согласно внешней степенной асимптотике. Последнее — результат ограниченного численного интервала, не всеобщая теорема.

**Не путать:** конечность носителя, collisionless linear stability (Antonov/DFB/Pérez–Aly) и полное quantum Bose collision equilibrium — три разные свойства.

Ключевые источники:
- Rein & Rendall, Compact support of spherically symmetric equilibria in non-relativistic and relativistic galactic dynamics, Math. Proc. Cambridge Philos. Soc. 128 (2000), https://doi.org/10.1017/S0305004199004193.
- Guo & Rein, Isotropic steady states in galactic dynamics revised, https://arxiv.org/abs/math-ph/0010031.
- Isotropic steady states in galactic dynamics, описание известного политропного диапазона, https://citeseerx.ist.psu.edu/document?doi=18beb63d2a558712573a8e7e0db71e2305545819&repid=rep1&type=pdf.
- Spherically symmetric equilibria for self-gravitating kinetic or fluid models — A simple proof for finite extension, https://arxiv.org/abs/1203.2590.

Следующий математический шаг: буквально сверить Rein–Rendall hypotheses для subcritical Bose F_k и строго доказать или ограничить радиальную compactness классификацию; особый случай точно критического центра fmax=∞ требует отдельного функционального анализа.
