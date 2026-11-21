# Устойчивость статического бозонного гало: кинематическая, термодинамическая, динамическая

## Почему необходимы три несмешиваемых уровня
(a) **Положительность** phase-space DF \(f\ge0\) и \(M<\infty\) — возможность описания некоторой популяции частиц.  
(b) **Стационарность** под Vlasov: \(f=F(E,L)\) при неподвижном сферическом потенциале; справедлива без столкновений.  
(c) **Устойчивость** к малым/большим perturbations, и (отдельно!) устойчивость к столкновительной гравотермальной эволюции — не следуют автоматически из (a),(b).

## Конкретный lowered Bose DF из continuation-02
Для \(\epsilon<E_t\):
\[
f(\epsilon)=\frac{1}{\exp[(\epsilon-\mu)/kT]-1}
-\frac{1}{\exp[(E_t-\mu)/kT]-1},\qquad
f=0,\;\epsilon\ge E_t .
\]
При \(\mu<E_{\min}\) знаменатели положительны; \(f\) неотрицательна. Для \(\epsilon<E_t\):
\[
\frac{df}{d\epsilon}=
-\frac{\exp[(\epsilon-\mu)/kT]}
{kT[\exp((\epsilon-\mu)/kT)-1]^2}<0.
\]
Функция непрерывна на \(E_t\) (обе стороны равны нулю); у неё конечная левая производная, правая 0, поэтому derivative имеет скачок, что требует аккуратности при проверке теорем/численных аппроксимаций. При \(\mu\to E_{\min}\) в центре возникает квантово критический IR-предел, который не обязательно лежит в классе регулярности, требуемом для используемых stability theorems.

Гладкая ветвь распределений с \(\psi_0<\eta\) уже дает 4 положительных конечномассовых решения ODE; новая прямая проверка \(df/d\epsilon<0\) создает перспективную связь с математической устойчивостью isotropic Vlasov systems. Тем не менее **отрицательная производная — не готовое доказательство всех форм нелинейной устойчивости**, нужно сверить условия compact support, class of perturbations и energy–Casimir.

## Конкретная prior art
1. Guo & Rein, *A Non-Variational Approach to Nonlinear Stability in Stellar Dynamics Applied to the King Model*, Communications in Mathematical Physics 271 (2007) 489–509, DOI https://doi.org/10.1007/s00220-007-0212-8. В статье доказана нелинейная устойчивость *King* против специального класса сферически-симметричных dynamically accessible perturbations; это НЕ автоматически доказательство для lowered Bose.
2. Guo & Rein, *Stable steady states in stellar dynamics*, arXiv:math-ph/9806006; energy–Casimir подход и устойчивость класса stationary states.
3. Chavanis, *Dynamical stability of collisionless stellar systems and barotropic stars: the nonlinear Antonov first law* (2006), https://arxiv.org/abs/astro-ph/0602429. Различает dynamical Vlasov stability от устойчивости соответствующей barotropic Euler fluid.
4. Kirejczyk, Müller & Chavanis, PhysRevE 105 034145 (2022), https://doi.org/10.1103/PhysRevE.105.034145, показывает сложную фазовую структуру самогравитирующего идеального Bose-газа, multiple macrostates и необходимость сравнения free energies и caloric curves.
5. Slepian & Goodman, *Ruling out bosonic repulsive dark matter in thermal equilibrium* (2012), https://doi.org/10.1111/j.1365-2966.2012.21901.x. Его отрицательный результат относится к определённой thermally equilibrated repulsive-core model, а не безоговорочно к любой bosonic halo. Но совместимость теплового равновесия и контактного сечения надо проверять прямо.
6. Связанная современная HF/HFB физика: *Finite temperature description of an interacting Bose gas*, PhysRevD 106 045025 (2022), https://doi.org/10.1103/PhysRevD.106.045025. Внутри condensed phase простое добавление \(2gn\) может быть недостаточно; нужны корреляции, корректная химическая связь компонент.

## Где сильная задача для новой теоремы
- Возьмем \(f(E)=A[n_{\rm BE}(E)-n_{\rm BE}(E_t)]_+\), \(\mu<E_{\min}\), с самосогласованным \(E=m\Phi+p^2/(2m)\) и плотностью из импульсного интеграла.
- Сформулируем корректно steady Vlasov-Poisson и компактность носителя, контроль зависимости \((\eta,\psi_0)\), массу \(M\), kinetic K, potential W и boundary conditions.
- Проверим потребные regularity / monotonicity conditions известных energy–Casimir теорем и построим соответствующую функциональную \(Q(f)\) на допустимой области.
- Проведем radial spectral perturbation test и затем rigorously bounded stability либо найдем instability ranges.
- **Отдельно** thermodynamics: конечное самогравитирующее Bose-облако требует ансамбля и полного T, entropy/free-energy treatment; отрицательная теплоемкость/gravothermal catastrophe не эквивалентна collisionless instability.
- **Отдельно** quantum critical: избегать подмены \(\mu\to E_{\min}\) корректной interacting HF/HFB critical line.

## Отдельная проблема числа конденсатов
В JCAP-модели Bose stars описаны собственной изотермической Jeans density \(x_{bs}\propto \exp[-\gamma y]\) с параметрами \(\alpha,\gamma\). Такая smooth density может представлять только **усреднение большого числа объектов**. Чтобы это действительно было моделью, нужно найти их распределение по массам, phase-space DF \(f_{bs}(E,L,M)\), velocity dispersion, tidal survival. Принятие \(\alpha=\gamma\) — дополнительное приближение source, не из Hamiltonian или равновесия.

## Доказательность
- \(df/dE<0\) и \(f\ge0\): **вывод непосредственно из выбранного анзаца**.
- Кинетическая стационарность: **следует из зависимости \(f=f(E)\) в фиксированном сферическом потенциале**, если Poisson самосогласовано выполнено.
- Нелинейная Vlasov stability: **не установлена**.
- Термодинамическая stability на космологических временах: **не установлена**.
