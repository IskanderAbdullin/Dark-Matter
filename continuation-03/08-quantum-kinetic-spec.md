# Минимальная кинетическая постановка: от холодного реликта к гало с near-critical Bose occupation

## 1. Почему стационарный интеграл по импульсам недостаточен
Исходная работа JCAP 2021 вычисляет плотность как интеграл от \(f_{\rm BE}(\mathbf p,r)\). Она не выводит \(f_{\rm BE}\) из evolution of initial cosmological \(f(q)\). Разделим отдельно:
- collisionless gravitational transport of \(f(\mathbf x,\mathbf p,t)\);
- quantum \(2\leftrightarrow2\) scattering;
- possible inelastic \(1\leftrightarrow2\), \(2\leftrightarrow3\) reactions and condensate transfer;
- separate Bose-star population phase-space DF \(f_{\rm BS}(\mathbf x,\mathbf v,M,t)\).

## 2. Минимальный quantum kinetic Boltzmann–Poisson
\[
\partial_t f+\frac{\mathbf p}{m}\cdot\nabla_x f
-m\nabla\Phi\cdot\nabla_p f=C_{22}[f]+C_{12}[f,n_c]+C_{\rm inelastic}[f].
\]
\[
\nabla^2\Phi=4\pi G(mn_{\rm gas}+\rho_c+\rho_{\rm BS}+\rho_{\rm baryon}),\quad
n_{\rm gas}=\int\frac{d^3p}{(2\pi\hbar)^3} f.
\]
Здесь \(f\) следует понимать как **occupation number**, в отличие от source JCAP, где функция \(f_{\rm source}\) уже включает фактор \((2\pi\hbar)^{-3}\). При численной реализации необходимо зафиксировать эту нормировку без двойного учета quantum cell.

Схематическая collision integral \(2\leftrightarrow2\) (без нормировочных констант матричного элемента):
\[
C_{22}[f_1]\propto\int d\Pi_{234}
(2\pi)^4\delta^{(4)}(p_1+p_2-p_3-p_4)|\mathcal M|^2
[(1+f_1)(1+f_2)f_3f_4-f_1f_2(1+f_3)(1+f_4)] .
\]
Знак зависит от выбранной convention вход/выход; важны Bose enhancement \((1+f)\), энергия, импульс и conservation of particle number в чистом \(2\leftrightarrow2\). Для BEC self-consistency необходим \(C_{12}\) с exchange of particles between zero-mode and thermal excitations. Если реальное поле имеет number-changing channels, \(\mu\) не может объявляться произвольным conserved chemical potential во все эпохи.

Проверка стационарности: \(f=[\exp((\epsilon-\mu)/kT)-1]^{-1}\) обращает \(C_{22}=0\) вследствие detailed balance, если кинетическая энергия и химпотенциал совместимы с инвариантами. Но **равенство collision integral нулю не доказывает**, что соответствующий initial condition достигнет этого состояния за время галактики.

## 3. Реальные времена и выбор режима
- Параметр \(\lambda_{\rm th}\) на 86 eV галактическом гало — десятки микрометров, и центральное \(n\lambda^3\sim2.6\). Отсюда **не** следует, что весь газ имеет \(f\gg1\).
- Для выбранной upper-limit-scale constant \(\sigma/m=1.25\) cm²/g \(\ell_0/R\approx64\) и single-pass optical depth до 10R лишь ~0.033 (без stimulation). Это неудобный Knudsen regime для обычной локальной hydrodynamic closure.
- В идеальном \(z\to1\) низкоимпульсное \(f\propto1/p^2\) требует не считать \((1+f)\) малыми при всех p. **Но доля по числу частиц с mode occupancy >10 только ~0.263**, >100 ~0.086; оценка одного multiplication factor неверна.
- Гравитационный Bose-star kinetic mechanism Levkov et al. 2018 опирается на режимы, подробно разработанные для fuzzy/wave DM. При 86 eV использование его \(t_{\rm nuc}\) требует явного checking occupancy, coherence, granularity and interaction kernels.
- Быстрая violent relaxation может поменять coarse-grained \(f(E)\), но не заменяет доказательства появления квантового Bose equilibrium.

## 4. Набор обязательных математических тестов
1. Collision conservation: \(\int C_{22}\,d^3p=0\), \(\int p_i C_{22}\,d^3p=0\), \(\int p^2 C_{22}/(2m)\,d^3p=0\) для закрытых упругих реакций в однородном поле.
2. Positivity of \(f\), correctness of stimulated factors and detailed balance, entropy production (quantum H-theorem) when applicable.
3. Time-dependent \(n_c\) and \(n_{\rm gas}\) including charge/number constraints; no hidden loss of particle mass to a Bose-star population.
4. Vlasov limit \(C\to0\) must retain phase-space transport and possible stationary \(f(E)\) without asserting collisionally established \(\mu\).
5. Kinetic collisionless validation and asymptotic \(x\to0\): recovery of original Bose–Poisson potential equation, if the initial data are **already** its exact stationary \(f(E)\).
6. Bose-star nucleation/evaporation in a bath: use threshold rates and total energy conservation; literature explicitly includes both signs of mass transfer (Chan et al. 2024).

## 5. Concrete fail/pass
- PASS for thermal interpretation only if \(\tau_{\rm relax}(r)\ll t_{\rm halo}\) across **the region whose equilibrium Bose distribution is used** and solution remains consistent with cosmological initial \(f(q)\).
- If only center thermalizes while outer halo collisionless, use a radially matched kinetic/phase-space model rather than one isothermal chemical potential to infinity.
- If no local equilibration occurs but an energy-dependent \(f(E)\) is dynamically admissible, label result **kinematic phenomenological halo**, not a thermally equilibrated Bose DM gas.
- Compact BEC fraction: demand independent \(t_{\rm nuc}\), \(dn_{\rm BS}/dM\), spatial distribution and observational test before attributing dwarf-satellite scarcity to Bose stars.

## 6. Recommended foundations
- Levkov, Panin & Tkachev (2018), https://doi.org/10.1103/PhysRevLett.121.151301
- Bar-Or, Fouvry & Tremaine (2021), https://doi.org/10.3847/1538-4357/abfb66
- Chan et al. (2024), https://doi.org/10.1007/JHEP01(2024)071
- Kinetic relaxation with interactions (2024), https://doi.org/10.1103/PhysRevD.109.016002
- Finite-temperature HFB (2022), https://doi.org/10.1103/PhysRevD.106.045025
- Stimulated production Sakurai & Yin (2025), https://doi.org/10.1007/JHEP03(2025)202

**Статус:** математическая спецификация и необходимые проверки, **не реализованная** нелинейная kinetic solver и не доказательство реального срока.
