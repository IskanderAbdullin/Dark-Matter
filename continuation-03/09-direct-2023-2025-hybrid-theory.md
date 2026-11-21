# Самый близкий к исходной гипотезе готовый математический фундамент: Proukakis–Rigopoulos–Soto (2023–2025)

**Статус проверки:** издательские аннотации/описания и открытая arXiv html версия статьи II изучены; это направленный теоретический cross-paper разбор, не независимо воспроизведенная Schwinger–Keldysh функциональная интеграция.

## Почему эта линия более релевантна 86 eV, чем копирование soliton/fuzzy DM
JCAP 2021 описывает доминирующие *неконденсированные* Bose частицы вместе с compact Bose stars. Главная слабость — отсутствие динамических уравнений образования и обмена компонент. Proukakis, Rigopoulos & Soto (PRD108.083513, 2023) независимо строят из единого nonrelativistic quartic bosonic action Gross–Pitaevskii для coherent component, Boltzmann-like kinetic equation для fast modes, Poisson с обеими плотностями, collision terms. PRD110.023504 (2024) изучает coupled linear cosmological fluid perturbations. PRD111.023505 (2025) добавляет stochastic forces и dissipation, через fluctuation-dissipation relation, включая nucleation from initially zero condensate coherent mean field through fluctuations.

Эти статьи **не являются** прямым продолжением работы Abdullin & Popov и не показывают, что при \(m=86\) eV реализовано нужное равновесие. Но они почти точно заполняют структуру отсутствующего dynamical sector и дают организованный порядок получения модели из первых принципов.

## Три уровня переменных и общие физические сохраняемые величины
1. Noncondensed Wigner distribution \(f(\mathbf x,\mathbf p,t)\) с нормировкой \(\int d^3p\,f/(2\pi\hbar)^3\).
2. Coherent condensate field \(\Psi(\mathbf x,t)\), \(n_c=|\Psi|^2\).
3. Discrete/mesoscopic Bose stars: функционал распределения \(F_{\rm BS}(\mathbf x,\mathbf v,M,t)\). **Этот третий уровень не следует автоматически** из двухполосного разложения coherent/incoherent: требуется модель selfbound fragmentation, condensation, mergers, tidal evolution.

Переменные \((\Psi,f)\) подвергаются Poisson gravity вместе с baryons; в простом number-conserving quartic model
\[
\frac{d}{dt}\left[\int n_c d^3x + \int n_n d^3x\right]=0.
\]
При добавлении индивидуальных Bose stars баланс должен включать \(N_{\rm BS}\), transfer terms и escaping particles. Для точного real scalar с number-changing interactions particle number не является фундаментальным U(1) charge. Нельзя требовать одновременно долговечного строго сохраненного \(\mu\) и физических inelastic channels, которые стирают particle number.

## Перевод на исходные JCAP параметры
- В source JCAP \(V_0=\sqrt{kT/m}\), \(\lambda_{\rm th}=\sqrt{2\pi}\hbar/(mV_0)\), \(\nu=n_0\lambda^3\), \(R^2=V_0^2/(4\pi Gmn_0)\), \(\alpha\) condensed-star smooth-density ratio и \(\gamma=V_0^2/v_{\rm BS}^2\).
- В kinetic theory \(\nu\) есть **локальная начальная или эволюционирующая phase-space degeneracy**, а не безусловно константа в пространстве. Если достигается local BE, \(\nu(r)=Li_{3/2}(z(r))\), \(z(r)\) следует из \(\mu,T,\Phi\) и mean field, причем source \(\nu\) обозначает центральный \(\nu(0)\).
- Dynamic boson star population \(\alpha,\gamma\) должны быть вычислены из \(F_{\rm BS}\) и transport, а не принудительно \(\alpha=\gamma\). Простая \(\Psi\)-полоса может хранить много coherent modes и не должна автоматически интерпретироваться как один сферический compact star.
- В идеальном критическом центре source \(z\to1\), \(f(p)\sim p^{-2}\), интегрируемый в 3D, но Bose-enhanced scattering and stochastic noise могут стать важны; необходим IR regulator and consistent self-interaction.

## Сравнение работы 2024 с thermal relic free streaming
PRD110.023504 study uses a *hydrodynamic limit of coupled fluids* with assumption on equation of state and reports parameter ranges with localized enhancement in \(P(k)\). Это не значит, что взаимодействующий конденсат создается в достаточно раннюю эпоху при \(m=86\) eV или компенсирует подавление \(\sim\) Mpc scales. Сначала derive production fractions, condensate history, m,g and assess validity of hydro closure (central mean-free-path benchmark \(\ell/R\sim64\) at F583-4 illustrative maximum cross-section suggests local fluid closure may be problematic in halo environment, though early Universe has different density). Следовательно нельзя без проверки использовать их графики вместо CLASS transfer for actual \(f(q)\).

## Что в 2025 очень перспективно, но не доказано для нас
Статья II указывает, что stochastic coherent fluctuations могут индуцировать конденсат даже при \(\Psi(t=0)=0\), и выводит fluctuation-dissipation form in appropriate thermal limit. Это гораздо сильнее статического предположения «звезды уже есть». **Однако** темп нуклеации, amplitude of noise, correlation length and wave occupation in m=86eV halo не известны. Quantum state (coherent vs incoherent vs number) должен быть задан физически.

## Постановка минимального проекта на основе их уравнений
A. Строго выписать лагранжиан complex scalar, covariant or nonrelativistic limit и сохраненный U(1), зафиксировать source normalization \(g\) and scattering cross-section.
B. Отдельно реализовать two-component (GP+kinetic) with number-conserving \(C_{12},C_{22}\) and conserved mass/energy before introducing stochastic terms.
C. Проверить при \(g=0\), \(n_c=0\) точное восстановление stationary source BE–Poisson, только если это заданный initial distribution; при nonzero \(g\) проверить HF boundary layer.
D. Включить stochastic nucleation based on PRD111.023505, проверяя positivity, noise spectra, FDT, N/BEC fraction, and finite-size Bose-star fragmentation.
E. Параллельно ранняя Вселенная: production Yin (2023) vs Sakurai–Yin (2025), Boltzmann transfer vs constraints; единую пару \((m,g)\) и symmetry/production couplings использовать в обеих эпохах.
F. Отслеживать observable quantities: \(P(k,z)\), halo Vc(r), \(\sigma_{\rm DM,1D}(r)\), kurtosis, \(\rho_{\rm BS}(r)\), Bose-star function \(dN_{\rm BS}/dM\), mass loss, lensing.

## Итог для приоритетов
Не начинать с пятой альтернативной гидродинамической closure и не переписывать готовую статическую BE модель просто ради еще одной кривой вращения. **Первоочередной теоретический проект:** честное перенесение общей GP–kinetic–Poisson схемы Proukakis/Rigopoulos/Soto на конкретные \(m\sim86\) eV и coupling, с проверкой режима validity и отдельной фазовой космологией холодного производства.

**Первоисточники:** https://doi.org/10.1103/PhysRevD.108.083513 ; https://doi.org/10.1103/PhysRevD.110.023504 ; https://doi.org/10.1103/PhysRevD.111.023505 ; https://arxiv.org/abs/2407.08860 ; источник модели https://doi.org/10.1088/1475-7516/2021/11/055 .
