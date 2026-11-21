# 13. Exact integrability threshold for the critically occupied six-dimensional Bose halo

**Mathematical consequence of the phase-volume calculations; 21 November 2026.** Not an observational bound or proof that a critical distribution forms from cosmological conditions.

Consider the *exactly centrally saturated* Bose distribution (either untruncated locally or finite-mass lowered Bose with an outer escape cutoff). At the center \(z_0=1\), \(y(r)\sim(r/R)^2/6\), and the dimensionless occupation behaves
\[
f(\mathbf x,\mathbf p)\sim [q^2+\xi^2/6]^{-1},\qquad q=p/(\sqrt2mV_0),\quad \xi=r/R.
\]
This is a **six-dimensional** singularity at the isolated phase-space point \((\mathbf x,\mathbf p)=(0,0)\). For its superlevel volume \(V(F)=\int d^3x d^3p\,\mathbf 1_{f>F}\), paper 03 derived
\[
\boxed{V(F)\sim 4\sqrt3\,\pi^3R^3(mV_0)^3 F^{-3}}\quad(F\to\infty).
\]
The positive *lowered* Bose cutoff shifts \(F\to F+c_\eta\) but does not change this leading exponent. The central phase-space cell factor \((2\pi\hbar)^{-3}\) may be included consistently and does not change integrability.

Use the layer-cake identity for any \(p>0\):
\[
\int f^p\,d\Gamma=p\int_0^\infty F^{p-1}V(F)\,dF.
\]
At high occupancy the integrand is proportional to \(F^{p-4}\), which is integrable at infinity exactly when \(p<3\). For a **finite-mass, finite-support lowered Bose halo**, the low-occupation contribution is finite. Therefore:
\[
\boxed{f\in L^p(d^3x\,d^3p)\ \text{for}\ 0<p<3,\qquad f\notin L^p\ \text{for}\ p\ge3.}
\]
At \(p=3\) the divergence is logarithmic; for \(p>3\) it is a power law. Pointwise \(f\notin L^\infty\). This distinction is important when applying global nonlinear Vlasov theorems requiring bounded distribution functions or weighted \(L^p\) control.

**Local vs global collision enhancement:** At precisely \(r=0\), \(\int d^3p\,f^2\) diverges in an ideal critical Bose distribution because \(f\sim p^{-2}\). But integrating the same \(f^2\) over a three-dimensional neighborhood of the **single** critical center is finite: the local six-dimensional measure behaves \(s^5ds\), while \(f^2\sim s^{-4}\), yielding \(s\,ds\). The third phase-space moment behaves \(s^{-1}ds\) and diverges logarithmically. Thus the previous divergence of a local occupation-weighted susceptibility at the central point does **not** imply that every globally integrated quadratic occupation measure is infinite. Nor does finite \(\int f^2\) establish a finite actual Bose collision rate: energy-momentum constraints, gain-loss cancellations and interaction kernels matter.

**Initial relic comparison:** A homogeneous, zero-\(\mu\), relativistically frozen Bose relic has \(f_i(q)\sim p_0/q\) at small q within a finite spatial patch. Its high-\(F\) phase-volume tail also scales \(F^{-3}\), hence the same \(L^p\) threshold 3. This does **not** guarantee a collisionless match between the complete distributions; papers 03 and 10 compute the subleading six-dimensional tail and the Dehnen excess-mass inequalities.

**Physical qualifications:** The ideal critical point has infinite f only in a vanishing phase-space set, and is expected to be modified by finite volume, interactions, coherent condensation, quantum fluctuations or infrared regularization. These mathematical \(L^p\) results apply to the **specified classical occupation ansatz**, not automatically to a quantum Wigner function of a real interacting field. For the original untruncated halo, mass diverges in the far outer tail, so the global \(L^p\) statement above explicitly assumes a *finite outer cutoff*; the near-center divergence statement does not.

Sources: Abdullin & Popov, JCAP 2021 055 https://doi.org/10.1088/1475-7516/2021/11/055; Tremaine–Hénon–Lynden-Bell mixing https://doi.org/10.1093/mnras/219.2.285; Dehnen phase-space excess mass https://doi.org/10.1111/j.1365-2966.2005.09099.x; our papers 03,10.
