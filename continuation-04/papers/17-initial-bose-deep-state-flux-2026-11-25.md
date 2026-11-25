# Initial Bose C22 redistribution into deeply bound states (2026-11-25)

Companion to paper 16. Use dimensionless center-of-mass C=(q1+q2)/2, d=(q1-q2)/2, A=C²+d², C=sqrt(A)cosθ, d=sqrt(A)sinθ, B=A sin2θ. Incoming squared speeds A±Bu; outgoing A±Bv. Angular averaging and relative speed yield weak-form measure 8 A^(5/2) cos²θ sin³θ dA dθ du dv on 0≤u≤v≤1 after folding angular signs and pairing forward/reverse events. With P_u=f1 f2, S_u=f1+f2, the reversible factor is P_u(1+S_v)−P_v(1+S_u); quartic terms cancel algebraically. Exact two-body particle number and kinetic energy conservation follows independently of quadrature. For cumulative number below q²=P, only u<t<v, t=|A-P|/B, contributes; count change sign(A-P). Initial support requires A≤ψ and u≤min(1,(ψ-A)/B). Mapping θ to [.5 asin(min(1,|A-P|/A)), π/4] with reflected weight eliminates zero kinematic region. Split outgoing integral at support boundary. Same dimensional prefactor as paper 16.

At bound linear lowered BE, L=ln(1+1/f) as a function of q² has strictly positive second derivative: for b=1/expm1(η−ψ+q²), c=1/expm1(η), 
L''=b c(b+1)[2b²−2bc+2b−c+1]/[(b−c)²(b−c+1)²]>0.
Thus greater energy spread increases L1+L2, giving positive reversible factor toward the wider pair; 4th moment of speed increases initially (pair increment 2B²(v²-u²)). Above-cutoff outgoing pairs have zero reverse rate. This does not mean every cumulative threshold flux has the same sign.

Define deepest particles by η−s<d0. Critical halo center ψ(0)=η. The fraction of all halo particles in this region asymptotically scales as (6^(3/2)π/(64 Ic Qt))*d0², obtained from ψ=η−r²/6+o(r²) and f~1/(r²/6+q²), integrating r²dr q²dq. Coefficients .0251749234/.0272427331 for the two models. Same physical parameters as paper 16.

| d0 | linear fraction | linear initial rate per deep population /Gyr (mapped n24,nr40) | quadratic fraction | quadratic rate /Gyr (mapped n24,nr40) |
|---:|---:|---:|---:|---:|
| .2 | .001383669 | .11788635 | .001534494 | .16690623 |
| .1 | .0003172621 | .13943062 | .0003495787 | .22158450 |
| .03 | 2.582282e-5 | .16354118 | 2.823451e-5 | .29238399 |
| .01 | 2.718300e-6 | .17491193 | 2.959625e-6 | .33013757 |
| .003 | 2.364023e-7 | .18081212 | 2.566914e-7 | .35136832 |

Negative numerical finding: unmapped n24 yields linear .01373134/Gyr instead of .18081 at d0=.003. Mapped n24 compared with independent log-mapped n16/n24 (same physical integrals, different quadrature) agrees within ~0.04% across linear rows. Mapped n32 independently reproduces paper 16 cutoff at ψ=1 and 3.9 to <1.2e-6 relative; mapped n24 quadratic agrees to 7.1e-6. Finite-radius integrated very deep flux tiny, linear 4.2744e-8 of all particles/Gyr at d0=.003. Do not equate fractional deep-population growth to macroscopic transfer.

## Independent soft-mode coefficient
At critical center, q1→0 and f1→∞; for fixed other squared speed S, q3²=zS, q4²=(1−z)S, z uniform. C[f]/f1 tends pointwise to Dσsqrt(2)V0 times integral (S/2)dS dz [f3 f4−f2(1+f3+f4)] over S≤2η. For S≤η and linear lowering, the exact BE identity b3 b4=b2(1+b3+b4) simplifies this bracket to c(1+2b2)−c², independent of z. For η<S<2η, f2=0 and inverse gain from two occupied modes must be retained; integrate z from 1−η/S to η/S. For quadratic lowering evaluate the full bracket. Values of dimensionless integrals over S≤η / S>η: linear .1022043231/.0185961092, quadratic .2161453437/.00601752835. The corresponding fractional pointwise soft-limit coefficients are .18491053 and .36747410/Gyr, consistent with the independently computed deep-population initial rates as d0 decreases. This establishes a finite positive diagnostic, not a rigorous theorem of uniform IR limit or condensate nucleation; coupled transport, changing potential, condensate feedback and correlations remain uncomputed.

Source continuation-04/code/lowered_bose_weak_flux.py; local SciPy/NumPy and symbolic derivative checks. Remote CI not verified. No new fitted observational conclusion.
