"""Independent analytic angular spin-trace control for finite vector exchange.

Only the vacuum spin-trace angular moment is integrated here. The full
finite-temperature five-dimensional Pauli heat cannot factor this moment:
boosted final-state Fermi factors and energy transfer depend on both angles.
"""
from __future__ import annotations
import math


def moments(b: float) -> tuple[float, float, float]:
    """J_n(b)=integral_0^2 y**n/(1+b*y)**2 dy, n=0,1,2.

    b=2*pstar**2/mV**2. For small b use the uniformly convergent
    geometric derivative series, avoiding catastrophic cancellation.
    """
    if not math.isfinite(b) or b < 0:
        raise ValueError("b must be nonnegative and finite")
    if b < 0.02:
        out = []
        for n in range(3):
            total = 0.
            for k in range(100):
                term = (-1.)**k*(k+1)*b**k*2.**(n+k+1)/(n+k+1)
                total += term
                if abs(term) < 2e-16*max(1., abs(total)):
                    break
            out.append(total)
        return tuple(out)
    z = 2*b
    log = math.log1p(z)
    return (2/(1+z),
            (log+1/(1+z)-1)/b**2,
            (z-2*log-1/(1+z)+1)/b**3)


def trace_angular(s: float, M: float, m: float, mV: float) -> float:
    """Integral_{-1}^{1} dc [D²+U²+2(M²+m²)t] F(t).

    D=s-M²-m², U=-D-t, t=-2*pstar²*(1-c);
    F(t)=[mV²/(mV²-t)]². Excludes overall 2G² and Pauli.
    """
    if min(s, mV) <= 0 or M < 0 or m < 0 or s < (M+m)**2:
        raise ValueError("invalid on-shell inputs")
    ps2 = max(0., (s-(M+m)**2)*(s-(M-m)**2)/(4*s))
    a = 2*ps2
    D = s-M*M-m*m
    B = M*M+m*m
    j0,j1,j2 = moments(a/mV**2)
    return 2*D*D*j0-2*a*(D+B)*j1+a*a*j2
