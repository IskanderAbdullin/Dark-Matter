"""Independent five-dimensional CM validation of the massive chi-lepton q kernel.

Conditional point-vector four-Dirac EFT only; no photon/nu plasma, no FRW,
no kinetic eigenvalue, no Bose relic. Energies GeV, G GeV^-2, Q GeV^5.
CM calculation does NOT call the transfer kernel. Its five integration
coordinates are incoming |p|,|k|, their angle, CM scattering polar angle,
and scattering azimuth. The second calculation uses the separately coded
invariant ell-cone kernel from finite_mass_chie_kernel.py.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit

from kinetics.finite_mass_chie_kernel import kernel


def gauss(a: float, b: float, n: int):
    x, w = leggauss(n)
    return a + (b-a)*(x+1)/2, w*(b-a)/2


def cm_heat(M, m, G, Te, Tchi, *, pmax, kmax, nr=24, nang=24,
            nphi=24, gchi=4., ge=4., mediator_mass=None):
    """Independent relativistic full 5D thermal weak Q_chi in GeV^5.

    The 1/2 multiplies the forward-minus-reverse integrand because the
    integration covers the ENTIRE initial/final domain, not one orientation.
    This independent code uses CM Lorentz boosts, not transfer cones.
    """
    if mediator_mass is not None and (not math.isfinite(mediator_mass) or mediator_mass<=0):
        raise ValueError('mediator_mass must be finite and positive')
    p, wp = gauss(0., pmax, nr)
    k, wk = gauss(0., kmax, nr)
    mu, wm = gauss(-1., 1., nang)
    c, wc = gauss(-1., 1., nang)
    phi, wf = gauss(0., 2*math.pi, nphi)
    k, wk = k[:,None,None,None], wk[:,None,None,None]
    mu, wm = mu[None,:,None,None], wm[None,:,None,None]
    c, wc = c[None,None,:,None], wc[None,None,:,None]
    phi, wf = phi[None,None,None,:], wf[None,None,None,:]
    w = np.sqrt(k*k+m*m)
    accum = 0.
    for pi, wi in zip(p, wp):
        E = math.hypot(pi, M)
        Et = E+w
        P2 = pi*pi+k*k+2*pi*k*mu
        P = np.sqrt(np.maximum(P2, 0.))
        s = Et*Et-P2
        root = np.sqrt(s)
        ps2 = np.maximum(0., (s-(M+m)**2)*(s-(M-m)**2)/(4*s))
        ps = np.sqrt(ps2)
        beta = P/Et
        gamma = Et/root
        parallel = np.where(P>0, pi*(pi+k*mu)/np.maximum(P,1e-100),0.)
        cb = np.clip(gamma*(parallel-beta*E)/np.maximum(ps,1e-100), -1., 1.)
        cbp = np.clip(cb*c+np.sqrt(np.maximum(0.,1-cb**2))*
                      np.sqrt(np.maximum(0.,1-c**2))*np.cos(phi),-1.,1.)
        Estar = (s+M*M-m*m)/(2*root)
        Ep = gamma*(Estar+beta*ps*cbp)
        wprime = Et-Ep
        t = -2*ps2*(1-c)
        D = s-M*M-m*m
        U = -D-t
        amplitude = 2*G*G*(D*D+U*U+2*(M*M+m*m)*t)
        if mediator_mass is not None:
            # G is the zero-transfer contact coefficient gchi*glep/mV^2.
            # No width is needed for spacelike t<=0; this is tree-level only.
            amplitude *= (mediator_mass**2/(mediator_mass**2-t))**2
        fi = expit(-E/Tchi)
        gi = expit(-w/Te)
        fo = expit(-Ep/Tchi)
        go = expit(-wprime/Te)
        forward = fi*gi*(1-fo)*(1-go)
        reverse = fo*go*(1-fi)*(1-gi)
        dot = E*w-pi*k*mu
        vm = np.sqrt(np.maximum(0.,dot*dot-M*M*m*m))/(E*w)
        terms = (wi*pi*pi*wk*k*k*wm*wc*(wf/(2*math.pi))*
                 (vm/s)*amplitude*(Ep-E)*(forward-reverse))
        accum += float(np.sum(terms))
    return accum*gchi*ge/(512*math.pi**5)


def transfer_heat(M, m, G, Te, Tchi, *, xmax=30., n=64, nq=24,
                  gchi=4., ge=4.):
    """Independent 3D energy collocation over gain-only E'-E>0 orientation.

    Uses the proposed ell-cone K, with no bin-interpolation or rate averaging.
    """
    E, wE = gauss(M, M+xmax*Te, n)
    omega, ww = gauss(m, m+xmax*Te, n)
    z, wz = gauss(0., 1., n)
    result = 0.
    for i, Ei in enumerate(E):
        fi = expit(-Ei/Tchi)
        for j, wj in enumerate(omega):
            gi = expit(-wj/Te)
            delta = (wj-m)*z
            Ep = Ei+delta
            wp = wj-delta
            fo = expit(-Ep/Tchi)
            go = expit(-wp/Te)
            forward = fi*gi*(1-fo)*(1-go)
            reverse = fo*go*(1-fi)*(1-gi)
            K = kernel(M,m,G,float(Ei),float(wj),Ep,wp,nq=nq,
                       gchi=gchi,glep=ge)
            result += wE[i]*ww[j]*(wj-m)*float(
                np.dot(wz,K*delta*(forward-reverse)))
    return result


def benchmark(*, temperature_ratio=1., chi_temperature_ratio=.9,
              nr=(24,32,40), ne=(64,96,160), nq=24, xmax=30.,
              electron_mass=.00051099895, G=1.):
    m = electron_mass
    M = .1*m
    T = temperature_ratio*m
    Tc = chi_temperature_ratio*T
    pmax = math.sqrt((M+xmax*T)**2-M*M)
    kmax = math.sqrt((m+xmax*T)**2-m*m)
    reference = {str(n):cm_heat(M,m,G,T,Tc,pmax=pmax,kmax=kmax,
                                nr=n,nang=n,nphi=n) for n in nr}
    q = {str(n):transfer_heat(M,m,G,T,Tc,xmax=xmax,n=n,nq=nq)
         for n in ne}
    denom = G**2*m**9
    return {"T_over_me":temperature_ratio, "Tchi_over_Te":chi_temperature_ratio,
            "M_over_me":.1, "G_GeV_minus2":G, "me_GeV":m,
            "incoming_kinetic_cut_over_T":xmax,
            "fiveD_CM_Q_over_G2_me9":{n:v/denom for n,v in reference.items()},
            "ell_Q_over_G2_me9":{n:v/denom for n,v in q.items()},
            "finest_relative_difference":q[str(ne[-1])]/reference[str(nr[-1])]-1,
            "limits":"Tree vector EFT; common zero chemical potentials; not a real coupled SM/FRW relic"}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--temperature-ratio", type=float, default=1.)
    p.add_argument("--chi-temperature-ratio", type=float, default=.9)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--output", default="")
    o = p.parse_args()
    if o.temperature_ratio<=0 or o.chi_temperature_ratio<=0:
        p.error("positive temperatures required")
    data = benchmark(temperature_ratio=o.temperature_ratio,
                     chi_temperature_ratio=o.chi_temperature_ratio,
                     nr=(12,16) if o.quick else (24,32,40),
                     ne=(24,32) if o.quick else (64,96,160))
    if o.output:
        path = Path(o.output)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(data,indent=2))
