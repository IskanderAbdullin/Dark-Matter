"""Check a detailed-balance failure when broadening only the shell constraint."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit

from early_pole_decay_mode_clock_reference import Q, intervals


def evaluate(es, b0, temp, minus, plus, slope):
    # This is a manufactured positive mismatch measure, not a conjecture
    # about the QED line shape, decay rate, or plasma width.
    n = 1/np.expm1(es/temp)
    ny = expit(-b0/temp)
    x0 = b0+es

    def turnover(d, lam=None):
        if lam is None:
            ep, eb=x0+d,b0
        else:
            ep,eb=x0+(1-lam)*d,b0+(1-lam)*d
        p,b=expit(-ep/temp),expit(-eb/temp)
        g=(1+n)*p*(1-b)
        l=n*b*(1-p)
        return g+l,g-l

    def phase(d):
        return 1+slope*d/(minus+plus)

    p0=expit(-x0/temp)
    H0=turnover(0.)[0]
    Hp0=-(1-p0)*p0*((1+n)*(1-ny)-n*ny)/temp
    exact_D0=-H0/(2*temp)
    analytic_D2=-(Hp0+slope*H0/(minus+plus))/(2*temp)

    def net(d):
        return -phase(d)*turnover(d)[0]*np.tanh(d/(2*temp))

    def divided(d):
        if abs(d)<2e-7*temp:
            return analytic_D2
        return (net(d)-d*exact_D0)/d**2

    # Parsed in the same positive bath-frame energy domain on both sides.
    det=quad(divided,-minus,0.,epsabs=7e-8,epsrel=4e-10)[0]+quad(
        divided,0.,plus,epsabs=7e-8,epsrel=4e-10)[0]
    linear_limit=(det+exact_D0*np.log(plus/minus))/(np.pi*H0)
    checks=[]
    for gamma in (2e-4, 2e-5, 2e-6, 2e-7, 2e-8):
        al,ar=np.arctan(-minus/gamma),np.arctan(plus/gamma)
        radial=sorted({0.,*[a for x in (1.,10.,100.,1000.,10000.,100000.,1000000.) for a in (-np.arctan(x),np.arctan(x)) if al<a<ar]})
        i0=quad(lambda z:phase(gamma*np.tan(z))*turnover(gamma*np.tan(z))[0],
                al,ar,points=radial,epsabs=1e-9,epsrel=4e-11)[0]/np.pi
        # Integrate the linear term exactly; the remainder is smooth
        # on the ordinary energy coordinates without a narrow peak.
        il=exact_D0*gamma/(2*np.pi)*np.log((plus*plus+gamma*gamma)/(minus*minus+gamma*gamma))
        ir=gamma/np.pi*sum(quad(lambda d:divided(d)*d*d/(d*d+gamma*gamma),
                                a,b,epsabs=1e-8,epsrel=3e-10,limit=160)[0]
                           for a,b in [(-minus,-gamma),(-gamma,0.),(0.,gamma),(gamma,plus)])
        i1=il+ir
        # Independent direct evaluation in the natural Lorentzian angle.
        direct=quad(lambda z:net(gamma*np.tan(z)),al,ar,
                    epsabs=1e-12,epsrel=7e-10,points=radial,limit=160)[0]/np.pi
        if abs(direct-i1)>2e-11*max(i0,1.):
            raise ArithmeticError(('independent Lorentzian convolution',gamma,direct,i1))
        defects=[]
        for lam in (0.,.25,.75,1.):
            residual=quad(lambda z:phase(gamma*np.tan(z))*turnover(gamma*np.tan(z),lam)[1],
                          al,ar,epsabs=3e-9,epsrel=1e-9,points=radial)[0]/(np.pi*i0)
            defects.append(abs(residual))
        if max(defects)>3e-14:
            raise ArithmeticError(('off-shell conserving frequencies did not balance',defects))
        checks.append(dict(gamma=gamma, pole_frozen_net_over_turnover=i1/i0,
                           pole_frozen_net_over_turnover_per_gamma=i1/i0/gamma,
                           independent_net_difference_over_turnover=abs(i1-direct)/i0,
                           off_shell_matched_max_net_over_turnover=max(defects)))
    end=checks[-1]['pole_frozen_net_over_turnover_per_gamma']
    if abs(end-linear_limit)>4e-4*max(abs(linear_limit),1.):
        raise ArithmeticError(('finite-window linear limit',linear_limit,end))
    pointwise=[]
    for d in [-minus,-minus/2,-1e-6,1e-6,plus/2,plus]:
        w,net0=turnover(d)
        x=-np.tanh(d/(2*temp))
        if abs(net0/w-x)>2e-14:
            raise ArithmeticError('Independent FD odds calculation did not close')
        pointwise.append(dict(energy_mismatch=d, raw_ratio=float(net0/w), expected=float(x)))
    return dict(T=temp,window=[-minus,plus],linear_phase_slope=slope,
                extrapolated_small_gamma_ratio_per_gamma=float(linear_limit),finite_width=checks,
                points=pointwise)


def run():
    fixtures=[]
    for g in ('fine','third'):
        with gzip.open(Path(__file__).parent/'results'/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:
            d=json.load(f)
        M=d['x_end']
        m=M*d['mass_ratio']
        s=d['checkpoints'][-1]['scalar_pole_over_T']
        E,lo,hi,span,_=intervals(M,s,m,Q)
        eb=lo+span/2-E
        rows=[evaluate(E,eb,t,lo0,hi0,w) for t in (.1,1.) for lo0,hi0,w in
              ((.02,.02,0.),(.02,.02,.7),(.006,.024,0.))]
        fixtures.append(dict(grid=g,scalar_E=E,fixed_reference_daughter_E=eb,controls=rows))
    return dict(scope='Manufactured positive Lorentzian mismatch test: comparing occupation weights evaluated on conserving frequencies to those frozen at peak energies. Not an archive width or a QED spectral ansatz.',data=fixtures)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
