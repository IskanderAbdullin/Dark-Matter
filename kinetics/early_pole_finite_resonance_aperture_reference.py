"""Static broadened angular aperture with a thermal detailed-balance audit."""
import argparse
import json

import numpy as np
from scipy.integrate import quad
from scipy.special import erf, erfinv, expit, roots_legendre

from early_pole_decay_mode_clock_reference import Q, intervals

X, MS, RB, G = .05, .0008672, .06, .005
MP, MB, ES = np.hypot(X, G/2), np.hypot(RB*X, G/2), np.hypot(MS, Q)
_, EL, EH, _, _ = intervals(MP, MS, MB, Q)
PL, PH = np.sqrt((np.array([EL, EH])-MP)*(np.array([EL, EH])+MP))
PC = .5*(PL+PH)


def defects(p):
    ep = np.hypot(p, MP)
    b = np.sqrt((p-Q)**2+MB**2), np.hypot(p+Q, MB)
    # Rationalize Ep-Eb; avoids subtracting two hard energies on the soft scale.
    n = MP*MP-MB*MB-Q*Q
    return np.array([(n-2*p*Q)/(ep+b[1]), (n+2*p*Q)/(ep+b[0])])-ES


def dcdf(d, shape):
    if shape == 'lorentz':
        return np.arctan(d)/np.pi
    return erf(d/np.sqrt(2))/2


def quantile(t, shape):
    if shape == 'lorentz':
        return np.tan(np.pi*t)
    return np.sqrt(2)*erfinv(2*t)


def calc(shape, width, t, gap, cut, n):
    epc = np.hypot(PC, MP)
    z, w = roots_legendre(n)
    # The p^2 dp dcos/(Ep Eb) measure reduces exactly to p dp dd/(q Ep).
    # Return q times the integral to avoid a pure common factor.
    def row(p, moment):
        ep = np.hypot(p, MP)
        lo, hi = defects(p)
        ua, ub = dcdf(np.array([lo, hi])/width, shape)
        if ub == ua:
            return 0.
        u = (ua+ub)/2+(ub-ua)*z/2
        d = np.clip(width*quantile(u, shape), lo, hi)
        a = expit(-(ep-gap)/t)
        f = a*expit((ep-ES-d)/t)
        if moment == 'phase':
            v = np.ones_like(d)
        elif moment == 'repaired':
            v = f*np.exp(d/(2*t))
        else:
            # equilibrium at gap=ES. Retaining central-pole occupations at
            # d!=0 while broadening only delta fails pointwise detailed balance.
            v = f*(-np.expm1(d/t))
        return p/ep*(ub-ua)*(.5*np.dot(w, v))
    grid = [0., cut]+[j for j in (.1, .3, 1., 3., PC, PL, PH) if 0<j<cut]
    if cut>PC:
        for j in [PL, PC, PH]:
            for h in (0.1,1,10,100):
                for sign in (-1,1):
                    p = j+sign*h*(PH-PL)*width/(defects(PC)[1]-defects(PC)[0])
                    if 0<p<cut:
                        grid.append(p)
    grid = sorted(set(grid))
    ans = [sum(quad(lambda p:row(p,k), a,b,epsabs=2e-13,
                    epsrel=2e-7,limit=120)[0] for a,b in zip(grid[:-1],grid[1:],strict=True))
           for k in ('phase','repaired','naive')]
    delta = quad(lambda p:p/np.hypot(p,MP),PL,min(PH,cut),epsabs=2e-13)[0] if cut>PL else 0.
    delta_seed = quad(lambda p:(p/np.hypot(p,MP))*expit(-(np.hypot(p,MP)-gap)/t)
                      *expit((np.hypot(p,MP)-ES)/t),PL,min(PH,cut),epsabs=2e-13)[0] if cut>PL else 0.
    # Detailed balance factor is independent of p and d at common T. Use a
    # representative d to check the separate direct forward/reverse expression.
    d = .002*np.linspace(-1,1,17)
    a, b = expit(-(epc-gap)/t), expit(-(epc-ES-d)/t)
    f, r = a*(1-b)*np.exp(d/(2*t)), b*(1-a)*np.exp(-d/(2*t))
    thermal_defect = np.max(np.abs(np.log(r/f)-(ES-gap)/t))
    assert thermal_defect<3e-14
    return dict(shape=shape,width=width,temperature=t,gap=gap,pmax=cut,
                center_angular_fraction=float(np.diff(dcdf(defects(PC)/width,shape))[0]),
                radial_phase=ans[0],radial_seed_thermal=ans[1],
                radial_naive_stimulus=ans[2],
                radial_thermal_stimulus=ans[1]*(-np.expm1((ES-gap)/t)),
                delta_phase=delta,delta_seed=delta_seed,pointwise_log_balance_error=float(thermal_defect))


def run(n):
    split = np.diff(defects(PC))[0]
    ans = [dict(ratio=r,**calc(s,r*split, .8, ES, cut,n))
           for s in ('gaussian','lorentz') for r in (.015,.1,.5,1.,2.,5.)
           for cut in (1.,4.,12.)]
    for a in ans:
        assert a['radial_thermal_stimulus']==0.
        if a['ratio']==.015 and a["shape"]=="gaussian" and a["delta_seed"]>0.:
            assert abs(a['radial_seed_thermal']/a['delta_seed']-1)<1e-3
    gate = [calc(s,split,.8,ES+d*ES,4.,n)
            for s in ('gaussian','lorentz') for d in (-.1,.1)]
    assert all(np.sign(a['gap']-ES)==np.sign(a['radial_thermal_stimulus']) for a in gate)
    return dict(x=X,scalar_mass=MS,parent_mass=MP,daughter_mass=MB,q=Q,scalar_energy=ES,
                resonant_parent_momentum=[PL,PH], angular_defect_span=split,
                rows=ans,thermal_sign_checks=gate)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nodes',type=int,default=96)
    print(json.dumps(run(p.parse_args().nodes),indent=2,allow_nan=False))
