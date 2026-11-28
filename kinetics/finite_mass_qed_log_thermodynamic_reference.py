"""Two independent quadratures of the finite-mass O(g²) QED log term."""
import json

import numpy as np
from scipy.integrate import quad
from scipy.special import expit, roots_legendre


def nodes(n):
    x,w=roots_legendre(n)
    return (x+1)/2,w/2


def pressure(x, n=128, cap=48.):
    """Dimensionless p^(2)/g²T⁴ and (rho-3p)^(2)/g²T⁴; x=m/T."""
    t,w=nodes(n)
    edges=np.array([0.,.001,.004,.016,.064,.256,1.,4.,12.,24.,cap])
    wid=np.diff(edges)
    p=(edges[:-1,None]+wid[:,None]*t).ravel()
    wp=(wid[:,None]*w).ravel()
    ep=np.hypot(p,x)
    fp=expit(-ep)
    # Symmetric triangle r<p, using 1-r/p=t^4 at the log endpoint.
    r=p[:,None]*(1-t**4)
    wr=p[:,None]*(4*t**3*w)
    er=np.hypot(r,x)
    fr=expit(-er)
    ell=np.log1p(1-t**4)-4*np.log(t)
    q=(wp*(p/ep)*fp)[:,None] * wr*(r/er)*fr * ell
    J=wp @ ((p*p/ep)*fp)
    Jprime=wp @ (-(p*p/ep)*fp*(x/(ep*ep)+x*(1-fp)/ep))
    ip0=-(np.pi**2*J/6+J**2/2)/np.pi**4
    iplog=x*x*q.sum()/(2*np.pi**4)
    # rho-3p = -x dp/dx; differentiation at fixed T, including FD weights.
    tr0=(np.pi**2/6+J)*Jprime*x/np.pi**4
    qp=(2-x*x*(1/ep**2+(1-fp)/ep)[:,None]
        -x*x*(1/er**2+(1-fr)/er))
    trlog=-x*x*(q*qp).sum()/(2*np.pi**4)
    return np.array((ip0,iplog,tr0,trlog))


def shell_integral(p, x, cap):
    """Independent weighted quadrature: log(p+r) minus log|p-r|."""
    def amp(r):
        en=np.hypot(r,x)
        return (r/en)*expit(-en) if en else .5
    if p==0:
        # Return the p→0 limit of the log integral divided by p.
        return 2*quad(lambda r: expit(-np.hypot(r,x))/np.hypot(r,x),
                      0,cap,epsabs=1e-13,epsrel=1e-13)[0]
    def regular(r):
        return amp(r)*np.log(p+r)
    plus=quad(regular,0,cap,points=[p],epsabs=2e-13,epsrel=5e-13)[0]
    minus=0.
    # Use analytic log singular weights at the lower endpoint of both segments.
    for sign,end in ((-1,p),(1,cap-p)):
        if end>0:
            minus+=quad(lambda t: amp(p+sign*t),0,end,weight="alg-loga",
                        wvar=(0.,0.),epsabs=2e-13,epsrel=5e-13)[0]
    return (plus-minus)/p


def from_poles(x,cap=48.):
    """Use the independently integrated fermion/photon on-shell energy shifts."""
    def measure(p):
        en=np.hypot(p,x)
        return p*p*expit(-en)/en
    J=quad(measure,0,cap,epsabs=2e-13,epsrel=5e-13)[0]
    # Remove O(g²) coupling factor, and integrate 1/2 ∫ state ΔE.
    photon_m2=2*J/np.pi**2
    fermion_m2=(np.pi**2/6+J)/np.pi**2
    nonlog=-(photon_m2*(np.pi**2/6) + 2*fermion_m2*J)/(4*np.pi**2)
    logged=0. if x==0 else quad(lambda p: measure(p)*shell_integral(p,x,cap),0,cap,
                epsabs=2e-12,epsrel=2e-11)[0]*x*x/(4*np.pi**4)
    return np.array((nonlog,logged))


def run():
    out=[]
    for x in (0.,.003,.05,.25,.5,1.,2.,5.):
        p=pressure(x)
        delta=np.max(np.abs(p-pressure(x,192)))
        tail=np.max(np.abs(p-pressure(x,192,56.)))
        if delta>6e-13 or tail>6e-13:
            raise ValueError('Triangle mesh or momentum tail not resolved.')
        v=from_poles(x)
        defect=np.max(abs(p[:2]-v))
        if defect>5e-13:
            raise ValueError('Independent on-shell pressure mismatch.')
        if x>0:
            e=2e-4
            samples=[pressure(x*np.exp(j*e))[:2] for j in (-2,-1,1,2)]
            trace=-(samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*e)
            traceerr=np.max(abs(trace-p[2:]))
        else:
            traceerr=np.max(abs(p[2:]))
            if abs(p[0]+5/288)>2e-14 or p[1]!=0.:
                raise ValueError('Incorrect massless pressure limit.')
        if traceerr>2e-11:
            raise ValueError('Thermodynamic trace identity failed.')
        # Profile the log relative to the total factorizable fermion pole shift.
        small=(0. if x==0 else shell_integral(0.,x,48.))
        candidate=[0.05,.25,1.,4.]
        J=quad(lambda k: k*k*expit(-np.hypot(k,x))/np.hypot(k,x),0,48.)[0]
        shifts=[-x*x*shell_integral(k,x,48.)/(2*(np.pi**2/6+J)) if x else 0. for k in candidate]
        out.append(dict(m_over_T=x,pressure_factorizable=float(p[0]),
                        pressure_log=float(p[1]),trace_factorizable=float(p[2]),
                        trace_log=float(p[3]),mesh_defect=float(delta),tail_defect=float(tail),
                        pressure_vs_integrated_poles_defect=float(defect),
                        trace_derivative_defect=float(traceerr),
                        pole_momentum_over_T=candidate,
                        pole_log_over_total_factorizable=shifts,
                        pole_rest_log_over_total_factorizable=-x*x*small/(2*(np.pi**2/6+J))))
    return dict(scope='One unpolarized Dirac species and one photon at equilibrium, charge symmetric; O(g²).',
                rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
