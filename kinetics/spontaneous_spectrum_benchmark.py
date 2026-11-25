"""Analytic scalar spontaneous-production spectrum for a massless scalar.

Parent FD at theta=1, empty daughter fermions and scalar, vacuum isotropic
decay. E1min=q/eta+eta*x^2/(4*q), eta=1-r^2.
Integrating FD(E1) from E1min to infinity gives log(1+exp(-E1min)).
This is a local independent benchmark, not an abundance calculation.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre,expit
from parent_angle_kernel import ParentAngle,SeparatedRectangle
from resolved_parent_bath import ResolvedParent


def reference(m,x,theta=1.,order=64):
    u,w=roots_legendre(order)
    l,h=m.edges[:-1],m.edges[1:]
    q=(l[:,None]+h[:,None])/2+(h-l)[:,None]*u/2
    eta=1-m.fermion_ratio**2
    E=q/eta+eta*x*x/(4*q)
    integrand=theta*np.logaddexp(0.,-E/theta)
    return m.g_fermion*x/(2*np.pi**2)*m.kappa*x*x/eta*(h-l)/2*(integrand@w)


def screen(m,x,theta=1.):
    m.geometry(x);g=expit(-m.e1/theta)
    c=m.dark(np.r_[np.zeros(2*m.n),g])[:m.n]
    actual=m.vol*c;ref=reference(m,x,theta);rate=ref/m.vol
    sel=rate>.1*np.max(rate)
    # Scalar number-weighted L1 and growth-rate errors in stimulated region.
    rel=np.abs(c[sel]-rate[sel])/rate[sel]
    q=m.mean
    return dict(x=x,total_relative_error=float(actual.sum()/ref.sum()-1),
        L1_number=float(np.abs(actual-ref).sum()/ref.sum()),
        max_relative_growth_error_peak_band=float(rel.max()),
        peak_growth_relative_error=float(c.max()/rate.max()-1),
        peak_q_actual=float(q[np.argmax(c)]),peak_q_reference=float(q[np.argmax(rate)]),
        actual_peak=float(c.max()),reference_peak=float(rate.max()),
        equilibrium_dropped_fraction=m.max_dropped_equilibrium_fraction)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True)
    a=p.parse_args();rows=[]
    for n in [80,120,240,320,480]:
        kw=dict(n=n,nparent=n,kappa=.006,scalar_ratio=1e-12,fermion_ratio=.9,qmin=1e-8,
                parent_floor=0.,parent_grid='uniform')
        for name,model in [('legacy-shared',ResolvedParent(**kw)),
                           ('rectangle-separated',SeparatedRectangle(fermion_grid='balanced',**kw)),
                           ('angle-separated',ParentAngle(fermion_grid='balanced',**kw))]:
            for x in [.1,.2,.3,.5,1.,2.,5.]:
                rows.append(dict(n=n,kernel=name,**screen(model,x)))
            selected=[r for r in rows if r['n']==n and r['kernel']==name]
            print(name,n,'max peak growth bias',max(abs(r['peak_growth_relative_error']) for r in selected),
                'max L1',max(r['L1_number'] for r in selected),flush=True)
    Path(a.output).write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
