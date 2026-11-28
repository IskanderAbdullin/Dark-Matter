"""Regular negative-frequency single-cut integrals in a bare scalar vertex."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import intervals,Q


def one(d,grid,cut):
    M,m=d['x_end'],d['mass_ratio']*d['x_end']
    s=d['checkpoints'][-1]['scalar_pole_over_T']
    Es,lo,hi,*_=intervals(M,s,m,Q)
    ep=(lo+hi)/2
    P=np.sqrt((ep-M)*(ep+M))
    z=(ep*Es-(M*M+s*s-m*m)/2)/Q
    if cut=='parent':
        mass,other,e,qz,sg=M,m,ep,z,1.
        name='parent'
    else:
        mass,other,e,qz,sg=m,M,ep-Es,z-Q,-1.
        name='fermion'
    p=np.sqrt((e-mass)*(e+mass))
    c=qz/p
    z0=(M*M+m*m-s*s)/2
    C=(M+m)**2/(z0+M*m)
    D=4-2*(M+m)*mass/(z0+M*m)
    dif=mass*mass+s*s-other*other
    ce=dif+2*sg*mass*s
    if sg*ce<=0:
        raise ValueError('The bare negative-frequency gap is not protected.')
    eta=np.arcsinh(p/mass)

    def quantities(r,mu):
        er=np.hypot(mass,r)
        pr=p*r
        a=2*mass*mass*np.cosh((np.arcsinh(r/mass)-eta)/2)**2
        # Negative cut R=(-E_r,r*direction). This is a future timelike photon.
        Df=dif+2*sg*(Es*er+Q*r*mu)
        disc=a*a+2*a*pr*(1+c*mu)+(pr*(mu+c))**2
        ph=1/(2*np.sqrt(disc))
        if (np.min(sg*Df)<sg*ce*(1-1e-13)
                or np.min(disc)<4*mass**4*(1-1e-13)):
            raise ValueError('A would-be negative-frequency denominator became singular.')
        return (C+4*z0*ph)/Df+D*ph/2

    def angle(r,n):
        t,w=roots_legendre(n)
        pr=p*r
        a=2*mass*mass*np.cosh((np.arcsinh(r/mass)-eta)/2)**2
        aa=np.sqrt(2*a*(1-c*c)/pr)
        L=np.arcsinh((-1+c)/aa)
        H=np.arcsinh((1+c)/aa)
        z=(L+H)/2+(H-L)*t/2
        u=aa*np.sinh(z)-c
        wt=w*(H-L)/2*aa*np.cosh(z)
        return np.sum(wt*quantities(r,u),axis=-1)

    az=[]
    for r,u in ((.2,-.7),(1.43,.18),(3.,.55),(8.,-.4)):
        er=np.hypot(mass,r)
        shift=dif+2*sg*(Es*er+Q*r*u)
        base=2*(mass*mass+e*er+p*r*c*u)
        dev=2*p*r*np.sqrt((1-c*c)*(1-u*u))
        def phi(v):
            inv=1/(base+dev*np.cos(v))
            return (C+4*z0*inv)/shift+D*inv/2
        val=quad(phi,0,2*np.pi,epsabs=2e-10,epsrel=2e-12)[0]/(2*np.pi)
        az.append(abs(quantities(np.asarray(r),np.asarray(u))/val-1))
    if max(az)>2e-12:
        raise ValueError('Negative-cut azimuth trace did not close.')

    defects=[]
    for r in (.08,.3,1.43,3.,9.):
        x=np.asarray(r)
        a=angle(x,52)
        width=mass/np.sqrt(p*x)
        pts=[v for v in [-c,*[-c+j*width for j in (-4,-1,1,4)]] if -1<v<1]
        b=quad(lambda u:float(quantities(x,u)),-1.,1.,points=pts,epsabs=2e-9,epsrel=1e-11)[0]
        defects.append(abs(a/b-1))
    if max(defects)>1e-9:
        raise ValueError('Negative cut angle integration disagrees with adaptive quadrature.')

    edges=np.asarray(d[f'endpoint_{name}_momentum_edges'])
    f=np.asarray(d[f'endpoint_{name}_occupations'])
    def radial(n,o,occupied):
        t,w=roots_legendre(n)
        dx=np.diff(edges)/2
        p=(edges[:-1]+dx)[:,None]+dx[:,None]*t
        mu=p.reshape(-1,1)
        a=angle(mu,o).reshape(p.shape)
        raw=np.sum(w*dx[:,None]*p*p/np.hypot(p,mass)*a/2,axis=-1)
        return (raw*occupied).sum()
    ref=radial(14,52,f)
    a=radial(8,32,f)
    if abs(a/ref-1)>3e-8:
        raise ValueError('The frozen negative cut did not reach the radial accuracy target.')
    parts=[radial(14,52,f*((edges[1:]<=1.)|(edges[:-1]>=2.))),
           radial(14,52,f*((edges[1:]>1.)&(edges[:-1]<2.)))]
    if abs(sum(parts)-ref)>1e-9:
        raise ValueError('The diagnostic cell grouping did not recover the full cut.')
    return dict(grid=grid,cut=cut,absolute_charged_gap_floor=float(sg*ce),
                positive_photon_mass_sq_floor=4*mass*mass,negative_cut_integral_no_prefactors=float(ref),
                cells_intersecting_one_to_two_fraction=float(parts[1]/ref),
                independent_azimuth_relative_defect=float(max(az)),ordinary_angle_relative_defect=float(max(defects)),radial_mesh_relative_defect=float(abs(a/ref-1)))


def run():
    out=[]
    for grid in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        with gzip.open(path) as f:d=json.load(f)
        for cut in ('parent','daughter'):
            out.append(one(d,grid,cut))
    return dict(scope='Regular negative-frequency bare single cut with charge-symmetric frozen occupation comparator. No thermal real-time matching or full rate.',rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
