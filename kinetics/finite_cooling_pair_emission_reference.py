"""Positive scalar recoil moments and selected CDF bins from vacuum pairs."""
import numpy as np
from scipy.special import roots_legendre

from archived_scalar_pair_reaction_reference import amplitude
from three_species_decay_onsager_reference import m,mp


def arcsine_bins(edges,center,half,weight,mass):
    """Clip a uniform azimuth to energy intervals and integrate p² analytically."""
    b=len(edges)-1
    cut=np.hypot(mass,edges)
    c,h,w=(np.ravel(x) for x in (center,half,weight))
    if not np.all(h>=0.) or not np.all(w>=0.):
        raise ValueError("Negative angular width or quadrature weight")
    bins=np.zeros((2,b))
    fixed=h<=np.finfo(float).eps*c
    if fixed.any():
        pp=np.maximum(0.,(c[fixed]-mass)*(c[fixed]+mass))
        i=np.searchsorted(edges,np.sqrt(pp),side="right")-1
        bins[0]+=np.bincount(i,w[fixed],minlength=b)
        bins[1]+=np.bincount(i,w[fixed]*pp,minlength=b)
    c,h,w=(x[~fixed,None] for x in (c,h,w))
    a=np.clip((cut-c)/h,-1.,1.)
    t=np.arccos(a)
    s=np.sqrt(np.maximum(0.,(1-a)*(1+a)))
    delta=-np.diff(t,axis=1)
    d_s=-np.diff(s,axis=1)
    d_s2=-np.diff(s*a,axis=1)
    moment=((c-mass)*(c+mass)+h*h/2)*delta+2*c*h*d_s+(h*h/2)*d_s2
    assert np.all(moment>=-1e-12*np.maximum(1.,c*c))
    bins[0]+=(w*delta).sum(axis=0)/np.pi
    bins[1]+=(w*np.maximum(0.,moment)).sum(axis=0)/np.pi
    return bins


def kinematics(p,w,order):
    i,j=np.tril_indices(len(p))
    c,g=roots_legendre(order)
    eps=np.hypot(m,p)
    V=p*p/(eps+m)
    pp=p[i,None]*p[j,None]
    E,F=eps[i,None],eps[j,None]
    S=4*m*m+2*(m*(V[i]+V[j])[:,None]+V[i,None]*V[j,None]-pp*c)
    P2=np.maximum(0.,(p[i]**2+p[j]**2)[:,None]+2*pp*c)
    B2=P2/S
    pcm2=np.maximum(0.,S/4-m*m)
    mu=np.clip((E-F)/2/np.sqrt(np.maximum(1e-300,B2*pcm2)),-1.,1.)
    pk2=np.maximum(0.,S/4-mp*mp)
    weight=(w[i]*w[j]*np.where(i==j,1.,2.))[:,None]*g/2
    weight*=np.sqrt(pk2/S)/(64*np.pi*E*F)
    return c,g,S,(E+F)/2,B2,pk2,mu,weight


def kernel_and_lab_p2(p,w,order,edges=None,red=1.):
    c,g,s,center,B2,k2,mu,weight=kinematics(p,w,order)
    moments=np.zeros(2)
    bins=None if edges is None else np.zeros((2,len(edges)-1))
    recoil=np.sqrt(B2*k2)
    for x,gx in zip(c,g):
        blocked=weight*amplitude(s,x)*gx
        C=center+recoil*mu*x
        D=recoil*np.sqrt((1-x*x)*np.maximum(0.,1-mu*mu))
        p2=(C-mp)*(C+mp)+D*D/2
        moments+=np.array([np.sum(blocked),np.sum(blocked*p2)])
        if bins is not None:
            # All weights are nonnegative; both daughters have the same
            # marginal after the even incoming-axis angular integration.
            b,bw,bp=(v.ravel() for v in (C,D,blocked))
            for i in range(0,len(b),4096):
                sel=slice(i,i+4096)
                bins+=arcsine_bins(edges,b[sel]*red,bw[sel]*red,bp[sel],mp*red)
    if bins is not None and moments[0]>0.:
        assert np.max(np.abs(bins.sum(axis=1)/(moments*np.array([1.,red*red]))-1.))<2e-10
    return moments,bins


def checks():
    # Exact count-trace splitting of a single uniform azimuth.
    edges=np.array([0.,.5,.8,1.,1.5,2.5,np.inf])
    c,h,mass=1.2,.6,.2
    v=arcsine_bins(edges,np.array([c]),np.array([h]),np.ones(1),mass)
    q0=np.arccos(np.clip((np.hypot(mass,edges)-c)/h,-1.,1.))
    from scipy.integrate import quad
    ref=np.array([quad(lambda x:(c+h*np.cos(x))**2-mass**2,b,a,
                      epsabs=1e-13,epsrel=1e-13)[0]/np.pi for a,b in zip(q0[:-1],q0[1:])])
    z=np.max(np.abs(v[1]-ref))
    assert z<1e-13

    x,w=roots_legendre(40)
    az=np.cos((np.arange(80)+.5)*np.pi/80)
    error=[]
    for p1,p2,c0 in ((.001,.004,-.7),(.023,.015,.6),(.09,.017,.8),(.18,.14,-.5)):
        E,F=np.hypot(m,[p1,p2])
        s=2*(m*m+E*F-p1*p2*c0)
        boost=np.sqrt(((E+F)**2-s)/s)
        ci=(E-F)/(2*boost*np.sqrt(s/4-m*m))
        d=boost*np.sqrt(s/4-mp*mp)
        C=(E+F)/2
        same=(C+d*ci*x)**2+(d*d)*(1-ci*ci)*(1-x*x)/2-mp*mp
        raw=amplitude(s,x)*w
        one=raw@same/raw.sum()
        rotated=ci*x[:,None]+np.sqrt(1-ci*ci)*np.sqrt(1-x*x)[:,None]*az
        q=amplitude(s,rotated).mean(axis=1)*w
        two=q@((C+d*x)**2-mp*mp)/q.sum()
        error.append(abs(one/two-1.))
    assert max(error)<1e-12
    return dict(azimuth_bin_primitive_absolute_defect=float(z),boosted_p2_rotated_angle_relative_defect=float(max(error)))
