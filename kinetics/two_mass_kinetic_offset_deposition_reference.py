"""Positive event-level energy deposition and a two-species trace error bound."""
import json

import numpy as np
from scipy.special import roots_legendre

DT=np.longdouble


def kinetic(p2,m):
    return p2/(np.sqrt(m*m+p2)+m)


def event(e,m,b,P,mu):
    e,m,b,P,mu=map(DT,(e,m,b,P,mu))
    if min(e,m,b)<=0 or P<0 or abs(mu)>1:
        raise ValueError("Positive daughter masses, parent excess and admissible boost/angle required.")
    M=m+b+e
    p2=e*(2*(m+b)+e)*(2*m+e)*(2*b+e)/(4*M*M)
    p=np.sqrt(p2)
    K=kinetic(P*P,M)
    bg=P/M
    g=1+K/M
    ds=g*p*mu+bg*np.sqrt(m*m+p2)
    df=-g*p*mu+bg*np.sqrt(b*b+p2)
    pt=p2*(1-mu)*(1+mu)
    z=np.array([ds*ds+pt,df*df+pt],dtype=DT)
    return z,kinetic(z,np.array([m,b],dtype=DT)),e+K


def shell(p,m,order=20):
    x,w=roots_legendre(order)
    lo,hi=np.asarray(p,dtype=DT)[:-1,None],np.asarray(p,dtype=DT)[1:,None]
    h=(hi-lo)
    if np.any(h<=0) or np.any(lo<0):
        raise ValueError("Strict momentum shells are required.")
    q=lo+h*(1+np.asarray(x,dtype=DT))/2
    W=q*q*np.asarray(w,dtype=DT)
    W/=W.sum(axis=1,keepdims=True)
    K=kinetic(q*q,DT(m))
    mean=np.sum(W*K,axis=1)
    curvature=np.sum(W*K*K/(DT(m)+K),axis=1)
    variance=np.sum(W*(K-mean[:,None])**2,axis=1)
    return mean,curvature,variance


def deposit(e,m,b,P,mu,grids):
    _p2,K,available=event(e,m,b,P,mu)
    arms=[]
    out=[]
    for v,h,z in zip(K,map(DT,(m,b)),grids,strict=True):
        mean,curvature,var=shell(z,h)
        i=int(np.searchsorted(mean,v,side='right')-1)
        if not 0<=i<len(mean)-1:
            raise ValueError("Outside the convex shell representatives; do not silently deposit.")
        frac=(v-mean[i])/(mean[i+1]-mean[i])
        w=np.array([1-frac,frac],dtype=DT)
        if min(w)<0 or abs(sum(w)-1)>1e-17:
            raise ValueError("Nonpositive or nonconservative deposition.")
        interp=w@mean[i:i+2]
        diff=w@curvature[i:i+2]-v*v/(h+v)
        bound=((v-mean[i])*(mean[i+1]-v)+w@var[i:i+2])/h
        if not 0<=diff<=bound*(1+1e-7) or abs(interp/v-1)>2e-16:
            raise ValueError(f"Convex stress defect outside bound: {diff}, {bound}.")
        arms.append(w)
        out.append(dict(mean_pair_K=list(map(float,mean[i:i+2])),
                        convex_weights=list(map(float,w)),
                        trace_bias=float(diff),trace_bias_certified_upper=float(bound),
                        projected_K=float(interp)))
    paired=np.outer(*arms)
    projected=sum(DT(z['projected_K']) for z in out)
    err=float(abs(projected-available)/available)
    if err>2e-15 or min(paired.ravel())<0 or abs(paired.sum()-1)>2e-17:
        raise ValueError("Paired number or stabilized event-energy check failed.")
    return dict(event_daughter_excess=list(map(float,K)),parent_energy_above_daughter_rest=float(available),
                energy_relative_error=err,paired_bin_number_weights=np.asarray(paired,float).tolist(),
                per_daughter=out)


def mp_kinetic(e,m,b,P,mu):
    import mpmath as mp
    with mp.workdps(85):
        E,s,f,p,c=map(lambda x:mp.mpf(float(x)),(e,m,b,P,mu))
        M=s+f+E
        parent=mp.sqrt(M*M+p*p)
        Es=(M*M+s*s-f*f)/(2*M)
        Ef=(M*M+f*f-s*s)/(2*M)
        q=mp.sqrt((Es-s)*(Es+s))
        return np.array([float((parent*Es+p*q*c)/M-s),float((parent*Ef-p*q*c)/M-f)])


def test():
    cases=[('ordinary',.093,.061,.122,.33,.23),
           ('near_threshold',2e-23,.061,.122,4e-12,-.31),
           ('backward_near_rest',.031,.065,.018,None,-1+2e-16)]
    for name,e,m,b,P,mu in cases:
        if P is None:
            p2,*_=event(e,m,b,0.,0.)
            P=float(np.sqrt(p2[0])*(DT(e)+DT(m)+DT(b))/DT(m))
        p2,K,E=event(e,m,b,P,mu)
        p=np.sqrt(p2)
        z=[np.geomspace(s/5,s*5,13) for s in p]
        a=deposit(e,m,b,P,mu,z)
        confirm=mp_kinetic(e,m,b,P,mu)
        error=max(abs(np.array(K,float)/confirm-1))
        if error>6e-10:
            raise ValueError(f"Independent high-precision boosted collision failed for {name}.")
        a.update(name=name, inputs=[e,m,b,P,mu], raw_absolute_daughter_energy_rounded_to_rest=[
            float(m+K[0])==m,float(b+K[1])==b],high_precision_reference_relative_error=float(error))
        yield a


if __name__=='__main__':
    print(json.dumps(dict(scope="Independent positive event deposition; no collision kernels or evolving source.",
                          cases=list(test())),indent=2,allow_nan=False))
