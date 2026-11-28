"""Disentangle cold scalar weight integration from bare kernel node sampling."""
import gzip,json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

D=np.longdouble


def audit(d,w):
    M=D(d['x_end']);m=M*D(d['mass_ratio']);s=D(d['checkpoints'][-1]['scalar_pole_over_T'])
    A=(M*M-m*m+s*s)/2
    W=np.sqrt((A-M*s)*(A+M*s))
    Z=(M*M+m*m-s*s)/2
    C=(M+m)**2/(Z+M*m)
    P=-np.pi*C*W/(4*s*s)
    L=np.pi*Z/(2*W)
    p=np.hypot(M,np.asarray(d['endpoint_parent_momentum_edges'],dtype=D))
    daughter_edges=np.hypot(m,np.asarray(d['endpoint_fermion_momentum_edges'],dtype=D))
    fp=np.asarray(d['endpoint_parent_occupations'],dtype=D)
    fb=np.asarray(d['endpoint_fermion_occupations'],dtype=D)
    qe=np.asarray(d['scalar_momentum_edges'])
    fq=np.asarray(d['endpoint_scalar_occupations'])
    qmin,qmax=1e-6,1e-5
    def get(e,o,x):
        j=np.searchsorted(e,x,side='right')-1
        if not 0<=j<len(o): raise ArithmeticError('Not on the archived occupation domain.')
        return o[j]

    def pair(q):
        q=D(q)
        es=np.hypot(s,q)
        e=es*A/s**2
        half=q*W/s**2
        pp,bb=(p-e)/half,(daughter_edges+es-e)/half
        sep=sorted({-1.,0.,1.,*[x for x in (*pp,*bb) if -1<x<1]})
        ans=[]
        for lo,hi in zip(sep[:-1],sep[1:],strict=True):
            z=(lo+hi)/2
            f=get(pp,fp,z)-get(bb,fb,z)
            ans.append(f*(P*(hi-lo)+L*(-1 if hi<=0 else 1)*np.log((abs(hi) or 1)/(abs(lo) or 1))))
        return fsum(ans)
    err=max(abs(pair(x['q'])-x['common_aperture_scaled_finite_part']) for x in w['momentum_nodes'])
    if err>2e-8:raise ArithmeticError('New cold integration did not reproduce the Q scan.')

    def parts(a,b,n,weighted):
        X,W=roots_legendre(n)
        knots=sorted({a,b,*[x for x in qe if a<x<b],
                      *[z['q'] for z in w['first_aperture_crossings'] if a<z['q']<b]})
        norm=[]; area=[]
        for x,y in zip(knots[:-1],knots[1:],strict=True):
            occ=get(qe,fq,(x+y)/2)
            e=x+(y-x)/2*(1+X)
            masses=(y-x)/2*W*occ*e*e/np.hypot(e,float(s))
            norm.append(fsum(masses))
            if weighted:area.append(fsum(m*pair(z) for m,z in zip(masses,e,strict=True)))
        return fsum(norm),fsum(area)
    a,b=parts(qmin,qmax,17,True),parts(qmin,qmax,31,True)
    ref=b[1]/b[0]
    if max(abs(a[0]/b[0]-1),abs(a[1]/b[1]-1))>5e-10:
        raise ArithmeticError('Split control did not converge.')
    out=[]
    for n in (8,16,32,64,128):
        e=np.geomspace(qmin,qmax,n+1)
        nodes=np.cbrt((e[1:]**3+e[:-1]**3)/2)
        kernel=[pair(q) for q in nodes]
        # Integrate q²/E without occupation for the ordinary node sample.
        z,wgt=roots_legendre(12)
        S=[]
        T=[]
        for l,r,q,k in zip(e[:-1],e[1:],nodes,kernel,strict=True):
            u=l+(r-l)/2*(1+z)
            S.append(get(qe,fq,q)*(r-l)/2*(wgt@(u*u/np.hypot(u,float(s)))))
            T.append(parts(l,r,12,False)[0])
        norm1=fsum(S)
        norm2=fsum(T)
        r1=fsum(x*y for x,y in zip(S,kernel,strict=True))/norm1
        r2=fsum(x*y for x,y in zip(T,kernel,strict=True))/norm2
        if not (min(kernel)<=r1<=max(kernel) and min(kernel)<=r2<=max(kernel)):
            raise ArithmeticError('A positive nodal weighted average left the convex hull.')
        if abs(norm2/b[0]-1)>5e-12:
            raise ArithmeticError('Exact scalar shell quadrant weights did not preserve full intensity.')
        out.append(dict(cold_cells=n,node_weighted_average=r1,exact_cell_scalar_weight_average=r2,
                        relative_integrated_scalar_weight_error=norm1/b[0]-1,
                        default_node_relative_error=r1/ref-1,exact_cell_weight_relative_error=r2/ref-1))
    return dict(full_split_weighted_average=ref,scalar_weight_integral=b[0],mesh_defect=abs(a[1]/a[0]-ref),
                input_pointwise_replay_defect=err,cold_grid_ladder=out)


def run():
    p=Path(__file__).parent/'results'
    with open(p/'finite-mass-qed-shared-polar-momentum-scan-2026-11-27.json') as f:w=json.load(f)
    out=[]
    for g in ('fine','third'):
        with gzip.open(p/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        out.append(dict(grid=g,**audit(d,next(x for x in w['rows'] if x['grid']==g))))
    return dict(scope='Numerical sampling audit of a convention-dependent aperture-scale finite part. Cold scalar weight q² f_s/E_s on 1e-6 to 1e-5. Not a physical rate.',rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
