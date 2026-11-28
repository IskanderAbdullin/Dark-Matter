"""Compare one moment-matched node per cold cell against a fully split bare average."""
import gzip,json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

D=np.longdouble


def audit(d,which,scan,prior):
    M=D(d['x_end']);m=M*D(d['mass_ratio']);s=D(d['checkpoints'][-1]['scalar_pole_over_T'])
    A=(M*M-m*m+s*s)/2
    H=np.sqrt((A-M*s)*(A+M*s))
    v=(M*M+m*m-s*s)/2
    P=-np.pi*(M+m)**2/(v+M*m)*H/(4*s*s)
    L=np.pi*v/(2*H)
    pa=np.hypot(M,np.asarray(d['endpoint_parent_momentum_edges'],dtype=D))
    ba=np.hypot(m,np.asarray(d['endpoint_fermion_momentum_edges'],dtype=D))
    fp=np.asarray(d['endpoint_parent_occupations'],dtype=D)
    fb=np.asarray(d['endpoint_fermion_occupations'],dtype=D)
    sk=np.asarray(d['scalar_momentum_edges'])
    so=np.asarray(d['endpoint_scalar_occupations'])
    change=sorted(z['q'] for z in scan['first_aperture_crossings'])
    def kernel(q):
        q=D(q)
        es=np.hypot(q,s)
        center=es*A/s**2
        scale=q*H/s**2
        pe,be=(pa-center)/scale,(ba+es-center)/scale
        z=sorted({-1.,0.,1.,*[p for p in (*pe,*be) if -1<p<1]})
        cols=[]
        for u,v in zip(z[:-1],z[1:],strict=True):
            a=(u+v)/2
            f=fp[np.searchsorted(pe,a)-1]-fb[np.searchsorted(be,a)-1]
            cols.append(float(f*(P*(v-u)+L*(-1 if v<=0 else 1)*np.log((abs(v) or 1.)/(abs(u) or 1.)))))
        return fsum(cols)
    def moments(l,r,n):
        x,w=roots_legendre(n)
        knots=sorted({l,r,*[q for q in sk if l<q<r]})
        a=[]
        for u,v in zip(knots[:-1],knots[1:],strict=True):
            q=u+(v-u)*(x+1)/2
            f=so[np.searchsorted(sk,(u+v)/2)-1]
            g=(v-u)/2*w*f*q*q/np.hypot(q,float(s))
            a.append((fsum(g),fsum(g*q),fsum(g*np.log(q))))
        return [fsum(p) for p in zip(*a,strict=True)]
    ref=prior['full_split_weighted_average']
    den=prior['scalar_weight_integral']
    rows=[]
    for n,old in zip((8,16,32,64,128),prior['cold_grid_ladder'],strict=True):
        e=np.geomspace(1e-6,1e-5,n+1)
        def solve(edges,mode,order):
            a=[];b=[]
            for l,r in zip(edges[:-1],edges[1:],strict=True):
                w,x,y=moments(l,r,order)
                if w<=0: raise ArithmeticError('Expected a positive occupied cold cell.')
                if mode=='default': q=np.cbrt((l**3+r**3)/2)
                elif mode=='linear':q=x/w
                else:q=np.exp(y/w)
                if not l<=q<=r:raise ArithmeticError('A matched centroid escaped its physical cell.')
                a.append(w);b.append(w*kernel(q))
            if abs(fsum(a)/den-1)>2e-12:raise ArithmeticError('The cold scalar weight was not preserved.')
            return fsum(b)/fsum(a)
        z=solve(e,'default',12)
        if abs(z-old['exact_cell_scalar_weight_average'])>3e-10:
            raise ArithmeticError('Did not reproduce the original exact-weight center rule.')
        m1=solve(e,'linear',12)
        m2=solve(e,'log',12)
        split=sorted({*e,*change})
        s1=solve(split,'linear',12)
        s2=solve(split,'linear',26)
        if abs(s1-s2)>3e-10:
            raise ArithmeticError('Independent scalar weights did not resolve the split placement test.')
        rows.append(dict(cold_cells=n,subcell_evaluations=len(split)-1,
                         ordinary_center_relative_error=z/ref-1,weighted_linear_centroid_error=m1/ref-1,
                         weighted_log_centroid_error=m2/ref-1,split_weighted_centroid_error=s1/ref-1))
    return dict(grid=which,fully_split_reference=ref,scalar_weight_integral=den,controls=rows)


def run():
    p=Path(__file__).parent/'results'
    with open(p/'finite-mass-qed-shared-polar-momentum-scan-2026-11-27.json') as f:d=json.load(f)
    with open(p/'finite-mass-qed-shared-polar-cold-mesh-2026-11-27.json') as f:k=json.load(f)
    out=[]
    for grid in ('fine','third'):
        with gzip.open(p/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as f:r=json.load(f)
        out.append(audit(r,grid,next(z for z in d['rows'] if z['grid']==grid),
                         next(z for z in k['rows'] if z['grid']==grid)))
    return dict(scope='Moment-matched placement only for a convention-dependent bare aperture diagnostic; not a kinetic integrator or a physical correction.',grids=out)


if __name__=='__main__':print(json.dumps(run(),indent=2,allow_nan=False))
