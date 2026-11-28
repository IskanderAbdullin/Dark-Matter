"""Positive occupation envelope at fixed daughter two moments and frozen aperture."""
import argparse
import json

import numpy as np
from scipy.optimize import linprog
from scipy.special import expit, roots_legendre

from early_pole_decay_mode_clock_reference import G,KAPPA
from early_pole_daughter_archive_independent_aperture import bounds
from early_pole_daughter_local_fraction_curvature_attribution import initial


def receiver_packets(C,e,j,es,fs,q,width):
    M=C.M
    B=M.cells[j]
    result=[]
    for i in range(M.r[1].start,M.r[1].stop):
        A=M.cells[i]
        intervals,geom=bounds(A,B,e[i],e[j],es,q,129)
        for lo,hi in intervals:
            l,h=geom(lo)[0],geom(hi)[0]
            l,h=max(l,B.lo),min(h,B.hi)
            if h>l:
                result.append((l,h,i,es,fs,width))
    return result


def cell_objects(C,j,loge,logf,oversample,poles,n,fix_gauge,mesh_breaks=(),min_half_width=2e-15):
    M=C.M
    B=M.cells[j]
    q=C.q
    es,_,_,fs,ms=poles
    x0=.05*np.exp(getattr(M,'eta',0.))
    coeff=C.sw*(G*KAPPA*x0*((1.06*x0)**2-ms*ms))/(q*es)
    parts=[]
    for k in range(len(n)):
        parts.extend(receiver_packets(C,loge,j,es[k],fs[k],q[k],coeff[k]))
    cut=np.unique(np.r_[B.lo,B.hi,mesh_breaks,[v for u in parts for v in u[:2]]])
    fine=np.unique(np.concatenate([np.linspace(a,b,oversample+1) for a,b in zip(cut[:-1],cut[1:])]))
    dx=np.diff(fine)/2
    mid=(fine[:-1]+fine[1:])/2
    # Do not create bins narrower than a handful of stored spacings.
    good=dx>min_half_width
    dx=dx[good]
    mid=mid[good]
    def sample(order):
        xi,w=roots_legendre(order)
        pp=mid[:,None]+dx[:,None]*xi
        weight=dx[:,None]*w
        storage=np.zeros_like(pp)
        for lo,hi,i,ee,ss,sc in parts:
            A=M.cells[i]
            a0,a1=loge[i]
            sel=(pp>=lo)&(pp<=hi)
            if not np.any(sel):
                continue
            v=pp[sel]
            kd=loge[j,0]+loge[j,1]*B.basis(v)
            wp=kd+ee
            energy=A.ec+A.de*(wp-a0)/a1
            p=A.pfrombare(energy)
            zp=A.basis(p)
            if np.any(p<A.lo-2e-12) or np.any(p>A.hi+2e-12):
                raise ValueError('Parent remap left an admissible aperture interval.')
            parent=expit(logf[i,0]+logf[i,1]*zp)
            vp=(a1/A.de)*p/energy
            storage[sel]+=sc*(parent+ss)*p*v/(wp*kd*vp)
        density=2*pp*pp/np.pi**2
        u=B.basis(pp)
        components=[density,density*u]
        if fix_gauge:
            # 1/E lies in the span of 1, u, and u^2/(1+de/ec*u).
            components.append(density*u*u/(1+(B.de/B.ec)*u))
        A=np.asarray(components)
        return weight,A,storage,pp
    sample.radial_bounds=(mid-dx,mid+dx)
    w,A,k,pp=sample(12)
    Aeq=np.einsum('jnm,nm->jn',A,w)
    objective=np.sum(w*k,axis=-1)
    den=Aeq[0].sum()
    scale=max(objective.sum(),1.e-9)
    return sample,Aeq,objective,den,scale,len(parts)


def extrema(m,objects):
    sample,Aeq,W,total,scale,count=objects
    A=Aeq/total
    c=W/scale
    if count==0:
        return dict(piecewise_min=0.,piecewise_max=0.,dual_min=0.,dual_max=0.,
                    max_abs_moment_residual=0.,max_dual_order_delta=0.)
    # Construct a strictly interior feasible anchor with an affine bin occupation.
    ansatz=(Aeq/Aeq[0]).T
    lift=ansatz@np.linalg.inv(A@ansatz)
    anchor=lift@(m/total)
    if not np.all((anchor>0.)&(anchor<1.)):
        raise ValueError('Moment set not reached by the interior feasibility anchor.')
    values=[]
    for direction in (1.,-1.):
        res=linprog(direction*c,A_eq=A,b_eq=m/total,bounds=(0.,1.),method='highs',
                    options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
        if not res.success or np.any(res.x<0.) or np.any(res.x>1.):
            raise ValueError(f'Inadmissible finite LP: {res.message}')
        for slack in (0.,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2):
            primal_vector=(1-slack)*res.x+slack*anchor
            primal_vector+=lift@((m/total)-A@primal_vector)
            if np.all((primal_vector>=0.)&(primal_vector<=1.)):
                break
        else:
            raise ValueError('Could not certify a positive fitted-bin witness after bounded repair.')
        defect=np.max(abs(Aeq@primal_vector-m))
        if defect>5e-13*total:
            raise ValueError('The finite positive fit did not preserve its moments.')
        lam=res.eqlin.marginals*scale/total
        dual=[]
        for N in (32,64):
            w,A0,k,_=sample(N)
            residual=direction*k-np.einsum('j,jab->ab',lam,A0)
            dual.append(float(lam@m+np.sum(w*np.minimum(residual,0.))))
        primal=float(W@primal_vector)
        if direction*primal<dual[-1]-5e-12*scale:
            raise ValueError('Dual does not enclose the discrete positive witness.')
        values.append((primal,direction*dual[-1],float(defect),abs(dual[-1]-dual[0]),slack))
    return dict(piecewise_min=values[0][0],piecewise_max=values[1][0],
                dual_min=values[0][1],dual_max=values[1][1],
                max_abs_moment_residual=max(v[2] for v in values),
                max_dual_order_delta=max(v[3] for v in values),
                max_witness_slack=max(v[4] for v in values))


def run(archive,oversample,fix_gauge=False):
    C,n,(y0,_,y1),_,depth=initial(archive)
    M=C.M
    e,f=M.weights(y0,True)[:2]
    poles=C.scalar(n)
    a=np.asarray(M.edge_sets[2])
    sel=np.flatnonzero((a[:-1]>=1.4138888888888888-2e-14)&(a[1:]<=1.4666666666666666+2e-14))
    if len(sel)!=2*depth:
        raise ValueError('Wrong physical window.')
    parts=[]
    fitted=[]
    for j in M.r[2].start+sel:
        O=cell_objects(C,int(j),e,f,oversample,poles,n,fix_gauge)
        B=M.cells[j]
        w,A,coeff,points=O[0](24)
        fits=[B.fit(q[j])[0] for q in (y0,y1)]
        aa,bb=(expit(z[0]+z[1]*B.basis(points)) for z in fits)
        fitted.append(float(np.sum(w*coeff*(aa-bb))))
        m0,m1=y0[j],y1[j]
        if fix_gauge:
            m0=np.r_[m0,np.sum(w*A[2]*aa)]
            m1=np.r_[m1,np.sum(w*A[2]*bb)]
        parts.append(dict(cell=int(j),microbins=len(O[2]),packets=O[-1],
                          native=extrema(m0,O),archive_local=extrema(m1,O)))
    def results(tag):
        def col(k,t):
            return sum(v[t][k] for v in parts)
        return [col(tag+'_min','native')-col(tag+'_max','archive_local'),
                col(tag+'_max','native')-col(tag+'_min','archive_local')]
    return dict(scope='A fixed-moment bounded-occupation envelope for one held initial aperture; numerical LP and dual quadrature, not an evolved result.',
                hard_depth=depth,cold_modes=len(n),subcells_per_packet_interval=oversample,
                fixed_gauge_moment=fix_gauge,independent_fitted_rate=sum(fitted),
                piecewise_witness_range=results('piecewise'),integrated_dual_range=results('dual'),
                cells=parts)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives',nargs='+')
    p.add_argument('--oversample',type=int,default=4)
    p.add_argument('--fix-gauge',action='store_true')
    a=p.parse_args()
    if a.oversample<1:
        p.error('Positive subcells required.')
    print(json.dumps([run(f,a.oversample,a.fix_gauge) for f in a.archives],indent=2,allow_nan=False))
