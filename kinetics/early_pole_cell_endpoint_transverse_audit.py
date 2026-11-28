"""Transverse frozen pole on fixed-moment charge endpoint populations."""
import argparse
import json
from pathlib import Path
from math import fsum

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit,logit,roots_legendre

from early_pole_cell_endpoint_collision_pilot import CellChart,split
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import G,Species
from finite_mass_qed_timelike_pole_angular_reference import angle


def charge_field(p,N,m,beta,order=36):
    # G^2 dN times the by-parts derivative divided by p^2. Written
    # with 2/E instead of 2*v/p so it is finite at a zero endpoint.
    E=np.hypot(m,p)
    v=p/E
    dv=(m*m)/E**3
    u,w=roots_legendre(order)
    h=(1-u*u)*w/4
    Y=1-beta*v[:,None]*u
    H=(1/Y)@h
    D=(1/Y**2)@(h*u)
    DD=(2/Y**3)@(h*u*u)
    z=(2/E+dv)*H + v*beta*dv*D
    s=(2/E+dv)*v*D + v*dv*(D+beta*v*DD)
    return G*G*np.array([fsum(N*z),fsum(N*s)])


def root(field,k):
    gap=np.sqrt(field(0.)[0])
    if not gap>0.:
        raise ArithmeticError('No positive plasma gap')
    lo=np.hypot(gap,k)
    hi=np.sqrt(k*k+3*gap*gap)
    def z(w):
        return (w-k)*(w+k)-field(k/w)[0]
    w=brentq(z,np.nextafter(lo,0.),hi,xtol=5e-17)
    beta=k/w
    derivative=field(beta)[1]
    h=2e-6
    check=(field(beta+h)[0]-field(beta-h)[0])/(2*h)
    if not derivative>0 or abs(check/derivative-1)>2e-7:
        raise ArithmeticError('Residue derivative failed the independent centered check')
    D=2*w+k/w**2*derivative
    if D<=0 or abs(z(w))>2e-12*gap**2:
        raise ArithmeticError('Transverse root or positive residue failed')
    return np.array([gap,w,1/D])



def continuum_moments(edges,m,count,order):
    # Enforce the physical number and first energy of each charge cell;
    # use a Fermi-logit affine ansatz rather than placing probability at its
    # two nodes. Coordinates u are normalized separately in every cell.
    e=np.hypot(m,edges)
    target=np.c_[count.reshape(-1,2).sum(1),np.diff(count.reshape(-1,2),axis=1)[:,0]]
    if np.any(target[:,0]<=0):
        raise ArithmeticError('Moment reconstruction here requires occupied charged cells')
    l,r=edges[:-1,None],edges[1:,None]
    q,w=roots_legendre(order)
    p=(l+r)/2+(r-l)*q/2
    u=2*(np.hypot(m,p)-e[:-1,None])/np.diff(e)[:,None]-1
    W=4*p*p*(r-l)*w/(4*np.pi**2)
    N=target[:,0]
    alpha=logit(N/W.sum(axis=1))
    beta=np.zeros_like(alpha)
    for j in range(70):
        f=expit(alpha[:,None]+beta[:,None]*u)
        F=W*f
        q0=F.sum(axis=1)/N-1
        q1=np.sum(F*u,axis=1)/N-target[:,1]/N
        defect=max(np.max(abs(q0)),np.max(abs(q1)))
        if defect<3e-13:
            return p.ravel(),F.ravel(),float(defect),float(np.max(abs(beta))),alpha,beta
        V=F*(1-f)/N[:,None]
        a=np.sum(V,axis=1)
        b=np.sum(V*u,axis=1)
        c=np.sum(V*u*u,axis=1)
        det=a*c-b*b
        if min(det)<=1e-50:
            raise ArithmeticError('Unresolved or infeasible continuum two-moment fit')
        da=(q0*c-q1*b)/det
        db=(q1*a-q0*b)/det
        scale=np.minimum(1.,5/np.maximum(5,np.maximum(abs(da),abs(db))))
        alpha-=scale*da
        beta-=scale*db
    raise ArithmeticError(('No two-moment continuum fit for these cells',defect))

def piecewise_logit(m,edges,alpha,beta,count,n):
    if n==1:
        volume=4*np.diff(edges**3)/(6*np.pi**2)
        return Species(m,edges,count.reshape(-1,2).sum(1)/volume)
    a=split(edges,n)
    x,w=roots_legendre(28)
    l,r=a[:-1,None],a[1:,None]
    p=(l+r)/2+(r-l)*x/2
    E=np.hypot(m,p)
    h=np.hypot(m,edges)
    u=2*(E-np.repeat(h[:-1],n)[:,None])/np.repeat(np.diff(h),n)[:,None]-1
    f=expit(np.repeat(alpha,n)[:,None]+np.repeat(beta,n)[:,None]*u)
    V=4*(a[1:]**3-a[:-1]**3)/(6*np.pi**2)
    num=np.sum((r-l)/(4*np.pi**2)*4*w*p*p*f,axis=1)
    if max(abs(num.reshape(-1,n).sum(1)/count.reshape(-1,2).sum(1)-1))>5e-11:
        raise ArithmeticError('Step discretization failed its parent continuum cell numbers')
    return Species(m,a,num/V)


def model_callback(m,edges,alpha,beta):
    energy=np.hypot(m,edges)
    def at(E):
        i=np.searchsorted(energy,E,side='right')-1
        valid=(i>=0)&(i<len(alpha))
        z=np.zeros_like(E)
        j=i[valid]
        u=2*(E[valid]-energy[j])/(energy[j+1]-energy[j])-1
        z[valid]=expit(alpha[j]+beta[j]*u)
        return z
    return at


def state(d,x,y):
    ps=(d.qp,d.qd)
    p,f=(d.wp*y[d.iP:],d.wd*y[d.ns:d.iP])
    masses=(x,x*d.ratio)
    X=[(p0,n,m) for p0,n,m in zip(ps,(p,f),masses,strict=True)]
    def effective(order):
        return lambda beta:sum((charge_field(*a,beta,order) for a in X))
    v=root(effective(36),.0015)
    exact=root(effective(52),.0015)
    if np.max(abs(v/exact-1))>3e-13:
        raise ArithmeticError('Independent polar order did not resolve root and residue')
    # A constant step fit preserves number per native bin; compare this
    # alternate reconstruction to the raw two-node charge occupations.
    a=[]
    for mass,shell,data in zip(masses,(d.p_edges,d.d_edges),(p,f),strict=True):
        count=data.reshape(-1,2).sum(1)
        volume=4*np.diff(shell**3)/(6*np.pi**2)
        a.append(Species(mass,shell,count/volume))
    step=root(lambda t:sum(np.array(c.field(t)) for c in a),.0015)
    return v,step


def run(args):
    observations=[]
    for g in ('fine','third'):
        src=load(g)[0]
        x=src['x_end']
        for levels in (1,2,4,8):
            d=CellChart(src,SPECS[g],levels)
            v,s=state(d,x,d.initial)
            # Replay against the original discontinuous-edge product.
            a=[Species(m,src[f'endpoint_{p}_momentum_edges'],src[f'endpoint_{p}_occupations'])
               for m,p in ((x,'parent'),(x*src['mass_ratio'],'fermion'))]
            direct=root(lambda t:sum(np.array(c.field(t)) for c in a),.0015)
            if np.max(abs(direct/s-1))>4e-12:
                raise ArithmeticError('Frozen archive step moment was not recovered')
            m=src['checkpoints'][-1]['scalar_pole_over_T']
            aperture=[]
            if levels in (1,8):
                q=3.5e-6
                Es=np.hypot(m,q)
                controls=[]
                for poles in (v,s):
                    low=angle(a,.0015,poles[1],poles[2],Es,'T',3,(3,8),q)[0]
                    high=angle(a,.0015,poles[1],poles[2],Es,'T',6,(6,12),q)[0]
                    if max(abs(low/high-1))>2e-7:
                        raise ArithmeticError('Cold sensitivity test not resolved')
                    controls.append(high.tolist())
                aperture=np.divide(controls[0],controls[1]).tolist()
            observations.append(dict(grid=g,kind='initial_archive_step',x=x,subcells=levels,
                                     raw_gap_frequency_residue=v.tolist(),constant_reference=s.tolist(),
                                     relative_discrepancy=(v/s-1).tolist(),cold_A_B_D_factor_at_q_3p5e_6=aperture))
    for f in args.paths:
        j=json.loads(Path(f).read_text())
        g=j['grid']
        src=load(g)[0]
        d=CellChart(src,SPECS[g],j['subcells_per_archived_cell'],reference_x=j['reference_x'])
        x=j['x_end']
        y=np.array(j['endpoint_occupation'])
        if j['counts']!=[d.ns,d.nf,d.np] or y.shape!=(d.n,):
            raise ArithmeticError('Bad endpoint snapshot for transverse audit')
        v,s=state(d,x,y)
        m=d.state(x,y)[0]
        charge=[]
        for ed,m1,n in ((d.p_edges,x,d.wp*y[d.iP:]),(d.d_edges,x*d.ratio,d.wd*y[d.ns:d.iP])):
            d0=continuum_moments(ed,m1,n,24)
            d1=continuum_moments(ed,m1,n,36)
            charge.append((m1,d0,d1))
        fits=[]
        for ix in (1,2):
            fits.append(root(lambda z:sum(charge_field(*u[ix][:2],u[0],z,36) for u in charge),.0015))
        if np.max(abs(fits[0]/fits[1]-1))>4e-11:
            raise ArithmeticError('Continuum optical closure underresolved')
        probe=[]
        exact_probe=None
        if args.rate:
            for n in (2,4,8,16):
                targets=[]
                for (mass,_,fit),(ed,cnt) in zip(charge,
                    ((d.p_edges,d.wp*y[d.iP:]),(d.d_edges,d.wd*y[d.ns:d.iP])),strict=True):
                    targets.append(piecewise_logit(mass,ed,*fit[-2:],cnt,n))
                ES=np.hypot(m,3.5e-6)
                low=angle(targets,.0015,fits[1][1],fits[1][2],ES,'T',4,(4,8),3.5e-6)[0]
                high=angle(targets,.0015,fits[1][1],fits[1][2],ES,'T',6,(6,12),3.5e-6)[0]
                probe.append(dict(subintervals=n,raw_A_B_D=high.tolist(),
                                  angle_mesh_difference=float(max(abs(low/high-1)))))
            targets=[]
            callbacks=[]
            for (mass,_,fit),(ed,cnt) in zip(charge,
                ((d.p_edges,d.wp*y[d.iP:]),(d.d_edges,d.wd*y[d.ns:d.iP])),strict=True):
                targets.append(piecewise_logit(mass,ed,*fit[-2:],cnt,1))
                callbacks.append(model_callback(mass,ed,*fit[-2:]))
            ES=np.hypot(m,3.5e-6)
            def integration(n,a,s,f):
                return angle(targets,.0015,fits[1][1],fits[1][2],ES,'T',n,(a,s),3.5e-6,
                             occupation_callbacks=f)[0]
            def stepwise(ch):
                E=np.hypot(ch.m,ch.e)
                return lambda p: np.r_[0.,ch.f,0.][np.searchsorted(E,p,side='right')]
            c=integration(6,6,12,None)
            v0=integration(6,6,12,[stepwise(f) for f in targets])
            match=float(np.max(abs(c/v0-1)))
            r0=integration(8,10,18,callbacks)
            r1=integration(12,16,24,callbacks)
            refine=float(np.max(abs(r0/r1-1)))
            if match>4e-12 or refine>2e-5:
                raise ArithmeticError(('Smooth radial test not resolved',match,refine))
            exact_probe=dict(raw_A_B_D=r1.tolist(),mesh_discrepancy=refine,
                             legacy_step_callback_discrepancy=match,
                             last_piecewise_step_difference_by_coefficient=(np.array(probe[-1]['raw_A_B_D'])/r1-1).tolist())
        observations.append(dict(grid=g,kind='evolved_endpoint_archive',x=x,subcells=j['subcells_per_archived_cell'],
                                 scalar_mass=m,raw_gap_frequency_residue=v.tolist(),constant_reference=s.tolist(),
                                 relative_discrepancy=(v/s-1).tolist(),
                                 two_moment_gap_frequency_residue=fits[1].tolist(),
                                 relative_two_moment_to_constant_difference=(fits[1]/s-1).tolist(),
                                 max_moment_fit_defect=max(u[2][2] for u in charge),
                                 max_logit_energy_slope=max(u[2][3] for u in charge),
                                 continuum_rate_discretization_probe=probe,
                                 direct_continuum_aperture=exact_probe))

    return observations


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('paths',nargs='*')
    p.add_argument('--rate',action='store_true')
    print(json.dumps(run(p.parse_args()),indent=2,allow_nan=False))
