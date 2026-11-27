"""Independent plasma-frame exchange-momentum quadrature for a two-flavor heat subset."""
import argparse
import json

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import expit, roots_legendre

from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_frozen_screened_elastic_amplitude_reference import tabulated_polarization
from early_pole_live_joint_hard_pair_packets_reference import gauss_each


def profile(grid, species, form):
    d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
    e=np.asarray(d[f'endpoint_{species}_momentum_edges'])
    f=np.asarray(d[f'endpoint_{species}_occupations'])
    mid=.5*(e[1:]+e[:-1])
    if form=='steps':
        vals=np.r_[0.,f,0.]
        return e,lambda q:vals[np.searchsorted(e,q)]
    if form=='pchip':
        h=PchipInterpolator(np.r_[e[0],mid,e[-1]],np.r_[f[0],f,f[-1]],extrapolate=False)
        return e,lambda q:np.nan_to_num(h(q))
    v=e**3/3
    cdf=PchipInterpolator(v,np.r_[0.,np.cumsum(np.diff(v)*f)],extrapolate=False)
    h=cdf.derivative()
    a,b,_,_=cdf.c
    u=np.clip(np.divide(-b,3*a,out=np.zeros_like(a),where=a!=0.),0.,np.diff(v))
    test=h(np.r_[v,v[:-1]+u])
    assert np.min(test)>-2e-12 and np.max(test)<1.
    return e,lambda q:np.maximum(np.nan_to_num(h(q**3/3)),0.)


def radial(e,Q,order,layer,bounds):
    lo,hi=bounds
    extra=np.array([])
    if layer:
        extra=np.r_[e-Q,e+Q,e-.25*Q,e+.25*Q]
    b=np.unique(np.r_[lo,hi,e[(e>lo)&(e<hi)],extra[(extra>lo)&(extra<hi)]])
    return gauss_each(b,order)


def traces(a,b,q,z,min_s,lower_tau,v):
    # Massless external fermions, fixed Q and omega = Q z.
    t=-q*q*(1-z*z)
    w=q*z
    ap=a*z+t/(2*q)
    bp=b*z-t/(2*q)
    at=np.maximum(a*a-ap*ap,0.)
    bt=np.maximum(b*b-bp*bp,0.)
    ab=np.sqrt(at*bt)
    den=np.maximum(2*ab,1e-300)
    s0=2*(a*b-ap*bp)
    cu=np.clip((s0-min_s)/den,-1.,1.)
    cl=np.clip((s0+t/lower_tau)/den,-1.,1.) if lower_tau else np.full_like(cu,-1.)
    good=cl<cu
    L=np.arccos(np.where(good,cu,-1.))
    R=np.arccos(np.where(good,cl,-1.))
    I0=(R-L)/np.pi
    I1=(np.sin(R)-np.sin(L))/np.pi
    I2=.5*I0+(np.sin(2*R)-np.sin(2*L))/(4*np.pi)
    L00=(4*a*(a+w)+t)*(4*b*(b-w)+t)*I0
    mix=4*(2*a+w)*(2*b-w)*ab*I1
    Tr=16*ab*ab*I2+(-4*t*(at+bt)+2*t*t)*I0
    D0=1/(q*q)
    DT=-1/t
    vacuum=L00*D0**2+Tr*DT**2-2*mix*D0*DT
    a0,a1,c0,c1=np.moveaxis(v(abs(z)),-1,0)
    D0=1/(q*q+a0+1j*c0)
    DT=1/(-t+a1-1j*c1)
    screened=L00*abs(D0)**2+Tr*abs(DT)**2-2*mix*np.real(D0*np.conj(DT))
    return vacuum,screened


def one(grid, form, level, lower_tau, layer, window, bounds=(.45,8.)):
    pl,ph=bounds
    if not (0.<pl<ph and 0.<window[0]<window[1]<=2*ph and 0.<=lower_tau<1.):
        raise ValueError('Require positive radial bounds, Q <= 2 p_max and 0 <= min_tau < 1')
    er,fr=profile(grid,'parent',form)
    ew,fw=profile(grid,'fermion',form)
    v,mD2=tabulated_polarization(grid)
    qs=np.array([5e-5,.001,.0025,.005,.01,.02,.04,.08,.16,.32,.64,1.28,2.56,5.12,10.24,16.])
    qs=np.unique(np.r_[window[0],qs[(qs>window[0])&(qs<window[1])],window[1]])
    nr,nq,nz=[(1,3,4),(2,4,6),(3,5,8)][level]
    xx,ww=roots_legendre(nz)
    o=np.zeros((len(qs)-1,4))
    # Isotropic q-axis measure, two distinguishable Dirac spin/state sums,
    # and half the net forward-minus-reverse phase space.
    factor=.005**4/(np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18)/(32*np.pi**5)
    for iq,(L,H) in enumerate(zip(qs[:-1],qs[1:],strict=True)):
        lq,dq=gauss_each(np.log([L,H]),nq)
        for q,wq in zip(np.exp(lq),dq,strict=True):
            pa,dpa=radial(er,q,nr,layer,bounds)
            pb,dpb=radial(ew,q,nr,layer,bounds)
            ff,gg=fr(pa),fw(pb)
            for p,dp,f0 in zip(pa,dpa,ff,strict=True):
                lo=np.maximum(np.maximum((abs(p-q)-p)/q,(pl-p)/q),
                              np.maximum(-1.,(pb-ph)/q))
                hi=np.minimum(np.minimum(1.,(ph-p)/q),
                              np.minimum((pb-abs(pb-q))/q,(pb-pl)/q))
                possible=hi>lo
                # Outgoing FD breakpoints: only edges a distance at most Q from each leg.
                ez=(er[(er>p-q)&(er<p+q)]-p)/q
                # Most entries of the daughter mass-edge matrix lie outside the allowed interval.
                daughter=(pb[:,None]-ew[None,:])/q
                candidates=np.r_[np.ones((1,len(pb)))*ez[:,None],daughter.T].T
                width=(q*q/mD2)
                extra=np.r_[0.,np.array([.15,.6,2.,8.])*width,-np.array([.15,.6,2.,8.])*width]
                c0=np.column_stack([candidates,np.broadcast_to(extra,(len(pb),len(extra)))])
                zb=np.sort(np.column_stack([lo,hi,np.clip(c0,lo[:,None],hi[:,None])]),axis=1)
                half=np.diff(zb,axis=1)/2
                j,k=np.nonzero((half>2e-14)&possible[:,None])
                if not j.size:
                    continue
                mid=(zb[j,k]+zb[j,k+1])/2
                z=mid[:,None]+half[j,k,None]*xx
                p1=p+q*z
                p2=pb[j,None]-q*z
                Va,Sc=traces(p,pb[j,None],q,z,1.,lower_tau,v)
                dW=factor*dp*dpb[j,None]*q*q*wq*half[j,k,None]*ww
                fo,go=fr(p1),fw(p2)
                gi=gg[j,None]
                net=f0*gi*(1-fo)*(1-go)-fo*go*(1-f0)*(1-gi)
                dE=-q*z
                o[iq,:2]+=np.r_[np.sum(dW*dE*net*Va),np.sum(dW*dE*net*Sc)]
                # Independent equal-temperature equilibrium and hotter-parent witness.
                ea,eb=expit(.1-p),expit(-.4-pb[j,None])
                ea1,eb1=expit(.1-p1),expit(-.4-p2)
                eq=ea*eb*(1-ea1)*(1-eb1)-ea1*eb1*(1-ea)*(1-eb)
                eb,eb1=expit(-.4-pb[j,None]/.65),expit(-.4-p2/.65)
                hot=ea*eb*(1-ea1)*(1-eb1)-ea1*eb1*(1-ea)*(1-eb)
                o[iq,2:]+=np.r_[np.sum(dW*dE*eq*Sc),np.sum(dW*dE*hot*Sc)]
    assert np.max(abs(o[:,2]))<5e-9 and np.min(o[:,3])>-1e-12 and sum(o[:,3])>0.
    return dict(grid=grid,occupation_model=form,order=[nr,nq,nz],q_edges=qs.tolist(),
                min_tau=lower_tau,radial_bounds=bounds,radial_edge_layers=layer,mD2=mD2,normalized_heat=o.tolist(),
                columns=['archive_vacuum','archive_screened','common_T_screened','parent_hot_screened'])

def azimuth_fixture():
    from early_pole_frozen_screened_elastic_amplitude_reference import matrix_trace, htl_massless
    rng=np.random.default_rng(41)
    pol=htl_massless(.005**2/3)
    u,w=roots_legendre(12)
    worst=np.zeros(2); num=0
    for p,b,q,z in np.c_[rng.uniform(.6,4.,(160,2)),rng.uniform(.025,.55,160),rng.uniform(-.7,.7,160)]:
        for cut in (0.,.15):
            t=-q*q*(1-z*z)
            ap=p*z+t/(2*q); bp=b*z-t/(2*q)
            en=2*np.sqrt((p*p-ap*ap)*(b*b-bp*bp))
            s0=2*(p*b-ap*bp)
            L=np.arccos(np.clip((s0-1)/en,-1,1))
            R=np.arccos(np.clip((s0+t/cut)/en,-1,1)) if cut else np.pi
            if R-L<1e-9: continue
            psi=(R+L)/2+(R-L)*u/2
            S=s0-en*np.cos(psi)
            exact=np.array([w@matrix_trace(S,*np.broadcast_arrays(t+np.zeros_like(S),p,b,q*z),(0.,0.),a)*(R-L)/(2*np.pi) for a in (None,pol)])
            test=np.asarray(traces(p,b,q,z,1.,cut,pol))
            worst=np.maximum(worst,abs(test/exact-1)); num+=1
    assert worst.max()<1e-9
    return dict(cases=num,max_relative_error=worst.tolist())

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grid',choices=('fine','third'),default='fine')
    p.add_argument('--form',choices=('steps','pchip','cell-conservative'),default='pchip')
    p.add_argument('--level',type=int,choices=range(3),default=0)
    p.add_argument('--lower-tau',type=float,default=0.)
    p.add_argument('--layer',action='store_true')
    p.add_argument('--window',type=float,nargs=2,default=[5e-5,.08])
    p.add_argument("--bounds",type=float,nargs=2,default=[.45,8.])
    p.add_argument("--check-azimuth",action="store_true")
    args=p.parse_args()
    result=azimuth_fixture() if args.check_azimuth else one(args.grid,args.form,args.level,args.lower_tau,args.layer,args.window,args.bounds)
    print(json.dumps(result,indent=2,allow_nan=False))
