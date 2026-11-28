"""Soft-cut transport from the frozen retarded and statistical responses."""
import argparse
import json
import numpy as np
from scipy.special import roots_legendre
from early_pole_frozen_dynamic_vlasov_reference import Snapshot,load,THERMAL,real_kernels
from early_pole_matched_wilson_geometry_reference import inputs
from early_pole_eikonal_soft_kick_sum_rules_reference import log_nodes
from early_pole_kick_characteristic_reference import sphere


def angular(n):
    a,w=roots_legendre(n)
    u,wu=log_nodes(2e-14,.25,n,55)
    edges=np.linspace(.25,.985,32)
    h=np.diff(edges)[:,None]/2
    nos=(edges[1:]+edges[:-1])[:,None]/2+h*a
    tip,tw=log_nodes(1e-10,.015,n,25)
    return np.r_[u,nos.ravel(),1-tip],np.r_[wu,(h*w).ravel(),tw]


def direct(U,uw,response,noise,g,extent,n,t):
    d=inputs(g)
    x=d['vP']*U
    k,kw=log_nodes(2e-5,extent,n,56)
    k=k[:,None]
    U=U[None,:]
    lr,tr,li,ti=(response/THERMAL).T
    nL,nT=((noise/THERMAL)*(2/np.pi)).T
    nL,nT=nL[None,:],nT[None,:]
    zT=d['vP']**2*(1-U*U)
    posL=(k*k+lr)**2+li*li
    posT=(k*k*(1-x*x)+tr)**2+ti*ti
    gauge=2*k*kw[:,None]*uw*(nL/posL+nT*zT/posT)
    p=k*U
    y=t*d['r']
    return np.array([np.sum(gauge*p*p),np.sum(gauge*k*k*(1-U*U)),
                     np.sum(gauge*(1-np.sinc(y*p/np.pi)))])


def run(n):
    out=[]
    g=.005
    u,uw=angular(n)
    d=inputs(g)
    s=d['vP']*u
    L,T,_,_=real_kernels(s,np.ones_like(s))
    thermal_R=THERMAL*np.column_stack([L,T,np.pi*s/2,np.pi*s*(1-s*s)/4])
    thermal_N=np.column_stack([thermal_R[:,2]/s,thermal_R[:,3]/s])
    spec=[('two_thermal',thermal_R,thermal_N,1.)]
    for grid in ('fine','third'):
        src=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
        for z in (.06,.18,.32):
            S=Snapshot(src,18,z)
            R,N=S.evaluate(s)
            assert min(R[:,2])>0. and min(R[:,3])>0. and min(N.ravel())>=0.
            mm=S.retarded(0.)[0]/THERMAL
            spec.append((grid+str(z),R,N,np.sqrt(mm)))
            if z==.18:
                zz=1e-5
                r,noise=S.evaluate(zz)
                tstar=noise[0]*zz/r[2]
                spec.append((grid+'_massless_matched_mD_and_noise',
                             thermal_R*mm,thermal_N*(mm*tstar),np.sqrt(mm)))
    for c in (8.,20.):
        hot=direct(u,uw,thermal_R,thermal_N,g,c,n,1.)
        old,pz,pp=sphere(c,n,d['vP'])
        independent=np.array([np.sum(old*pz*pz),np.sum(old*pp*pp),
                              np.sum(old*(1-np.sinc(d['r']*pz/np.pi)))])
        assert max(abs(independent/hot-1))<3e-5
        for name,R,N,m in spec:
            for id,cut in [('same_T',c),('same_mD',c*m)]:
                a=direct(u,uw,R,N,g,cut,n,1.)
                out.append(dict(case=name,cutoff_rule=id,cutoff_label=c,cutoff_mD0=cut,
                                mass_over_two_thermal=m,scattering_relative_to_two_thermal=(a/hot).tolist(),
                                normalized_outputs=a.tolist()))
    return dict(scope='Frozen collisionless isotropic dynamical polarization with statistically '
                      'independent fermion occupations; spacelike soft subset only, not a complete rate.',
                columns=['radial_second_moment','transverse_second_moment','kick_phase_shape'],rows=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nodes',type=int,default=8)
    print(json.dumps(run(p.parse_args().nodes),indent=2,allow_nan=False))
