"""True-event quantum rates from fixed-temperature cell-number reconstruction."""
import json
import numpy as np
from scipy.special import expit, roots_legendre
from massive_fv_fixed_bath_reference import bounds, kinematics
from thermal_cell_inversion_reference import shell_average, infer_shift


def rates(n,x,ratios=(.08,.18),dmu=.0,order=12):
    m,b=(a*x for a in ratios)
    K,D=kinematics(x,m,b)
    q=np.geomspace(2e-5,24.,n+1)
    e,r=np.hypot(q,m),np.hypot(q,b)
    # Recover local chemical potentials using cell-number measurements only.
    old_num, potential = [], []
    mus=(.9*m,-.015)
    for mass,mu,stat in zip((m,b),mus,("BE","FD")):
        data=np.array([shell_average(l,h,mass,1.,mu-mass,stat)
                       for l,h in zip(q[:-1],q[1:])])
        recovered=np.array([infer_shift(l,h,mass,1.,f,stat)+mass
                           for l,h,f in zip(q[:-1],q[1:],data)])
        old_num.append(data)
        potential.append(recovered)
    t,w=roots_legendre(order)
    L,H=np.broadcast_arrays(e[:-1,None]+np.zeros((1,n)),e[1:,None]+np.zeros((1,n)))
    rlo,rhi=r[None,:-1],r[None,1:]
    z=bounds(rlo,b,m,K,D)+bounds(rhi,b,m,K,D)
    pts=np.sort(np.stack([L,H,*[np.clip(v,L,H) for v in z]]),axis=0)
    S=np.zeros((2,7))
    for lo,hi in zip(pts[:-1],pts[1:]):
        dy=hi-lo
        for t0,w0 in zip(t,w):
            y=lo+dy*(t0+1)/2
            rl,rh=bounds(y,m,b,K,D)
            lo2=np.minimum(rhi,np.maximum(rlo,rl))
            hi2=np.maximum(lo2,np.minimum(rhi,rh))
            dr=hi2-lo2
            for t1,w1 in zip(t,w):
                z=lo2+dr*(t1+1)/2
                p=expit(sum(mus)+dmu-y-z)
                wt=dy*dr*w0*w1
                values=((1/np.expm1(y-potential[0][:,None]),expit(potential[1][None,:]-z)),
                        (old_num[0][:,None],old_num[1][None,:]))
                for j,(bi,fi) in enumerate(values):
                    g=wt*p*(1-fi)*(1+bi)
                    l=wt*(1-p)*fi*bi
                    net=g-l
                    S[j]+=np.array([g.sum(),l.sum(),abs(net).sum(),
                                    (net*y).sum(),(net*z).sum(),-(net*(y+z)).sum(),
                                    abs(g-l*np.exp(dmu)).sum()])
    pref=4/(2*np.pi**2)*x/(2*D/x**2)/4
    S*=pref
    def extract(z):
        gross=z[0]+z[1]
        return dict(relative_affinity_identity_residual=float(z[6]/gross),
                    abs_net_over_gross=float(z[2]/gross),gain=float(z[0]),loss=float(z[1]),
                    physical_energy_rates=[float(v) for v in z[3:6]],
                    energy_residual_over_gross=float(abs(z[3:6].sum())/gross))
    return dict(cells=n,M_over_T=x,mu_parent_shift=dmu,
                max_absolute_recovered_mu_error=max(float(np.max(abs(v-u)))
                                                    for v,u in zip(potential,mus)),
                reconstructed=extract(S[0]),constant_occupation=extract(S[1]))


def run():
    rows=[rates(n,x,dmu=d) for x in (.1,.005) for n in (24,48,96) for d in (0.,.2)]
    assert max(r["reconstructed"]["relative_affinity_identity_residual"] for r in rows)<3e-13
    assert max(r["reconstructed"]["energy_residual_over_gross"] for r in rows)<1e-13
    assert min(r["constant_occupation"]["abs_net_over_gross"] for r in rows if r["mu_parent_shift"]==0)>1e-6
    return dict(rows=rows)
if __name__ == "__main__":
    print(json.dumps(run(),indent=2,allow_nan=False))
