"""Scan the exact paired bare polar kernel over archived cold32 momenta."""
import gzip,json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

D=np.longdouble


def audit(d,nodes,witness):
    M=D(d['x_end']);m=M*D(d['mass_ratio']);s=D(d['checkpoints'][-1]['scalar_pole_over_T'])
    A=(M*M-m*m+s*s)/2
    W=np.sqrt((A-M*s)*(A+M*s))
    Z=(M*M+m*m-s*s)/2
    C=(M+m)**2/(Z+M*m)
    V=-np.pi*C*W/(4*s*s)
    L=np.pi*Z/(2*W)
    xp=np.hypot(M,np.asarray(d['endpoint_parent_momentum_edges'],dtype=D))
    xb=np.hypot(m,np.asarray(d['endpoint_fermion_momentum_edges'],dtype=D))
    fp=np.asarray(d['endpoint_parent_occupations'])
    fb=np.asarray(d['endpoint_fermion_occupations'])
    center_p=A/s
    center_b=center_p-s
    pstar=W/s
    threshold=[]
    # The two aperture boundaries move with the scalar rapidity.
    # Use their inverse rapidity map, then check with a lab-frame root.
    for key,edges,mass,e0 in (('parent',xp,M,center_p),('daughter',xb,m,center_b)):
        y0=np.arccosh(e0/mass)
        for i,edge in enumerate(edges[1:-1],1):
            y=np.arccosh(edge/mass)
            t=s*np.sinh(abs(y-y0))
            if 1e-8<t<max(nodes):
                sign=1 if edge>e0 else -1
                fn=lambda q:float(np.hypot(D(q),s)/s*e0+D(q)/s*pstar*sign-edge)
                direct=brentq(fn,0.,max(nodes),xtol=2e-17)
                tol=float(abs(direct/t-1))
                if tol>2e-10:
                    raise ArithmeticError('Inverse-rapidity aperture edge did not match the lab root.')
                threshold.append(dict(q=float(t),cut=key,archive_edge=i,entry_side=sign,relative_root_defect=tol))
    threshold.sort(key=lambda t:t['q'])

    def sh(q):
        q=D(q)
        es=np.hypot(s,q)
        E=es*A/(s*s)
        half=q*W/(s*s)
        pe=(xp-E)/half
        be=(xb+es-E)/half
        u=sorted({-1.,0.,1.,*[v for v in (*pe,*be) if -1<v<1]})
        def occ(t):
            i=np.searchsorted(pe,t)-1
            j=np.searchsorted(be,t)-1
            if not 0<=i<len(fp) or not 0<=j<len(fb):
                raise ArithmeticError('Selected cold aperture exceeds archived occupations.')
            return D(fp[i])-D(fb[j])
        parts=[]
        for l,h in zip(u[:-1],u[1:],strict=True):
            p=occ((l+h)/2)
            sgn=-1 if h<=0 else 1
            parts.append((l,h,p,p*V*(h-l),p*L*sgn*np.log((abs(h) or 1.)/(abs(l) or 1.))))
        vc,lc=([fsum(row[i] for row in parts) for i in (3,4)])
        raw=vc+lc
        around=parts[u.index(0.)-1:u.index(0.)+1]
        coef=L*sum(o[2] for o in around)
        residuals=[]
        for eps in (1e-5,1e-8):
            active=[]
            for l,h,f,*_ in parts:
                lo,hi=(l,min(h,-eps)) if h<=0 else (max(l,eps),h)
                if hi<=lo:continue
                a,b=sorted((float(np.log(abs(lo))),float(np.log(abs(hi)))))
                active.append(quad(lambda x:float(f*(V*np.exp(x)+L)),a,b,epsabs=3e-11,epsrel=3e-12)[0])
            direct=fsum(active)
            correct=direct-coef*np.log(1/eps)+V*eps*sum(o[2] for o in around)
            residuals.append(float(abs(correct-raw)))
        if max(residuals)>6e-9:
            raise ArithmeticError('Independent radial cutoff replay failed.')
        k=len(u)-3
        if k!=sum(q>o['q'] for o in threshold) and q<=max(nodes):
            raise ArithmeticError('Cold aperture junction count did not follow the crossing ledger.')
        if k==0 and abs(raw-2*V*parts[0][2])>2e-11:
            raise ArithmeticError('The constant-shell aperture did not recover its Q-independent finite part.')
        return dict(q=float(q),radial_half_width=float(half),aperture_shell_junctions=k,
                    common_aperture_scaled_finite_part=raw,logarithmic_coefficient=float(coef),
                    normalized_constant_and_log_terms=[vc,lc],independent_cutoff_replay_defect=max(residuals))
    result=[sh(q) for q in [*nodes,1e-8,1e-7]]
    origin=sh(4.1e-6)
    a=next(x for x in witness['controls'] if x['tuned_step'] is None and abs(x['parent_energy']-
           np.hypot(float(s),4.1e-6)*float(A/(s*s)))<1e-12)
    moved=a['chosen_finite_part']-a['full_aperture_log_coefficient']*np.log(origin['radial_half_width']/3e-4)
    err=abs(origin['common_aperture_scaled_finite_part']-moved)
    if err>5e-8:
        raise ArithmeticError('The Q=4.1e-6 center did not reproduce the declared previous scale.')
    return dict(momentum_nodes=result,first_aperture_crossings=threshold,
                central_prior_finite_part_change_of_scale_defect=float(err))


def run():
    p=Path(__file__).parent/'results'
    with gzip.open(p/'early-pole-unsplit-hard8-cold32-native-2026-11-27.json.gz') as f:
        q=json.load(f)['cold_momentum_nodes_over_T']
    with open(p/'finite-mass-qed-shared-polar-delta-2026-11-27.json') as f:
        w=json.load(f)
    out=[]
    for n in ('fine','third'):
        with gzip.open(p/f'early-pole-spectrum-{n}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        out.append(dict(grid=n,**audit(d,q,next(z for z in w['rows'] if z['grid']==n))))
    return dict(scope='Exact paired bare positive polar kernel at vacuum aperture centers, normalization at its half-width. No physical rate, width, or scalar dynamics.',
                grid_name='retained cold32 momentum grid',rows=out)


if __name__=='__main__':print(json.dumps(run(),indent=2,allow_nan=False))
