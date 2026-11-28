"""Nested daughter subdivisions holding the input histogram exactly fixed."""
import json
from math import fsum

import numpy as np
from scipy.special import roots_legendre

from early_pole_double_charge_stencil_pilot import chart,exact_step_on_common_span
from early_pole_endpoint_stencil_remap_pilot import LABELS,SPECS,archives,load,parts,remap,MassiveFermion
from early_pole_live_parent_energy_reference import state_moments


def subdivide_daughter(d,level):
    a=d.d.fed
    n=2**level
    b=np.concatenate([np.linspace(l,h,n+1)[:-1] for l,h in zip(a[:-1],a[1:],strict=True)]+[a[-1:]])
    if not np.array_equal(b[::n],a):
        raise ArithmeticError('Daughter refinement is not nested')
    d.d.fed=b
    d.d.ef=np.hypot(b,d.d.b)
    d.d.wf=4*np.diff(b**3)/(6*np.pi**2)
    d.nf=len(d.d.wf)
    d.iP=d.ns+d.nf
    d.n=d.iP+d.np
    d.w=np.r_[d.d.w,d.d.wf,d.pw]
    z,w=roots_legendre(16)
    step=np.diff(b)[:,None]
    d.daughter_q=b[:-1,None]+step*(z+1)/2
    d.daughter_w=4*step*w*d.daughter_q**2/(4*np.pi**2)
    return b


def affine_arms(ch,d,y,x):
    Es,Ed,Ep,i,j,k,a,b,W=ch
    E_d=(1-a)*Ed[j-1]+a*Ed[j]
    E_p=(1-b)*Ep[k-1]+b*Ep[k]
    fs,fd,fp=np.split(y,(d.ns,d.iP))
    Ip=np.searchsorted(np.hypot(d.pe,x),E_p)-1
    Id=np.searchsorted(np.hypot(d.d.fed,d.daughter_mass(x)),E_d)-1
    pg=d.parent_products(fp[k-1],fp[k],1-b,b,False)
    dg=d.parent_products(fd[j-1],fd[j],1-a,a,False)
    ps=fp[Ip],1-fp[Ip]
    ds=fd[Id],1-fd[Id]
    stats=[]
    for P,D in ((ps,ds),(pg,ds),(ps,dg),(pg,dg)):
        C=fsum(W*P[0]*(1+fs[i])*D[1]),fsum(W*P[1]*fs[i]*D[0])
        stats.append(C)
    return np.asarray(stats)


def run():
    out=[]
    for grid,shape in SPECS.items():
        archive=load(grid)[0]
        old=archives(archive)
        oldf=[archive['endpoint_'+k+'_occupations'] for k in LABELS]
        x=archive['x_end']
        levels=[]
        for l in (0,1,2,3):
            d=MassiveFermion(*shape,archive['mass_ratio'],archive['elastic_rate'],True)
            e=subdivide_daughter(d,l)
            parts0=parts(d)
            f=[remap(a,v,b)[0] for a,v,b in zip(old,oldf,parts0,strict=True)]
            np.testing.assert_allclose(f[1],np.repeat(oldf[1],2**l),atol=1e-29,rtol=4e-16)
            y=np.concatenate(f)
            _,dr,pr,m=state_moments(d,x,y)
            n=fsum(d.d.w*f[0])
            scale=np.sqrt(d.radiation/(d.radiation+dr+pr))/n
            A=[]
            for order in (4,8):
                ch,geo=chart(d,x,m,order,histogram_cuts=True)
                A.append(affine_arms(ch,d,y,x)*scale)
            exact=exact_step_on_common_span(d,x,m,y)*scale
            if np.max(abs(A[-1][0]/exact-1))>8e-13:
                raise ArithmeticError('Discontinuous histogram was not resolved by the shared cuts')
            defect=float(np.max(abs(A[-1]/A[-2]-1.)))
            if defect>5e-6:
                raise ArithmeticError(('Shared interpolation quadrature did not converge',grid,l,defect))
            clocks=A[-1][:,0]-A[-1][:,1]
            levels.append(dict(level=l,daughter_cells=len(e)-1,gross_forward_and_inverse=A[-1].tolist(),
                               net_clocks=clocks.tolist(),last_quadrature_difference=defect,**geo))
        out.append(dict(grid=grid,preconditioning_arm_order=['both_exact_steps','parent_smoothed',
                       'daughter_smoothed','both_smoothed'],levels=levels))
    return out


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
