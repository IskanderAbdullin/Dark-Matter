"""Re-evaluate retained cold coefficients on the same-restart D8 clock pair."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,model_callback,piecewise_logit,root
from early_pole_endpoint_peer_charge_shape_transport import unpack
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from finite_mass_qed_timelike_pole_angular_reference import angle

Q=(3.5e-6,1.5e-5)
K=.0015


def setup(d,y,params):
    x=.05
    raw=[]; ch=[]; cb=[]
    for mass,edges,n in ((x,d.p_edges,d.wp*y[d.iP:]),(x*d.ratio,d.d_edges,d.wd*y[d.ns:d.iP])):
        p,w,_,_,a,b=continuum_moments(edges,mass,n,36)
        raw.append((p,w,mass))
        cb.append(model_callback(mass,edges,a,b))
        ch.append(piecewise_logit(mass,edges,a,b,n,1))
    _,omega,R=root(lambda beta:sum(charge_field(*a,beta,36) for a in raw),K)
    m=d.state(x,y)[0]
    if np.max(abs(np.array([m,omega,R])/np.array(params)-1))>2e-13:
        raise ArithmeticError('Reconstructed optical clock parameters changed')
    return ch,cb,m,omega,R


def point(d,u,opts,q):
    ch,cb,m,w,R=opts
    E=np.hypot(q,m)
    f=u[np.searchsorted(d.ds.ed,q,side='right')-1]
    out=[]
    for i,j,l in ((4,10,18),(8,16,24)):
        a,b,s=angle(ch,K,w,R,E,'T',i,(j,l),q,occupation_callbacks=cb,mirror_azimuth=True)[0]
        if min(a,b,s)<=0 or abs(a-b-s)>4e-13*max(a,b,s):
            raise ArithmeticError('The cold clock test requires nonzero affine components')
        z=q*q/E*np.array([(1+f)*a,f*b,a+f*s])
        out.append(z)
    defect=float(np.max(abs(out[0]/out[1]-1)))
    if defect>2e-6:
        raise ArithmeticError('Clock point coefficient quadrature not resolved')
    return out[1],defect


def run(low,high,clock,q_points=Q):
    a,sha=low
    b,sig=high
    g=b['grid']
    d,z,_=unpack(b,g,8)
    if a['steps']!=1024 or any(a[k]!=b[k] for k in ('grid','counts','reference_x','x_end','subcells_per_archived_cell')):
        raise ArithmeticError('Clock states must be saved on exactly the same endpoint chart')
    y=np.asarray(a['endpoint_occupation'])
    if y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1:
        raise ArithmeticError('Invalid recovered occupation array')
    if (clock['grid'],clock['steps'],clock['start'],clock['end'],clock['k'])!=(g,[1024,2048],.04,.05,K):
        raise ArithmeticError('Only a matched last-segment clock is declared')
    vectors=(clock['optical_params_coarse'],clock['optical_params_fine'])
    opts=[setup(d,u,p) for u,p in zip((y,z),vectors,strict=True)]
    out=[]
    for q in q_points:
        (u,e),(v,f)=[point(d,h,o,q) for h,o in zip((y,z),opts,strict=True)]
        neutral=[w[1]/w[2] for w in (u,v)]
        out.append(dict(q=q,coarse=u.tolist(),fine=v.tolist(),coarse_fine_relative_difference=(u/v-1).tolist(),
                        neutrality_coarse_fine_relative=neutral[0]/neutral[1]-1,max_joint_mesh_defect=max(e,f)))
    return dict(grid=g,x=.05,k=K,clock_steps=[1024,2048],sha256_uncompressed_states=[sha,sig],rows=out,
                scope=('At the same fixed k two q points, independently recovered scalar gap and photon fits; same .04 restart only. No whole-history or band error bound.' if tuple(q_points)==Q else 'At fixed k and the declared q point set, independently recovered scalar gap and photon fits; same .04 restart only. No whole-history or band error bound.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('coarse','fine','clock'):
        p.add_argument(n,type=Path)
    p.add_argument('--q',type=float,choices=Q,action='append')
    a=p.parse_args()
    points=tuple(a.q) if a.q else Q
    if len(set(points))!=len(points):
        p.error('q points must be distinct')
    print(json.dumps(run(read(a.coarse),read(a.fine),json.loads(a.clock.read_text()),points),indent=2,allow_nan=False))
