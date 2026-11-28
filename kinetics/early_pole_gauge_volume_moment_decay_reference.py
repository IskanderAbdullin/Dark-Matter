"""Retain event-resolved energy and momentum as separate finite-volume moments."""
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre,expit

from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_snapshot_gauge_ledger_reference import cell_J,photon
from early_pole_gauge_projection_momentum_reference import curve


def one(d,g,Q,n):
    x=d['x_end'];B=d['mass_ratio']*x
    S=d['checkpoints'][-1]['scalar_pole_over_T']
    labs=['parent','fermion']
    h=[np.asarray(d[f'endpoint_{a}_momentum_edges']) for a in labs]
    fs=[np.asarray(d[f'endpoint_{a}_occupations']) for a in labs]
    J=[cell_J(e,f,m,.5) for e,f,m in zip(h,fs,[x,B],strict=True)]
    AP,AB=.5*g*g/np.pi**2*(photon(1.,0.,.5)+np.asarray(J))
    E=np.hypot(S,Q)
    lo=brentq(lambda p:curve(p,x,AP)-E-curve(p+Q,B,AB),.2,4.,xtol=5e-14)
    hi=brentq(lambda p:curve(p,x,AP)-E-curve(p-Q,B,AB),.2,4.,xtol=5e-14)
    def inverse(e,M,A):
        t=(e+np.sqrt(np.maximum(e*e-4*A,0.)))/2
        return np.sqrt(np.maximum((t-M)*(t+M),0.))
    zp=h[0]
    zb=inverse(E+curve(h[1],B,AB),x,AP)
    cuts=np.unique(np.r_[lo,zp[(zp>lo)&(zp<hi)],zb[(zb>lo)&(zb<hi)],hi])
    z,w=roots_legendre(n)
    theta=np.arcsin(np.clip(2*(cuts-lo)/(hi-lo)-1.,-1.,1.))
    dt=np.diff(theta)[:,None]/2
    ang=(theta[1:]+theta[:-1])[:,None]/2+dt*z
    p=lo+(hi-lo)*(.5+.5*np.sin(ang.ravel()))
    measures=(dt*w*np.cos(ang)).ravel()/2
    P=curve(p,x,AP);D=P-E
    k=inverse(D,B,AB)
    mu=(p*p+Q*Q-k*k)/(2*p*Q)
    assert max(abs(mu))<=1.
    scalar=np.column_stack([Q*np.sqrt(1-mu*mu),Q*mu])
    parent=np.column_stack([np.zeros_like(p),p])
    daughter=parent-scalar
    assert max(abs(np.linalg.norm(daughter,axis=1)-k))<1e-13
    ip,ib=(np.searchsorted(a,v)-1 for a,v in zip(h,[p,k],strict=True))
    Fp,Fb=fs[0][ip],fs[1][ib]
    # Only the aperture geometry is tested. No spin-matrix and resonance
    # residues are supplied; use a positive phase measure common to signs.
    V=np.hypot(k,B)
    dv=k/V*(1-AB/V**2)
    Jaco=p*k/(P*D*dv)
    flux=measures*Jaco*(Fp*(1-Fb)+1000*(Fp-Fb))
    mom=lambda q,e,p:np.column_stack([q,q*e,q[:,None]*p])
    origin=mom(-flux,P,parent)
    end=mom(flux,D,daughter)
    scal=mom(flux,E,scalar)
    storage=[]
    for a,b,ii in zip(h,[origin,end],[ip,ib],strict=True):
        grid=np.zeros((len(a)-1,4))
        np.add.at(grid,ii,b)
        storage.append(grid)
    scalar_tot=scal.sum(axis=0)
    a,b=[z.sum(axis=0) for z in storage]
    residual=a+b+scalar_tot
    ndef=max(abs(a[0]+b[0]),abs(a[0]+scalar_tot[0]))
    pnorm=max(np.sum(abs(flux)),1e-30)
    assert max(abs(residual[1:]))/pnorm<8e-15 and ndef/pnorm<8e-15
    fp,fb=expit(-(P-.12)/.7),expit(-(D-.12)/.7)
    fs_=1/np.expm1(E/.7)
    fwd=fp*(1-fb)*(1+fs_)
    bwd=(1-fp)*fb*fs_
    thermal=max(abs(fwd/bwd-1))
    assert thermal<2e-14
    return dict(g=g,q=Q,nodes=n,number_defect_per_norm=ndef/pnorm,
                energy_defect_per_norm=residual[1]/pnorm,
                momentum_defect_per_norm=(residual[2:]/pnorm).tolist(),
                total_scalar_moments=scalar_tot.tolist(),thermal_max_relative=thermal)


def run():
    rows=[]
    for grid in ('fine','third'):
        d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
        for q in (1e-6,4.1e-6,1e-5):
            a,b=one(d,.005,q,16),one(d,.005,q,32)
            defect=max(abs(np.asarray(a['total_scalar_moments'])-np.asarray(b['total_scalar_moments'])))
            assert defect<2e-8,defect
            rows.append(dict(grid=grid,flux_refinement_absolute_defect=defect,**b))
    return dict(scope='Frozen separate moment deposits, not a dynamical multi-moment closure or actual rate.',
                rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
