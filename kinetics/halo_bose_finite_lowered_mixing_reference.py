"""Finite lowered Bose halos: full-volume collisionless excess tests."""
import json
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.special import roots_legendre, expit


def solve(F,W,return_profile=False):
    d=brentq(lambda d:1/np.expm1(d)-1/np.expm1(d+W)-F,1e-12,4.,
             xtol=1e-15)
    h=1/np.expm1(d+W)
    xx,ww=roots_legendre(96)
    uu=(xx+1)/2
    def f(e):
        # Bose lowering, evaluated without cancellation near escape.
        return np.exp(d+e)*np.expm1(np.maximum(W-e,0.))/np.expm1(d+W)/np.expm1(d+e)
    def rho(y):
        u=np.sqrt(max(W-y,0.))*uu
        return float(np.sqrt(max(W-y,0.))*(ww@(u*u*f(y+u*u))))
    n=rho(0.)
    rg=np.unique(np.r_[0.,np.geomspace(d*1e-7,W,10_000),W-np.geomspace(1e-10,W/6,200)])
    fn=PchipInterpolator(rg,[rho(y)/n for y in rg])
    def diff(l,y):
        return [y[1],np.exp(2*l)*fn(np.clip(y[0],0.,W))-y[1]]
    def edge(l,y):
        return y[0]-W
    edge.terminal=True
    l0=-18.
    ep=np.exp(2*l0)
    s=solve_ivp(diff,(l0,20.),[ep/6,ep/3],events=edge,dense_output=True,
                rtol=3e-9,atol=1e-12,max_step=.04)
    assert s.success and len(s.t_events[0])==1
    bound=s.t_events[0][0]
    rt=np.exp(bound)
    norm=n*s.sol(bound)[1]*rt
    if return_profile:
        return d, h, n, norm, s.sol, bound
    # Use a separate quadrature over both full phase volume and lower excess.
    r,w=roots_legendre(18)
    v=np.linspace(l0,bound,450)
    lo,hi=v[:-1,None],v[1:,None]
    u=((hi+lo)/2+(hi-lo)/2*r).ravel()
    dw=((hi-lo)/2*w).ravel()
    Y=np.clip(s.sol(u)[0],0.,W)
    v3=dw*np.exp(3*u)
    et=W-Y
    vl=float(v3@((2/3)*et**1.5))
    chk=float(v3@np.array([rho(y) for y in Y]))
    assert abs(chk/norm-1)<3e-5, (F,W,chk,norm,rt)
    rows=[]
    for tau in (1e-10,1e-7,1e-4,.01,.1,.5,.9):
        level=F*tau
        emax=min(W,np.log1p(1/(level+h))-d)
        inner=np.maximum(emax-Y,0.)
        p=np.sqrt(inner[:,None])*uu
        up=np.sqrt(inner)*(np.maximum(f(Y[:,None]+p*p)-level,0.)*p*p@ww)
        z=float(v3@up)/norm
        eps=.045
        fzero=expit(1/eps)
        qmax=1+eps*np.log(1/(tau*fzero)-1)
        I=quad(lambda u:u*u*expit((1-u)/eps)/fzero,0.,4.,epsabs=1e-13)[0]
        J=quad(lambda u:u*u*expit((1-u)/eps)/fzero,0.,qmax,epsabs=1e-13)[0]
        pr=(J-tau*qmax**3/3)/I
        rows.append(dict(threshold_fraction=tau,halo_excess_per_particle=z,
                         primordial_excess_per_particle=pr,raw_retention_upper=pr/z,retained_fraction_upper=min(1.,pr/z)))
    # Sufficient analytic condition on the volume of primordial modes
    # above f*, to prove incompatibility with retention exactly one.
    qneed=(3*I*F*vl/norm)**(1/3)
    logthr=np.log(expit((1-qneed)/eps)/fzero) if qneed<20 else (1-qneed)/eps-np.log(fzero)
    return dict(cap=F,escape_depth=W,central_Bose_offset=d,
                tidal_radius_R=rt,total_number_dimensionless=norm,
                mass_closure_relative_defect=abs(chk/norm-1),
                finite_phase_volume_over_count=vl/norm,
                all_retained_obstructed_at_log_level_over_cap_below=logthr,
                sampled_levels=rows)


def run():
    return dict(scope="Single-component finite spherical lowered Bose stationary models, no baryons or condensate. "
                      "Necessary collisionless mixing tests against a uniform finite patch with compact-width "
                      "logistic primordial momentum shape, which has a positive infinite tail. No reachability, "
                      "stability, galaxy fit or actual retained fraction inferred.",
                configurations=[solve(F,W) for F,W in ((.512,5.),(7000.,2.),(7000.,5.),(7000.,8.))])


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
