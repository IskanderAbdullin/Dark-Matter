"""Bare-mass variation of the terminal scalar stress at frozen occupations."""
import argparse
import json
from types import SimpleNamespace

import numpy as np
from scipy.optimize import brentq

from early_pole_daughter_retained_kinetic_stress import build,r
from early_pole_live_pauli_gauge_moment_decay_reference import Tool


def state(q,n,lam,bare,scale=1.):
    q=q*scale
    n=n*scale**3
    total=float(sum(n))
    # At the root E>=m; a slightly enlarged this cubic bound brackets it.
    hi=2.*(bare+(.5*lam*total)**(1./3))
    eq=lambda m:m*m-bare*bare-.5*lam*(n@np.reciprocal(np.hypot(q,m)))
    gap=brentq(eq,max(1e-14,bare),hi,xtol=8e-18)
    E=np.hypot(q,gap)
    shift=gap*gap-bare*bare
    V=shift*shift/(2*lam)
    rho=float(n@E-V)
    P=float(n@(q*q/(3*E))+V)
    trace=float(bare*bare*(n@np.reciprocal(E)))
    if rho<=0. or not 0<P/rho<=1./3+2e-14 or abs(rho-3*P-trace)>5e-13*rho:
        raise ValueError('Massive stationary Hartree trace did not close.')
    return dict(bare_mass=float(bare),effective_gap=float(gap),pressure=P,energy=rho,
                pressure_to_energy=P/rho,trace_over_energy=trace/rho,counterterm=V)


def support(C,gap):
    """Reach only, using the original native hard geometry with a replaced scalar rest pole."""
    M=C.M
    e=M.weights(M.init,True)[0]
    inside=0
    any_h=0
    for q in C.q:
        model=SimpleNamespace(n=M.decay.n,cells=M.decay.cells,q=q,nod=M.z,w=M.w)
        active=list(Tool.events(model,dict(coeff=e[M.start:],es=float(np.hypot(q,gap)))))
        any_h+=bool(active)
        inside+=any(M.cells[k+M.start].lo>=1.4138888888888888-2e-14
                    and M.cells[k+M.start].hi<=1.4666666666666666+2e-14 for _,k,*_ in active)
    return dict(out_of=len(C.q),any_represented_decay_nodes=any_h,same_daughter_interval_nodes=inside)


def frozen_energy_cap(C,target_local):
    """Maximum kinematically available scalar energy, not a rate."""
    M=C.M
    e=M.weights(M.init,True)[0]
    ceilings=[]
    for q in C.q:
        cap=0.
        for i in range(M.r[1].start,M.r[1].stop):
            A=M.cells[i]
            for j in range(M.r[2].start,M.r[2].stop):
                B=M.cells[j]
                if target_local and (B.lo<1.4138888888888888-2e-14 or B.hi>1.4666666666666666+2e-14):
                    continue
                lo=max(A.lo,B.lo-q)
                hi=min(A.hi,B.hi+q)
                if hi<lo:
                    continue
                pts=np.array([lo,hi,np.clip(B.lo+q,lo,hi)])
                # For given p, daughter p' minimal is max(B.lo,p-q).
                # The constant daughter branch is increasing in p. On
                # the p-q branch, certify a conservative second derivative
                # lower bound on the complete interval.
                cut=max(lo,B.lo+q)
                if cut<hi:
                    lower=(e[i,1]/A.de)*A.m**2/np.hypot(hi,A.m)**3
                    upper=(e[j,1]/B.de)*B.m**2/np.hypot(cut-q,B.m)**3
                    if not lower>=upper:
                        raise ValueError('Convex support cap was not verified in this cell.')
                W=e[i,0]+e[i,1]*A.basis(pts)
                V=e[j,0]+e[j,1]*B.basis(np.maximum(B.lo,pts-q))
                cap=max(cap,float(max(W-V)))
        ceilings.append(np.sqrt(max((cap-q)*(cap+q),0.)))
    return float(max(ceilings)),float(min(ceilings))


def cap_diagnostic(C,n,local):
    last,first=frozen_energy_cap(C,local)
    q=np.r_[C.q,C.spq]
    E=np.hypot(q,last)
    bm2=last*last-.5*C.lam*float(n@np.reciprocal(E))
    if bm2<=0.:
        raise ValueError('The native support already excludes the scalar benchmark.')
    v=state(q,n,C.lam,np.sqrt(bm2))
    if abs(v['effective_gap']-last)>2e-17:
        raise ValueError('The exact kinematic cap did not reproduce in the scalar gap.')
    lo=support(C,last*(1.-1e-6))
    hi=support(C,last*(1.+1e-6))
    label='same_daughter_interval_nodes' if local else 'any_represented_decay_nodes'
    if lo[label]<=0 or hi[label]!=0:
        raise ValueError('Endpoint support calculation did not pass independent side checks.')
    return dict(last_node_effective_gap_ceiling=last,smallest_per_node_effective_gap_ceiling=first,
                bare_mass_at_last_node_ceiling=float(np.sqrt(bm2)),stress_at_ceiling=v,
                support_just_below=lo,support_just_above=hi)


def run(bundle,raw,angular):
    C,hard=build(bundle,angular)
    nodes=np.r_[C.q,C.spq]
    reference=C.scalar(C.ns)[-1]
    rows=[]
    ratios=(0.,.25,.5,1.,2.,3.,4.,8.,16.,32.,64.)
    for a,tag in zip(raw,('native','local-daughter'),strict=True):
        if a.get('prepared_charged_swap','native')!=tag or a.get('aperture_subdivision',1)!=hard:
            raise ValueError('Unmatched archived geometry for the frozen stress.')
        h=np.asarray(a['cold_mode_count_history'])
        np.testing.assert_array_equal(h[0],C.ns)
        np.testing.assert_array_equal(a['cold_momentum_nodes_over_T'],C.q)
        if len(h[-1])!=len(C.q) or not np.all(h[-1]>0.) or abs(a['energy_defect'])>=3e-12:
            raise ValueError('The terminal held scalar fixture does not pass its original gate.')
        counts=np.r_[h[-1],C.spect]
        def test(x):
            return state(nodes,counts,C.lam,x*reference)
        base=test(0.)
        if abs(base['energy']-C.scalar(h[-1])[1])>3e-17:
            raise ValueError('Massless control did not match the propagated chart.')
        sweep=[test(k) for k in ratios]
        if any(left['pressure_to_energy']<=right['pressure_to_energy'] for left,right in zip(sweep[:-1],sweep[1:],strict=True)):
            raise ValueError('Expected strictly decreasing stress ratio.')
        stage=[]
        for w in (.2,.1,.03,.01,.003,.001):
            hi=64.
            while test(hi)['pressure_to_energy']>w:
                hi*=2.
            m=brentq(lambda z:test(z)['pressure_to_energy']-w,0.,hi,xtol=5e-12)
            v=test(m)
            if abs(w-v['pressure_to_energy'])>2e-12:
                raise ValueError('Pressure threshold did not reproduce.')
            dd=1e-3
            t=[state(nodes,counts,C.lam,v['bare_mass'],np.exp(k*dd))['energy'] for k in (-2,-1,1,2)]
            P=(t[0]-8*t[1]+8*t[2]-t[3])/(36*dd)-v['energy']
            if abs(P-v['pressure'])>2e-7*v['pressure']:
                raise ValueError('Fixed bare mass pressure dilation failed.')
            stage.append(dict(target_pressure_to_energy=w,bare_to_initial_gap=float(m),state=v,
                              pressure_dilation=P,**support(C,v['effective_gap'])))
        rows.append(dict(arm=tag,massless=base,sweep=sweep,thresholds=stage,
                         local_aperture_cap=cap_diagnostic(C,counts,True),
                         whole_represented_hard_cap=cap_diagnostic(C,counts,False)))
    return dict(scope='Frozen terminal scalar counts, same scalar coupling, hypothetical bare masses; geometry reach uses the earlier native hard model.',
                hard_depth=hard,cold_modes=len(C.q),initial_massless_gap=reference,
                original_native_initial_reach=support(C,reference),arms=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('archive','native','local'):p.add_argument(k)
    p.add_argument('--angular-mesh')
    a=p.parse_args()
    print(json.dumps(run(a.archive,[r(a.native),r(a.local)],r(a.angular_mesh) if a.angular_mesh else None),
                     indent=2,allow_nan=False))
