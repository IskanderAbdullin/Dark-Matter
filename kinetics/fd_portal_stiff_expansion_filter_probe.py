"""Probe stiff FD splitting while the equilibrium state moves with expansion."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_tabulated_entropy_split_expansion import Stage
from fd_portal_tabulated_exponential_filter import Filtered
from sm_dark_entropy_exchange_reference import StandardModel


class FilteredExpansion(Stage):
    def closed_bath(self,T,u):
        return Filtered(T,u,self.nodes,self.eos)


def cell(eos,q,n,T,eta,r):
    ref=Stage(eos,q,n)
    fil=FilteredExpansion(eos,q,n)
    a=np.array([T,eta,np.log(T*r)])
    rate=fil.closed_bath(T,a[1:])
    modes=np.linalg.eigvalsh(rate.matrices(q))
    H=ref.rhs(0.,a)[0]
    horizon=min(.014/H,44/modes[1])
    S,E=ref.entropy(0.,a),ref.energy(0.,a)

    def rhs(_,v):
        x=v[0]
        state=v[1:4]
        h,u,p=ref.rhs(x,state)
        P=eos.state(state[0])[1]+ref.bath(state[1:])['pressure']
        return horizon*np.r_[h,u,p/S,3*h*np.exp(3*x)*P/E]

    def solve(tol):
        m=solve_ivp(rhs,(0.,1.),np.r_[0.,a,0.,0.],method='DOP853',rtol=tol,atol=.01*tol)
        if not m.success:
            raise ArithmeticError('No physical continuous solution')
        x=m.y[:,-1]
        if max(abs(ref.entropy(x[0],x[1:4])/S-1.-x[4]),abs(ref.energy(x[0],x[1:4])/E-1.+x[5]))>2e-11:
            raise ArithmeticError('Unmodified differential budgets did not close')
        return x

    z,y=[solve(t) for t in (2e-12,6e-13)]
    if np.max(abs(z-y))>3e-10:
        raise ArithmeticError('Independent tolerance control failed')
    rows=[]
    for label,o in (('unfiltered',ref),('filtered',fil)):
        for k in (4,8,16,32):
            x=0.
            u=a.copy()
            ent=work=0.
            for _ in range(k):
                x,u,p,e,w,d=o.step(x,u,horizon/k)
                ent+=p
                work+=w
            f=np.r_[x,u]
            err=(f-y[:4])/np.array([.01,T,1.,1.])
            h=[abs(o.entropy(x,u)/S-1.-ent/S),abs(o.energy(x,u)/E-1.+work/E)]
            if max(h)>2e-11 or ent<=0. or f[0]<=0.:
                raise ArithmeticError('Stiff step broke the discrete expansion ledgers')
            rows.append(dict(method=label,count=k,initial_fast_rate_times_step=horizon*modes[1]/k,
                             terminal_coordinates=f.tolist(),component_scaled_errors=err.tolist(),
                             maximum_scaled_state_error=np.max(abs(err)),
                             ledger_abs_errors_entropy_energy=h,relative_entropy_minus_physical=(ent/S-y[4])))
    return dict(initial_visible_T=T,initial_state=a.tolist(),horizon=horizon,
                initial_local_mode_rates=modes.tolist(),continuous_state_and_ledgers=y.tolist(),steps=rows)


def run(path):
    eos=StandardModel(path)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T,eta,r in ((.64135,.002,.9995),(.6882487043,-.002,1.0005)):
        a,b=[cell(eos,o,n,T,eta,r) for o,n in ((p,112),(q,168))]
        diff=max(np.max(abs(np.asarray(x['terminal_coordinates'])-y['terminal_coordinates']))
                 for x,y in zip(a['steps'],b['steps'],strict=True))
        if diff>2e-7:
            raise ArithmeticError('Different phase-space grids do not support this stiff comparison')
        rows.append(dict(maximum_mesh_endpoint_difference=diff,**b))
    return dict(rows=rows,table_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Stiff near-equilibrium short expansion with purely timelike vacuum FD kinetics. Not a claim of unconditional accuracy or abundance.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    print(json.dumps(run(p.parse_args().eos),indent=2,allow_nan=False))
