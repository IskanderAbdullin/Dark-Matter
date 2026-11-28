"""Retain early comoving distributions and test whether production has shut off."""
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

from early_pole_decay_mode_clock_reference import X
from early_pole_live_parent_energy_reference import state_moments
from early_pole_massive_daughter_energy_stencil_reference import MassiveFermion
from early_pole_resolved_parent_reference import S


def snapshot(d,xend):
    parent=expit(-d.parent_energies(X))
    init=np.r_[np.zeros(d.iP),parent]
    atol=np.r_[np.full(d.ns,1e-7),np.full(d.nf,8e-12),np.maximum(parent*5e-10,1e-19)]
    eta=np.log(xend/X)
    sol=solve_ivp(d.rhs,(0.,eta),init,method='BDF',jac=d.jac,
                  atol=atol,rtol=2e-8,max_step=.02,dense_output=True)
    assert sol.success,sol.message
    Q=d.d.p
    con=[]
    # Return flux across a cut based on whole cells; only 1e-5 and .1
    # are guaranteed constructed edges. Other choices use conservative brackets.
    def measure(x,y):
        E,D,P,m=state_moments(d,x,y)
        n=d.d.w*y[:d.ns]
        dot=d.d.w*d.rhs(np.log(x/X),y)[:d.ns]
        negative=-np.minimum(n,0.).sum()
        nn=np.maximum(n,0.)
        tot=nn.sum()
        assert negative<1e-10*tot
        cumulative=np.r_[0.,np.cumsum(nn)]/tot
        qtz={}
        for z in (.1,.5,.9,.99,.999):
            i=np.searchsorted(cumulative,z,side='left')
            qtz[str(z)]=d.d.ed[i-1:i+1].tolist()
        p_cut={}
        for p0 in (1e-6,1e-5,.01,.1):
            i=np.searchsorted(d.d.ed,p0,side='left')
            j=np.searchsorted(d.d.ed,p0,side='right')-1
            p_cut[str(p0)]=[float(cumulative[j]),float(cumulative[i])]
        bare=(x*1e-9)**2
        U=(m*m-bare)**2/(2*d.d.lam)
        eps=E[:d.ns]
        scalar_energy=nn@eps-U
        scalar_pressure=(nn@(Q*Q/(3*eps)))+U
        sum_d=d.d.wf@y[d.ns:d.iP]
        sum_p=d.pw@(y[d.iP:]-parent)
        con.append(max(abs(sum_d-n.sum()),abs(sum_p+n.sum())))
        return dict(x=x,scalar_yield=float(tot/S),negative_number_fraction=float(negative/tot),
                    scalar_pole_over_T=m,mean_q=float(nn@Q/tot),rms_q=float(np.sqrt(nn@(Q*Q)/tot)),
                    scalar_w=float(scalar_pressure/scalar_energy),
                    number_quantile_shell_brackets=qtz,below_cut_number_fraction_brackets=p_cut,
                    integrated_scalar_energy_over_T4=float(scalar_energy),
                    net_log_growth_rate=float(dot.sum()/tot),
                    total_absolute_cell_change_per_log_a=float(abs(dot).sum()/tot),
                    fraction_sub_1e_5_of_signed_net=float(dot[d.d.ed[1:]<=1e-5].sum()/dot.sum()))
    xx=np.unique(np.r_[np.array([.0075,.01,.02,.04,.05,.075,.1])[np.array([.0075,.01,.02,.04,.05,.075,.1])<xend],xend])
    rows=[measure(float(x),sol.sol(np.log(x/X))) for x in xx]
    assert max(con)<2e-12
    return dict(mass_ratio=d.ratio,elastic_rate=d.rate,x_end=xend,
                cells=[d.ns,d.nf,d.np],units='q=p/T_visible; visible T proportional to 1/a',
                maximum_absolute_sampled_pair_or_parent_number_defect_over_T3=float(max(con)),
                checkpoints=rows,
                scalar_momentum_edges=d.d.ed.tolist(),
                scalar_momentum_nodes=Q.tolist(),
                endpoint_scalar_signed_shell_number_over_T3=(d.d.w*sol.y[:d.ns,-1]).tolist(),
                endpoint_scalar_occupations=sol.y[:d.ns,-1].tolist(),
                endpoint_fermion_occupations=sol.y[d.ns:d.iP,-1].tolist(),
                endpoint_parent_occupations=sol.y[d.iP:,-1].tolist(),
                endpoint_fermion_momentum_edges=d.d.fed.tolist(),
                endpoint_parent_momentum_edges=d.pe.tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=65)
    p.add_argument('--fermion',type=int,default=90)
    p.add_argument('--parent',type=int,default=72)
    p.add_argument('--ratio',type=float,default=0.)
    p.add_argument('--rate',type=float,default=0.)
    p.add_argument('--xend',type=float,default=.05)
    a=p.parse_args()
    assert a.ratio>=0. and a.rate>=0. and .01<=a.xend<=.1
    d=MassiveFermion(a.scalar,a.fermion,a.parent,a.ratio,a.rate,True)
    print(json.dumps(snapshot(d,a.xend),indent=2,allow_nan=False))
