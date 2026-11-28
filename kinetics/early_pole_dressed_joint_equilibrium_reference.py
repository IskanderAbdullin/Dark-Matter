"""Common detailed-balance fixture and tangent test for the dressed hard network."""
import argparse
import json

import numpy as np
from scipy.special import expit
from scipy.linalg import null_space
from early_pole_live_joint_hard_packets_frw_reference import Expanding


def fixture(M,beta,mu):
    state=M.init.copy()
    Ns=M.Ns
    for it in range(150):
        E=M.weights(state,True)[2]
        nxt=[]
        for i,c in enumerate(M.cells):
            A,B=beta*(E[i]-[mu,0.])
            if i<M.r[0].stop:
                occ=1/np.expm1(A+B*c.u)
                v=c.moments(occ)
            else:
                occ=expit(-(A+B*c.u))
                v=c.I@occ
            nxt.append(v)
        nxt=np.array(nxt)
        es=M.scalar(Ns)[0]
        m=M.Ws/np.expm1(beta*es)
        diff=max(np.max(abs(nxt-state)),abs(m-Ns))
        state=.3*state+.7*nxt
        Ns=.3*Ns+.7*m
        if diff<3e-15:break
    assert diff<3e-14
    return np.r_[state.ravel(),Ns],it+1


def evaluate(M,y,P,D):
    M.pair_scale=P
    z=M.rhs_joint(np.r_[y,0.,0.,0.,0.],True,D)
    return z[:len(y)]


def grad(M,y):
    E,F,G,R,S=M.weights(y[:-1].reshape(-1,2),True)
    ex,_,_,f=M.scalar(y[-1])
    sigma=np.ones(len(F))
    sigma[M.r[0].stop:]*=-1
    return np.r_[(sigma[:,None]*F).ravel(),np.log1p(1/f)],np.r_[G.ravel(),ex]


def tangent(M,y,beta):
    F,E,G,_,_=M.weights(y[:-1].reshape(-1,2),True)
    S0,R0=grad(M,y)
    es,rho,entropy,fs=M.scalar(y[-1])
    k=np.r_[np.repeat(y[:-1].reshape(-1,2)[:,0],2),y[-1]]
    rng=np.random.default_rng(721)
    d=rng.uniform(-.75,.75,len(y))*k
    weights,scalar=M.weights,M.scalar
    def action(a,h):
        def at(s):
            fields=(S0+s*h*a)[:-1].reshape(-1,2)
            ff=fields.copy()
            ff[M.r[0].stop:]*=-1
            sc=1/np.expm1((S0+s*h*a)[-1])
            M.weights=lambda *_args:(F,ff,G,0.,0.)
            M.scalar=lambda *_args:(es,rho,entropy,sc)
            try:return evaluate(M,y,1.,1.)
            finally:M.weights=weights;M.scalar=scalar
        return (at(1)-at(-1))/(2*h)

    rows=[]
    C=np.r_[np.tile([1.,0.],len(E)),0.]
    for h in (2e-4,6e-5,2e-5):
        yp,ym=y+h*d,y-h*d
        pp,pm=grad(M,yp),grad(M,ym)
        Hs=(pp[0]-pm[0])/(2*h)
        Hr=(pp[1]-pm[1])/(2*h)
        truth=(evaluate(M,yp,1.,1.)-evaluate(M,ym,1.,1.))/(2*h)
        complete=action(Hs-beta*Hr,h)
        incomplete=action(Hs,h)
        scale=np.linalg.norm(truth)
        rows.append(dict(step=h,relative_l2_residual=np.linalg.norm(truth-complete)/scale,
                         omitted_curvature_relative_l2=np.linalg.norm(truth-incomplete)/scale,
                         energy_left_null_abs=abs(R0@truth),count_left_null_abs=abs(C@truth),
                         action_norm=scale))
    return rows


def spectrum(M,y,beta):
    sc=np.r_[np.repeat(y[:-1].reshape(-1,2)[:,0],2),y[-1]]
    S,R=grad(M,y)
    c=np.r_[np.tile([1.,0.],(len(y)-1)//2),0.]
    basis=null_space(np.vstack([R*sc,c*sc]))
    h=1e-5
    Js,Gs=[],[]
    for j in range(len(y)):
        d=np.eye(len(y))[j]*h*sc
        yp,ym=y+d,y-d
        gp,gm=grad(M,yp),grad(M,ym)
        Gs.append(-sc*((gp[0]-gm[0])-beta*(gp[1]-gm[1]))/(2*h))
        Js.append((evaluate(M,yp,1.,1.)-evaluate(M,ym,1.,1.))/(2*h*sc))
    J=np.column_stack(Js)
    G=np.column_stack(Gs)
    rec=np.linalg.norm(G-G.T)/np.linalg.norm(G)
    JP=basis.T@J@basis
    GP=basis.T@(G+G.T)@basis/2
    gp,U=np.linalg.eigh(GP)
    assert min(gp)>0.
    H=(U*np.sqrt(gp))@U.T
    Hinv=(U/np.sqrt(gp))@U.T
    L=-H@JP@Hinv
    defect=np.linalg.norm(L-L.T)/np.linalg.norm(L)
    ev=np.linalg.eigvalsh((L+L.T)/2)
    assert defect<1e-6 and min(ev)>-1e-7*max(ev)
    return dict(entropy_hessian_relative_symmetry_defect=rec,dimension=JP.shape[0],
                metric_jacobian_relative_symmetry_defect=defect,
                minimum_metric_eigenvalue=min(gp),relaxation_rates_sorted=ev.tolist(),
                rates_below_relative_1e_7=int(sum(ev<1e-7*max(ev))),
                conservation_left_defect=np.linalg.norm(np.vstack([R*sc,c*sc])@J))


def run(level):
    quad=((2,3,10,-3),(3,5,16,-4))[level]
    rows=[]
    for eta,beta,mu in [(0.,1.2,-.2),(.028,1.05,-.4)]:
        M=Expanding(12,quad,1.)
        M.eta=eta
        y,it=fixture(M,beta,mu)
        byfamily={}
        for P,D in [(0.,1.),(1.,0.),(1.,1.)]:
            byfamily[f'{P:g},{D:g}']=float(max(abs(evaluate(M,y,P,D))))
        s,r=grad(M,y)
        count=np.r_[np.tile([1.,0.],(len(s)-1)//2),0.]
        thermal=max(abs(s-beta*(r-mu*count)))
        assert max(byfamily.values())<1e-11 and thermal<2e-11
        a=tangent(M,y,beta)
        assert a[-1]['relative_l2_residual']<3e-7
        rows.append(dict(eta=eta,beta=beta,common_mu=mu,fixture_iterations=it,
                         stationarity_inf_by_pair_decay_strength=byfamily,thermal_tangent_inf=thermal,
                         tangent=a,projected_matrix=spectrum(M,y,beta) if level==0 else None))
    return dict(scope='Fixed-scale collision Jacobian, conditional on hard support and arbitrary clock.',
                quadrature=quad,rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,default=0)
    a=p.parse_args()
    print(json.dumps(run(a.level),indent=2,allow_nan=False))
