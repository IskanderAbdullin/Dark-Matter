"""Elastic steps in admissible FD coordinates with exact discrete energy."""
import argparse
import json

import numpy as np
from scipy.linalg import block_diag
from scipy.special import expit
from early_pole_archived_cold_screened_elastic_reference import Screened


def approximate_jacobian(m,G,F,k,dt,cell_sizes,F_mid=None,F_end=None):
    M=m.M
    ix,u,pb,Q,W=m.graph
    d=M.r[2]
    v0=G[ix,0]+G[ix,1]*u
    target=v0[0]+v0[1]-v0[2]
    j=d.start+np.searchsorted(G[d,0]+G[d,1],target)
    inside=j<d.stop
    j=np.minimum(j,d.stop-1)
    v=(target-G[j,0])/G[j,1]
    p=M.cells[d.start].pfrombare(M.shape[j,0]+M.shape[j,1]*v)
    take=inside&(abs(v)<=1.)&(abs(pb-p)<=Q)&(Q<=pb+p)
    v=np.where(take,v,0.)
    A,B,C=expit(F[ix,0]+F[ix,1]*u)
    D=expit(F[j,0]+F[j,1]*v)
    # Forward-minus-reverse derivatives in the four event logits.
    Delta=np.array([A*(1-A)*(B*(1-C)*(1-D)+C*D*(1-B)),
                    B*(1-B)*(A*(1-C)*(1-D)+C*D*(1-A)),
                   -C*(1-C)*(A*B*(1-D)+D*(1-A)*(1-B)),
                   -D*(1-D)*(A*B*(1-C)+C*(1-A)*(1-B))])
    pos=np.stack([*ix,j])-k
    basis=np.stack([*u,v])
    b=np.stack([np.ones_like(basis),basis],axis=1)
    deriv=(b*Delta[:,None,:]).reshape(8,-1)
    deposit=(b*np.array([-1.,-1.,1.,1.])[:,None,None]).reshape(8,-1)
    ids=(2*pos[:,None,:]+np.array([0,1])[None,:,None]).reshape(8,-1)
    fw=m.factor*W*take*m.strength
    n=(len(F)-k)*2
    J=np.zeros((n,n))
    for a in range(8):
        for b in range(8):
            J+=np.bincount(ids[a]*n+ids[b],weights=fw*deposit[a]*deriv[b],
                           minlength=n*n).reshape(n,n)
    def susceptibility(state):
        blocks=[]
        for j,cell in enumerate(M.cells[k:],start=k):
            o=expit(state[j]@cell.B)
            blocks.append((cell.I*(o*(1-o)))@cell.B.T)
        return block_diag(*blocks)
    end=susceptibility(F if F_end is None else F_end)
    # Remap and secant derivatives are omitted only from the trial Jacobian.
    response=J if F_mid is None else .5*J@np.linalg.solve(susceptibility(F_mid),end)
    return (end-dt*response)/cell_sizes.ravel()[:,None]



def step(model,Y0,dt,method="backward"):
    M=model.M
    k=M.r[1].start
    cells=M.cells[k:]
    scales=np.repeat([c.W for c in cells],2).reshape(-1,2)
    _,theta0,_,rho0,S0=M.weights(Y0,True)
    if dt<=0. or method not in ("backward","midpoint"):
        raise ValueError("positive step with backward or midpoint method required")
    symmetric=method=="midpoint"
    def currents(Y):
        j0=sum(M.cells[j].fit(Y[j])[2] for j in range(k))
        v=[sum(M.cells[j].fit(Y[j])[1] for j in range(a.start,a.stop)) for a in M.r[1:]]
        return M.background+np.r_[j0,v]
    J0=currents(Y0)
    def partner(j):
        return np.r_[sum(j[1:]),j[0]+j[1:]]
    def material(theta):
        f=theta0.copy()
        f[k:]=theta.reshape(-1,2)
        Y=Y0.copy()
        for j,c in enumerate(cells,start=k):
            Y[j]=c.I@expit(f[j]@c.B)
        return Y,f
    def evaluate(theta):
        Y,F=material(theta)
        mid=(Y+Y0)/2
        _,F_mid,G,_,_=M.weights(mid,True)
        dy=Y-Y0
        v=np.sum(dy*dy)
        if v>1e-26:
            J=currents(Y)
            er=np.sum(M.shape*dy)+.5*M.c*(J+J0)@partner(J-J0)
            G=G+(er-np.sum(G*dy))/v*dy
        F_step=F
        if symmetric:
            _,_,_,_,S1=M.weights(Y,True)
            F_step=F_mid if v<=1e-26 else F_mid+((-S1+S0)-np.sum(F_mid*dy))/v*dy
        D,S,_,heat=model.collision(G,F_step,G)
        return ((dy-dt*D)[k:]/scales).ravel(),Y,F_step,G,S,heat,F_mid,F
    if symmetric:
        predictor,_=step(model,Y0,dt,"backward")
        theta=M.weights(predictor,True)[1][k:].ravel().copy()
    else:
        theta=theta0[k:].ravel().copy()
    current=evaluate(theta)
    history=[]
    for _ in range(80 if symmetric else 56):
        r,Y,F,G,_,_,Fm,Fe=current
        error=float(np.linalg.norm(r,np.inf))
        history.append(error)
        if error<1e-13:
            break
        J=approximate_jacobian(model,G,F,k,dt,scales,Fm if symmetric else None,Fe)
        direction=np.linalg.solve(J,-r)
        c=1.
        for it in range(22):
            trial=theta+c*direction
            if np.max(abs(trial))<45.:
                try:
                    nxt=evaluate(trial)
                except (AssertionError,ValueError,np.linalg.LinAlgError):
                    pass
                else:
                    if np.linalg.norm(nxt[0],np.inf)<(1.-(.001 if symmetric else .05)*c)*error:
                        theta,current=trial,nxt
                        break
            c*=.5
        else:
            raise RuntimeError((dt,history))
    r,Y,F,G,edot,heat,_,_=current
    assert max(abs(r))<5e-13,(dt,history,max(abs(r)))
    _,_,_,r1,s1=M.weights(Y,True)
    count=[np.sum(Y[z,0]-Y0[z,0]) for z in M.r[1:]]
    assert abs(r1-rho0)<2e-12 and max(abs(x) for x in count)<2e-12,(r1-rho0,count)
    if symmetric:
        assert abs(s1-S0-dt*edot)<2e-12
    else:
        assert s1-S0>=dt*edot-2e-12
    return Y,dict(h=dt,history=history,rho_drift=r1-rho0,count_drift=count,
                  entropy_gain=s1-S0,entropy_lower_bound=dt*edot,end_daughter_heat=heat)


def run(grid,dt,method="backward"):
    model=Screened(grid,1,8,dt,(1,1,2),.08,1.)
    Y0=model.M.init
    Y1,r1=step(model,Y0,dt,method)
    Yh,r2=step(model,Y0,dt/2,method)
    Y2,r3=step(model,Yh,dt/2,method)
    return dict(grid=grid,method=method,one=r1,halves=[r2,r3],halving_maximum_moment_error=float(np.max(abs(Y1-Y2))))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("grid",choices=("fine","third"))
    p.add_argument("--dt",type=float,default=8e-8)
    p.add_argument("--method",choices=("backward","midpoint"),default="backward")
    a=p.parse_args()
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
