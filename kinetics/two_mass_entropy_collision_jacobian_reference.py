"""Independent derivatives and event-space solve for the entropy collision chart."""
import json

import numpy as np

from two_mass_entropy_extent_reference import DT, batch, exponents, rates


def derivative(n,g,bose,V):
    a=exponents(n,g,bose)
    f=n/g
    H=1/(g*f*(1+f*np.where(bose,1,-1)))
    F=np.zeros((len(n),len(n)),dtype=DT)
    for v in V:
        x=a[0]
        y=(v*a)@bose
        z=(v*a)@(~bose)-v[0]*a[0]
        B=np.exp(-y)/(-np.expm1(-y))
        Ff=np.exp(-np.logaddexp(DT(0),z))
        one_minus_Fp=np.exp(-np.logaddexp(DT(0),-x))
        alpha,flux=map(DT,rates(n,g,bose,v))
        prev=one_minus_Fp*Ff*B
        ell=(B*v*bose-Ff*v*~bose)*H
        ell[0]=H[0]*one_minus_Fp
        F+=np.outer(v,flux*ell-prev*v*H)
    return F

def check(e,m,b,P,T):
    K,g,bose,V=batch(e,m,b,P,(-.65,.21,.72))
    eta=np.where(bose,DT(-.4),DT(-.2))
    eta[0]=DT(-.4)+DT(-.2)-DT(e)/DT(T)
    f=1/(np.exp(K/DT(T)-eta)+np.where(bose,-1,1))
    n=f*g
    a=exponents(n,g,bose)
    H=1/(g*f*(1+f*np.where(bose,1,-1)))
    A=np.zeros((len(n),len(n)),dtype=DT)
    columns=[]
    for v in V:
        x=a[0]
        y=(v*a)@bose
        z=(v*a)@(~bose)-v[0]*a[0]
        M=np.exp(-np.logaddexp(DT(0),x)-np.logaddexp(DT(0),-z)-np.log(-np.expm1(-y)))
        A-=M*np.outer(v,v)*H
        columns.append(np.asarray(np.sqrt(M*H)*v,float))
    def rhs(n):
        return sum((r*rates(n,g,bose,r)[1] for r in V),np.zeros_like(n))
    poll=[]
    for i in range(len(n)):
        delta=DT(3e-5)*min(n[i],g[i]-n[i]) if not bose[i] else DT(3e-5)*n[i]
        z=np.zeros_like(n);z[i]=delta
        y=[rhs(n+k*z) for k in (-2,-1,1,2)]
        poll.append((y[0]-8*y[1]+8*y[2]-y[3])/(12*delta))
    err=float(np.max(abs(np.array(poll).T-A)/np.max(abs(A))))
    W=np.sqrt(H)
    B=np.asarray(W[:,None]*A/W,dtype=float)
    E=K.copy(); E[0]+=DT(e)
    scalar=np.asarray(bose,dtype=DT); fermi=(~bose)&(np.arange(len(n))>0)
    visible=[]
    invariants=np.array([E/max(E),scalar+np.eye(len(n),dtype=DT)[0],fermi+np.eye(len(n),dtype=DT)[0]])
    for u in invariants:
        visible.append(float(np.max(abs(u@A))/np.max(abs(A))))
    C=np.linalg.eigvalsh((B+B.T)/2)
    X=np.array(columns).T
    if np.max(abs(B+X@X.T))>3e-13:
        raise ValueError("Thermal Jacobian did not factor into the three event directions.")
    # Solve a manufactured implicit linear response in three event
    # coordinates rather than the 21 occupation coordinates.
    given=np.cos(.37*np.arange(len(n)))
    ledger=np.asarray(invariants/W,float)
    err3=0.
    maxdefect=0.
    for dt in (.01,1.,100.,1000.):
        u=given-dt*X@np.linalg.solve(np.eye(X.shape[1])+dt*X.T@X,X.T@given)
        expected=np.linalg.solve(np.eye(len(n))-dt*B,given)
        err3=max(err3,float(np.linalg.norm(u-expected)/np.linalg.norm(expected)))
        maxdefect=max(maxdefect,float(np.max(abs(ledger@(u-given)))/np.max(abs(ledger))))
        if np.linalg.norm(u)>np.linalg.norm(given)*(1+5e-14):
            raise ValueError("Positive detailed-balance operator increased entropy norm.")
    if err3>8e-12 or maxdefect>8e-13:
        raise ValueError("Event-space inverse did not match the dense conserved solve.")
    if err>2e-9 or max(visible)>2e-16 or max(C)>1e-13 or np.count_nonzero(C<-1e-9)!=3:
        raise ValueError('Collision Jacobian did not meet analytic/finite/invariant agreement.')
    off=n.copy()
    off[0]=g[0]*DT(.63)
    diff=derivative(off,g,bose,V)
    poll=[]
    for i in range(len(n)):
        v=DT(2e-5)*min(off[i],g[i]-off[i]) if not bose[i] else DT(2e-5)*off[i]
        z=np.zeros_like(off);z[i]=v
        y=[rhs(off+k*z) for k in (-2,-1,1,2)]
        poll.append((y[0]-8*y[1]+8*y[2]-y[3])/(12*v))
    err2=float(np.max(abs(np.array(poll).T-diff))/np.max(abs(diff)))
    if err2>3e-9: raise ValueError(f"Nonthermal Jacobian finite difference failed: {err2}.")
    off_ledger=float(np.max(abs(invariants@diff))/np.max(abs(diff)))
    if off_ledger>2e-16:
        raise ValueError("Nonthermal Jacobian did not preserve the left ledgers.")
    return dict(parent_excess=e,scaled_finite_difference_error=err,
                off_equilibrium_scaled_finite_difference=err2,
                normalized_left_ledger_defects=visible,negative_eigenvalues=C[:3].tolist(),
                off_equilibrium_left_ledger_defect=off_ledger,
                maximum_event_space_solve_relative_error=err3,maximum_event_space_solve_ledger_defect=maxdefect)


if __name__=='__main__':
    print(json.dumps(dict(scope="Manufactured Jacobian regression; detailed-balance spectral sign not assumed away from equilibrium.",
        cases=[check(*r) for r in ((.093,.061,.122,.33,.07),(2e-23,.061,.122,4e-12,4e-23))]),indent=2,allow_nan=False))
