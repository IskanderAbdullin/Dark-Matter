"""A manufactured unsplit backward step solved on conserved event extents."""
import json

import numpy as np
from scipy.optimize import least_squares

from two_mass_entropy_extent_reference import DT, batch, entropy, exponents, rates
from two_mass_entropy_collision_jacobian_reference import derivative


def event_slopes(n,g,bose,V):
    a=exponents(n,g,bose)
    f=n/g
    H=1/(g*f*(1+f*np.where(bose,1,-1)))
    rows=[]
    for v in V.T:
        alpha,flux=map(DT,rates(n,g,bose,v))
        x=a[0]
        y=(v*a)@bose
        z=(v*a)@(~bose)-v[0]*a[0]
        B=np.exp(-y)/(-np.expm1(-y))
        Ff=np.exp(-np.logaddexp(DT(0),z))
        unblocked=np.exp(-np.logaddexp(DT(0),-x))
        reverse=unblocked*Ff*B
        u=(B*v*bose-Ff*v*~bose)*H
        u[0]=unblocked*H[0]
        rows.append(flux*u-reverse*v*H)
    return np.array(rows,dtype=DT)


def solve(n,g,bose,V,h):
    records=[]
    x=np.zeros(V.shape[1])
    bounds=np.where(bose,np.inf,g)
    def flux(y):
        return np.array([rates(y,g,bose,v)[1] for v in V.T])
    def valid(y):
        return np.all(y>0) and np.all(y<bounds)
    # Homotopy is numerical only. Each rung solves the *same* start
    # state with a different h; it does not advance the time sequence.
    for dt in h*np.array([1/256,1/64,1/16,1/4,1.]):
        for j in range(48):
            Y=n+V@x
            f=x-dt*flux(Y)
            err=np.linalg.norm(f,np.inf)
            if err<3e-11:
                records.append(j)
                break
            M=np.eye(len(x))-dt*np.asarray(event_slopes(Y,g,bose,V)@V,float)
            dx=np.linalg.solve(M,-f)
            step=1.
            while step>=2.**-35:
                z=x+step*dx
                y=n+V@z
                if valid(y) and np.linalg.norm(z-dt*flux(y),np.inf)<err*(1-1e-4*step):
                    x=z
                    break
                step/=2
            else:
                raise ValueError("No residual-reducing feasible Newton step.")
        else:
            raise ValueError("The implicit nonlinear event extent was not resolved.")
    Y=n+V@x
    if not valid(Y):
        raise ValueError("The nonlinear coupled solver left its admissible domain.")
    return Y,x,records,float(np.linalg.norm(x-h*flux(Y),np.inf))


def case(title,e,m,b,P,T):
    K,g,bose,V=batch(e,m,b,P,(-.65,.21,.72))
    V=np.asarray(V,dtype=DT).T
    eta=np.where(bose,DT(-.4),DT(-.2))
    eta[0]=DT(-.4)+DT(-.2)-DT(e)/DT(T)
    n_eq=g/(np.exp(K/DT(T)-eta)+np.where(bose,-1,1))
    cases=[]
    for tag in ('forward','inverse','nearly_empty_daughters'):
        init=n_eq.copy()
        init[0]=g[0]*DT(.7 if tag=='forward' else .002)
        if tag=='nearly_empty_daughters':
            init[0]=g[0]*DT(.9999)
            init[1:]=g[1:]*DT(1e-8)
        M=event_slopes(init,g,bose,V)
        reference=derivative(init,g,bose,list(V.T))
        err=float(np.max(abs(V@M-reference))/np.max(abs(reference)))
        if err>4e-14:
            raise ValueError("Event gradient did not match the independent full-space Jacobian.")
        rec=[]
        for h in (.08,1.,50.):
            nxt,xi,it,res=solve(init,g,bose,V,h)
            delta=nxt-init
            kinetic=delta@K-DT(e)*sum(xi)
            increases=entropy(nxt,g,bose)-entropy(init,g,bose)
            ar=np.array([rates(nxt,g,bose,v) for v in V.T],dtype=DT)
            bound=h*(ar[:,0]@ar[:,1])
            fin=np.arange(len(nxt))
            fermi=(~bose)&(fin>0)
            conserved=max(abs(delta[bose].sum()+delta[0]),abs(delta[fermi].sum()+delta[0]))
            if abs(kinetic)>3e-15*max(K) or conserved>3e-15 or increases<bound-2.5e-10 or bound<0.:
                raise ValueError("Positive/entropy or split energy ledger failed.")
            err0=None
            if h==.08 and tag!='nearly_empty_daughters':
                def cost(y):
                    Y=np.asarray(y,dtype=DT)
                    f=np.array([rates(Y,g,bose,v)[1] for v in V.T])
                    return np.asarray(Y-init-h*V@f,float)
                opt=least_squares(cost,np.asarray(init,float),bounds=(
                    np.full(len(init),1e-18),np.asarray(np.where(bose,np.inf,g*(1-DT(1e-15))),float)),
                    xtol=4e-14,gtol=4e-14,ftol=4e-14)
                err0=float(max(abs(opt.x-nxt)))
                if not opt.success or max(abs(opt.fun))>2e-10 or err0>3e-10:
                    raise ValueError("Independent bounded occupation-space solution disagreed.")
            rec.append(dict(step=h,continuation_iterations=it,event_residual=res,
                            extensive_entropy_gain=float(increases),implicit_entropy_lower_bound=float(bound),
                            event_energy_relative_defect=float(abs(kinetic)/max(K)),
                            count_defect=float(conserved),independent_occupation_solve_absolute=err0))
        cases.append(dict(start=tag,steps=rec))
    return dict(name=title,capacity_count=len(g),event_count=V.shape[1],experiments=cases)


if __name__=='__main__':
    print(json.dumps(dict(scope="Manufactured nonlinear coupled implicit step, no physical source or convergence statement.",
        cases=[case('ordinary',.093,.061,.122,.33,.07),case('near-rest',2e-23,.061,.122,4e-12,4e-23)]),
                     indent=2,allow_nan=False))
