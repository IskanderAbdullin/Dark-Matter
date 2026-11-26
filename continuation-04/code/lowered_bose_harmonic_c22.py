"""Rotation-covariant weak C22 linearization on the paper-23 energy lattice.
Analytic P moments and exact global-rotation/relative-azimuth average.
Instantaneous finite-grid Jacobian, not nonlinear angular transport.
"""
import argparse,json
import numpy as np
from lowered_bose_quadrupole import occ

def matrices(N,Emax,psi,model,ls=(0,1,2)):
    de=Emax/N; E=(np.arange(N)+.5)*de;q=np.sqrt(E);D=q*de/2
    if model=='MB': f=np.exp(-E); fp=f.copy()
    elif model=='BE': f=1/np.expm1(.3+E);fp=f*(1+f)
    else:f,fp=occ(E,psi,model)
    mats={l:np.zeros((N,N)) for l in ls}; nquart=0
    for s in range(2*N-1):
        pair=np.arange(max(0,s-N+1),min(N-1,s//2)+1);other=s-pair
        aa,bb=np.triu_indices(len(pair),0); ids=np.array([pair[aa],other[aa],pair[bb],other[bb]]);a,b,c,d=f[ids];ei=E[ids];qi=q[ids];nquart+=len(aa)
        w=de**3*np.minimum(qi[0],qi[2])/4*np.where(ids[0]==ids[1],1,2)*np.where(ids[2]==ids[3],1,2)/2*np.where(aa==bb,.5,1.)
        if model=='MB': deriv=np.array([b,a,-d,-c])
        else:deriv=np.array([b*(1+c+d)-c*d,a*(1+c+d)-c*d,a*b-d*(1+a+b),a*b-c*(1+a+b)])
        S=ei[0]+ei[1];delta=ei-ei[[1,0,3,2]]
        p0=np.maximum(qi[1]-qi[0],qi[3]-qi[2]);p1=np.minimum(qi[0]+qi[1],qi[2]+qi[3])
        m2=(p0*p0+p0*p1+p1*p1)/3
        m4=(p0**4+p0**3*p1+p0*p0*p1*p1+p0*p1**3+p1**4)/5
        inv2=np.divide(1,p0*p1,out=np.zeros_like(p0),where=p0>0)
        inv4=np.divide(m2,p0**3*p1**3,out=np.zeros_like(p0),where=p0>0)
        ap2=(m2+2*delta+delta**2*inv2)/4
        signs=(-1,-1,1,1)
        for j in range(4):
            for i in range(4):
                if i==j:dot=ei[i];dot2=ei[i]**2
                elif i//2==j//2:dot=(m2-S)/2;dot2=(m4-2*S*m2+S*S)/4
                else:
                    x,y=delta[i],delta[j]
                    dot=(m2+x+y+x*y*inv2)/4
                    prod=(m4+2*(x+y)*m2+x*x+y*y+4*x*y+2*x*y*(x+y)*inv2+x*x*y*y*inv4)/16
                    dot2=1.5*prod+.5*(ei[i]*ei[j]-ei[i]*ap2[j]-ei[j]*ap2[i])
                cov={0:1.,1:dot/3,2:(3*dot2-ei[i]*ei[j])/10}
                idx=ids[j]*N+ids[i];v=w*signs[j]*deriv[i]
                for l in ls:mats[l]+=np.bincount(idx,weights=v*cov[l],minlength=N*N).reshape(N,N)
    return E,D,f,fp,mats,nquart

def run(N,Emax,psi,model):
    E,D,f,fp,K,nq=matrices(N,Emax,psi,model);out=dict(N=N,Emax=Emax,psi=psi,model=model,quartets=nq,modes={})
    x=-2*fp
    for l,M in K.items():
        mass=D*E**l/(2*l+1);L=M/mass[:,None]
        rate=-float(np.sum(M@x)/(mass@x))
        eig=np.linalg.eigvals(L)
        info=dict(moment_response=rate,max_real=float(eig.real.max()),min_real=float(eig.real.min()),max_imag=float(abs(eig.imag).max()),positive=int(np.sum(eig.real>1e-9)),sum_weak_residual=float(abs(M.sum(0)).max()),energy_weak_residual=float(abs(E@M).max()))
        if model in ('MB','BE'):
            chi=f if model=='MB' else f*(1+f);sym=M*chi[None,:];info.update(detailed_balance_asym=float(abs(sym-sym.T).max()/abs(sym).max()),largest_sym_eig=float(np.linalg.eigvalsh((sym+sym.T)/2).max()),boost_residual=float(abs(M@chi).max()))
        out['modes'][l]=info
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--N',type=int,default=80);p.add_argument('--Emax',type=float,default=12);p.add_argument('--psi',type=float,default=3.9);p.add_argument('--model',default='linear');a=p.parse_args();print(json.dumps(run(a.N,a.Emax,a.psi,a.model)))
