"""Fixed-background Dirac Floquet cycles. Finite initial-vacuum and momentum-cutoff diagnostics only."""
import math
import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.special import ellipj, ellipk
m=86e-9; M1=.1; M2=.02; y=.001

def cycles(l,fac,p,ns=(1,2,4,8,16,32,64,128),start='turn',rtol=3e-9):
    p=np.asarray(p); N=len(p); ph=fac*math.sqrt(M1*M2)/y
    A=ph*math.sqrt(l)/(math.sqrt(12)*m); ep=A*A/(1+2*A*A)
    K=ellipk(ep); scale=m*math.sqrt(1+2*A*A)
    H0=np.zeros((N,4,4)); H0[:,0,0]=M1; H0[:,1,1]=M2
    H0[:,2,2]=-M1; H0[:,3,3]=-M2
    H0[:,0,2]=H0[:,2,0]=H0[:,1,3]=H0[:,3,1]=p
    B=np.zeros((4,4)); B[0,1]=B[1,0]=1; B[2,3]=B[3,2]=-1
    lo,hi=(0,2*K) if start=='turn' else (-K,K)
    def h(z): return H0+y*ph*ellipj(z,ep)[1]*B
    E,V=np.linalg.eigh(h(lo)); pos=V[:,:,2:]; neg=V[:,:,:2]
    Y0=np.broadcast_to(np.eye(4,dtype=complex),(N,4,4)).copy()
    def f(z,u): return (-1j/scale*h(z)@u.reshape(N,4,4)).ravel()
    s=solve_ivp(f,[lo,hi],Y0.ravel(),rtol=rtol,atol=rtol/100,method='DOP853')
    if not s.success: raise RuntimeError(s.message)
    U=s.y[:,-1].reshape(N,4,4); S=np.diag([1,-1,1,-1]); C=S@U@S@U
    norm=float(np.max(abs(C.conj().transpose(0,2,1)@C-np.eye(4))))
    def probs(U):
        B=pos.conj().transpose(0,2,1)@U@neg
        return np.sum(abs(B)**2,axis=2)
    occ={str(n):probs(np.linalg.matrix_power(C,n)) for n in ns}
    eig=np.linalg.eigvals(C)
    return E[:,2:],occ,norm,float(np.max(abs(abs(eig)-1))),s.nfev

def run(l,fac,ng=65,pmax=.12,start='turn',rtol=3e-9):
    p=np.linspace(0,pmax,ng); pieces=[]
    for i in range(0,ng,12): pieces.append(cycles(l,fac,p[i:i+12],start=start,rtol=rtol))
    E=np.concatenate([x[0] for x in pieces])
    occ={k:np.concatenate([x[1][k] for x in pieces]) for k in pieces[0][1]}
    ph=fac*math.sqrt(M1*M2)/y; rho=m*m*ph**2/2+l*ph**4/24
    n={k:(4/(2*np.pi**2)*simpson(p[:,None]**2*f,x=p,axis=0)).tolist() for k,f in occ.items()}
    enr={k:float(4/(2*np.pi**2)*simpson(np.sum(p[:,None]**2*f*E,axis=1),x=p)/rho) for k,f in occ.items()}
    return dict(ng=ng,pmax=pmax,start=start,n=n,energy_fraction=enr,
                maximum_occupation=max(float(np.max(f)) for f in occ.values()),
                norm=max(x[2] for x in pieces),eigerr=max(x[3] for x in pieces))

if __name__=='__main__':
    import json
    for start in ['turn','zero']:
        for n,pmax in [(65,.12),(129,.12),(257,.12),(257,.24)]:
            d=run(3.82587e-8,1.1,n,pmax,start)
            assert d['norm']<1e-7 and d['eigerr']<1e-7
            assert d['maximum_occupation']<=1+1e-6
            print(json.dumps(d),flush=True)
