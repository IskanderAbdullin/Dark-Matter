"""Local nonlinear C22 initial Jacobian on conservative energy lattice."""
import argparse,json,numpy as np
from scipy.linalg import eigvals
from lowered_antonov_hessian import Halo
from lowered_bose_signed_source import occupation
from lowered_bose_exterior_source import R,V,GYR

def run(H,psi,N=120,Emax=12):
    de=Emax/N;E=(np.arange(N)+.5)*de;q=np.sqrt(E);D=q*de/2;f=occupation(H,psi,E);ids=[];weights=[]
    for s in range(2*N-1):
        i=np.arange(max(0,s-N+1),min(N-1,s//2)+1);j=s-i;aa,bb=np.triu_indices(len(i),1);i1,j1,k1,l1=i[aa],j[aa],i[bb],j[bb]
        if len(aa)==0:continue
        ids.append(np.array([i1,j1,k1,l1]));weights.append(de**3*np.minimum(q[i1],q[k1])/4*np.where(i1==j1,1,2)*np.where(k1==l1,1,2)/2)
    ids=np.concatenate(ids,axis=1);W=np.concatenate(weights);rho=V**2/(4*np.pi*6.6743e-11*R**2);W*=rho*.125*np.sqrt(2)*V*GYR/H.Ic
    def jac(f):
        a,b,c,d=f[ids];partials=[b*(1+c+d)-c*d,a*(1+c+d)-c*d,a*b-d*(1+a+b),a*b-c*(1+a+b)];J=np.zeros(N*N)
        for row,sign in zip(ids,[-1,-1,1,1]):
            for col,p in zip(ids,partials):J+=np.bincount(row*N+col,weights=sign*W*p,minlength=N*N)
        return J.reshape(N,N)/D[:,None]
    J=jac(f);ev=eigvals(J);eq=jac(1/np.expm1(.3+E));ee=eigvals(eq)
    return dict(model=H.model,psi=psi,N=N,Emax=Emax,max_real=float(ev.real.max()),min_real=float(ev.real.min()),max_imag=float(abs(ev.imag).max()),positive_count=int(np.sum(ev.real>1e-7)),number_left_residual=float(np.max(abs(D@J))),energy_left_residual=float(np.max(abs((D*E)@J))),equilibrium_max_real=float(ee.real.max()),equilibrium_min_real=float(ee.real.min()))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--psi',type=float,default=3.9);p.add_argument('--N',type=int,default=120);p.add_argument('--Emax',type=float,default=12);a=p.parse_args();print(json.dumps(run(Halo(model=a.model),a.psi,a.N,a.Emax)))
