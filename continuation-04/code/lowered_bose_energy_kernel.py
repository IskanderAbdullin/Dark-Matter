"""Analytically angle-integrated isotropic C22 and conservative energy lattice."""
import argparse,json,numpy as np
from lowered_antonov_hessian import Halo
from lowered_bose_signed_source import nodes,occupation
from lowered_bose_pointwise_c22 import c22
from lowered_bose_exterior_source import R,V,GYR

def exact_c22(H,psi,q,n=32):
    a=q*q;gap=max(H.eta-psi,1e-15);f1=float(occupation(H,psi,a));total=0.
    if a>=2*psi:return 0.
    x,w=nodes(n);pmax=2*psi-a;cuts=np.unique(np.r_[0,min(psi,pmax),pmax])
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        lg=np.log1p((hi-lo)/(gap+lo));b=lo+(gap+lo)*np.expm1(lg*x);wb=(gap+b)*lg*w;f2=occupation(H,psi,b)
        for bi,wi,fi in zip(b,wb,f2):
            S=a+bi;A=S/2;edge=max(0,S-psi);sub=np.unique([0,edge,min(a,bi),A]);sub=sub[(sub>=0)&(sub<=A)]
            for l,h in zip(sub[:-1],sub[1:]):
                lg2=np.log1p((h-l)/(gap+l));d=l+(gap+l)*np.expm1(lg2*x);wd=(gap+d)*lg2*w;c=S-d
                f3=occupation(H,psi,c);f4=occupation(H,psi,d);net=f3*f4*(1+f1+fi)-f1*fi*(1+f3+f4)
                ratio=np.minimum(1,np.sqrt(np.minimum(bi,d)/a)) if a>0 else np.ones_like(d)
                total+=wi*np.sum(wd*ratio*net)
    return float(total)

def lattice(H,psi,N=160,Emax=12,T=1.,dtmax=.02):
    de=Emax/N;E=(np.arange(N)+.5)*de;q=np.sqrt(E);D=q*de/2;f=occupation(H,psi,E);f0=f.copy();ids=[];weights=[]
    for s in range(2*N-1):
        i=np.arange(max(0,s-N+1),min(N-1,s//2)+1);j=s-i
        aa,bb=np.triu_indices(len(i),1);i1,j1,k1,l1=i[aa],j[aa],i[bb],j[bb]
        if len(aa)==0:continue
        ids.append(np.array([i1,j1,k1,l1]));weights.append(de**3*np.minimum(q[i1],q[k1])/4*np.where(i1==j1,1,2)*np.where(k1==l1,1,2)/2)
    i,j,k,l=np.concatenate(ids,axis=1);W=np.concatenate(weights);rho=V**2/(4*np.pi*6.6743e-11*R**2);W*=rho*.125*np.sqrt(2)*V*GYR/H.Ic
    def deriv(f):
        a,b,c,d=f[i],f[j],f[k],f[l];net=W*(a*b*(1+c+d)-c*d*(1+a+b));g=np.bincount(k,weights=net,minlength=N)+np.bincount(l,weights=net,minlength=N)-np.bincount(i,weights=net,minlength=N)-np.bincount(j,weights=net,minlength=N);return g/D
    def S(f):
        p=f>0;return D@(np.where(p,(1+f)*np.log1p(f)-f*np.log(np.maximum(f,1e-300)),0))
    def report(t,f):return dict(t=t,N=float(D@f),energy=float((D*E)@f),entropy=float(S(f)),f_first=float(f[0]),outside=float(D[E>=psi]@f[E>=psi]),min_f=float(f.min()),high_tail=float(D[E> .9*Emax]@f[E>.9*Emax]))
    eq=1/np.expm1(.3+E);eq_res=float(np.max(abs(deriv(eq))));c0=deriv(f);t=0.;rows=[report(t,f)];reject=0;steps=0;min_entropy_diff=0.
    for target in [.1,.5,1.]:
        if target>T+1e-12:continue
        while t<target-1e-12:
            g=deriv(f);neg=g<0;lim=np.min(-f[neg]/g[neg]) if neg.any() else np.inf;dt=min(dtmax,target-t,.4*lim)
            while True:
                p=f+dt*g
                if np.min(p)<0:dt*=.5;reject+=1;continue
                h=deriv(p);new=.5*f+.5*(p+dt*h)
                if np.min(new)<0 or S(new)<S(f)-1e-12:dt*=.5;reject+=1;continue
                break
            min_entropy_diff=min(min_entropy_diff,S(new)-S(f));f=new;t+=dt;steps+=1
        rows.append(report(t,f))
    return dict(model=H.model,psi=psi,N=N,Emax=Emax,events=len(W),steps=steps,reject=reject,initial_number_derivative=float(D@c0),initial_energy_derivative=float((D*E)@c0),initial_first_rate=float(c0[0]/f0[0]),initial_outside_rate=float(D[E>=psi]@c0[E>=psi]),min_entropy_diff=min_entropy_diff,BE_equilibrium_residual=eq_res,dtmax=dtmax,rows=rows)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--psi',type=float,default=3.9);p.add_argument('--N',type=int,default=120);p.add_argument('--Emax',type=float,default=12);p.add_argument('--T',type=float,default=1);p.add_argument('--compare',action='store_true');p.add_argument('--dtmax',type=float,default=.02);a=p.parse_args();H=Halo(model=a.model)
    if a.compare:
        for psi in [1.,3.,3.9,3.999999]:
            for q in [0.,.01,.3,1.,np.sqrt(psi)*1.1]:
                if q*q<2*psi:
                    x=exact_c22(H,psi,q,48);y=c22(H,psi,q,32);print(psi,q,x,y,(x-y)/max(abs(x),1e-30))
    else:print(json.dumps(lattice(H,a.psi,a.N,a.Emax,a.T,a.dtmax)))
