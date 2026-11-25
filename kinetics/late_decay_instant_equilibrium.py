"""Opposite limiting closure: immediate kinetic equilibrium of daughter B.

NR A->B+eta, eta held empty; B excited BE gas plus condensate at q=0.
Use number/kinetic energy closure, exact angular thermal Bose integrals,
positive finite-volume A loss and condensate stimulated source Nc*fA(b).
No actual C22/C12 equilibration proof or coherent condensate treatment.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.special import gamma,zeta,spence,roots_legendre
from scipy.integrate import BDF
from scipy.interpolate import PchipInterpolator
from two_state_halo_equilibrium import pl
from late_decay_bose_screen import mesh,vgrid
from direct_activation_consistency import T0
from late_activation_overlap import H
from resonant_mediator_inverse import M,C,SEC_GYR

C3=gamma(1.5)/(4*np.pi**2);C5=gamma(2.5)/(4*np.pi**2)
Z3=zeta(1.5,1);Z5=zeta(2.5,1)
a_grid=np.r_[0.,np.geomspace(1e-12,120.,2000)]
l3=np.log(pl(1.5,a_grid));l5=np.log(pl(2.5,a_grid))
logR=np.log(C3)+l3-.6*(np.log(C5)+l5)
def equilibrium(N,K):
    if N<=0:return 0.,np.inf,0.
    if K<=0:return 0.,0.,N
    T=(K/(C5*Z5))**.4;Ns=C3*Z3*T**1.5
    if N>=Ns:return T,0.,N-Ns
    R=np.log(N)-.6*np.log(K)
    if R<logR[-1]:
        # Boltzmann limit.
        T=2*K/(3*N);a=np.log(C3)+1.5*np.log(T)-np.log(N)
    else:
        a=np.interp(R,logR[::-1],a_grid[::-1])
        T=np.exp(.4*(np.log(K/C5)-np.interp(a,a_grid,l5)))
    return T,a,0.

def thermal(q,b,T,a):
    if T<=0 or not np.isfinite(a):return np.zeros_like(q),np.zeros_like(q)
    L=(q-b)**2/T;U=(q+b)**2/T
    yl=L+a;yu=U+a
    lo=np.log(np.maximum(-np.expm1(-yl),1e-300));hi=np.log(np.maximum(-np.expm1(-yu),1e-300))
    width=4*q*b/T
    # Average over angle: uniform in q_B^2.
    F=(hi-lo)/np.maximum(width,1e-300)
    G=T*(U*hi-L*lo-spence(-np.expm1(-yu))+spence(-np.expm1(-yl)))/np.maximum(width,1e-300)
    mid=(L+U)/2
    small=(width<1e-5*np.maximum(1.,mid))|(F<0)|(G<0)
    ff=np.exp(-mid-a)/np.maximum(-np.expm1(-mid-a),1e-300)
    F=np.where(small,ff,F);G=np.where(small,T*mid*ff,G)
    return np.maximum(F,0.),np.maximum(G,0.)

def run(row,recoil,lifetime=10.,n=100,cut=20.,rend=1.,rtol=2e-5,wall=150.):
    start=time.monotonic();p0=T0*(3.909/row['g_bath'])**(1/3)/row['final']['theta']*row['final']['peak_bin_q']
    rmax=min(1e8,.03*M/(cut*p0))
    e,_,f0,meta=mesh(row,n,cut,0.)
    vol,q=vgrid(e);N0=float(vol@f0);K0=float((vol*q*q)@f0)
    na=len(q);norm=np.maximum(f0,1.);rate=1/(lifetime*SEC_GYR)
    end=np.log(rmax/rend)
    u,w=roots_legendre(5)
    last_t=None
    def geometry(t):
        r=rmax*np.exp(-np.clip(t,0,end));b=recoil/(p0*r)
        # Split every parent cell at b; prevents log singular quadrature
        # node exactly at the condensate threshold. Integrate p^2 dq.
        lo,hi=e[:-1],e[1:];c=np.clip(b,lo,hi)
        L=np.stack([lo,c],axis=1);U=np.stack([c,hi],axis=1)
        nodes=(L[:,:,None]+U[:,:,None])/2+(U-L)[:,:,None]*u/2
        weights=(U-L)[:,:,None]*w/2*nodes*nodes/(2*np.pi**2)
        # Linear interpolation of f_A at the momentum selected by zero
        # recoil daughter. Use only nonnegative weights, conserve number.
        cp=np.zeros(na)
        if b<=e[-1]:
            if b<=q[0]:cp[0]=1.
            elif b>=q[-1]:cp[-1]=1.
            else:
                j=int(np.searchsorted(q,b)-1);v=(b-q[j])/(q[j+1]-q[j]);cp[j]=1-v;cp[j+1]=v
        return r,b,nodes,weights,cp,rate/H(r-1)
    def rates(t,y):
        r,b,qn,wt,cp,h=geometry(t)
        A=y[:na]*norm
        NB=max(y[-2]*N0,0.);KB=max(y[-1]*K0,0.)
        T,alpha,Nc=equilibrium(NB,KB)
        F,G=thermal(qn,b,T,alpha)
        # Collapsed intervals can yield q=0 and undefined intermediate
        # thermal values; their quadrature weights are identically zero.
        F=np.where(wt>0,F,0.);G=np.where(wt>0,G,0.)
        Wn=vol+np.sum(wt*F,axis=(1,2))
        Wk=vol*(q*q+b*b)+np.sum(wt*G,axis=(1,2))
        perA=h*(Wn+Nc*cp)
        return perA, h*Wk,(T,alpha,Nc,b,h,cp)
    def rhs(t,y):
        pn,pk,_=rates(t,y);A=y[:na]*norm
        return np.r_[-pn*A/vol/norm,(pn@A)/N0,(pk@A)/K0]
    def jac(t,y):
        pn,pk,_=rates(t,y)
        J=np.zeros((na+2,na+2));J[:na,:na]=np.diag(-pn/vol)
        J[-2,:na]=pn*norm/N0;J[-1,:na]=pk*norm/K0
        for i in [-2,-1]:
            # One-sided derivative at empty initial B, centered otherwise.
            step=max(1e-16,abs(y[i])*1e-5);p=y.copy();m=y.copy()
            p[i]+=step;m[i]=max(0.,y[i]-step)
            J[:,i]=(rhs(t,p)-rhs(t,m))/(p[i]-m[i])
        return J
    ini=np.r_[f0/norm,0.,0.]
    sol=BDF(rhs,0.,ini,end,jac=jac,rtol=rtol,atol=np.r_[np.full(na,1e-10),1e-16,1e-16],max_step=.035,first_step=1e-7)
    hist=[];first={};steps=0;msg='finished';last=start
    def stat(t,y):
        A=y[:na]*norm;N=y[-2]*N0;K=y[-1]*K0;T,a,Nc=equilibrium(N,K)
        r,b,_,_,cp,h=geometry(t)
        return dict(r=float(r),B_fraction=float(y[-2]),B_kinetic_over_initial=float(y[-1]),
          condensate_fraction_over_initial=float(Nc/N0),T_comoving_q2=float(T),alpha=float(a) if np.isfinite(a) else None,
          stimulated_condensate_source_over_H_per_condensate=float(h*(cp@A)),
          number_residual=float(vol@A/N0+y[-2]-1),negative_A_fraction=float(vol@np.maximum(-A,0.)/N0))
    hist.append(stat(0.,sol.y))
    while sol.status=='running':
        t0=sol.t;y0=sol.y.copy();note=sol.step();steps+=1
        v=stat(sol.t,sol.y)
        if v['negative_A_fraction']>1e-8 or sol.y[-2]<-1e-9 or sol.y[-1]<-1e-9:
            sol.t=t0;sol.y=y0;msg='negative accepted state';break
        if sol.status=='failed':msg=str(note);break
        for x in [1e-8,1e-6,.0001,.001,.01,.1,.5,.9]:
            if v['B_fraction']>=x and str(x) not in first:first[str(x)]=v
        if steps%5==0 or sol.status=='finished':hist.append(v)
        if time.monotonic()-last>12:
            print('progress',json.dumps(v),flush=True);last=time.monotonic()
        if time.monotonic()-start>wall:msg='wall budget';break
    if hist[-1]['r']!=float(rmax*np.exp(-np.clip(sol.t,0,end))):hist.append(stat(sol.t,sol.y))
    return dict(success=sol.status=='finished' and msg=='finished',message=msg,recoil_eV=recoil,kick_kms=recoil/M*C/1e5,
        lifetime_Gyr=lifetime,n=na,rmax=rmax,rend=rend,meta=meta,N0=N0,K0=K0,steps=steps,wall_seconds=time.monotonic()-start,
        first_crossings=first,history=hist,final=hist[-1],edges=e.tolist(),initial_A=f0.tolist(),final_A=(sol.y[:na]*norm).tolist(),
        scope='Immediate ideal BE+condensate closure, no eta occupation or inverse, no B mean field, no actual condensation proof, approximate zero-mode source interpolation; conditional limiting experiment.')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--input',required=True);a.add_argument('--output',required=True)
    a.add_argument('--recoil',nargs='+',type=float,default=[1e-9,1e-8,1e-7,1e-6,1e-5]);a.add_argument('--n',type=int,default=100);a.add_argument('--lifetime',type=float,default=10.);a.add_argument('--wall',type=float,default=120.);a.add_argument('--rend',type=float,default=1.);p=a.parse_args()
    row=json.loads(Path(p.input).read_text())[0];out=[]
    for k in p.recoil:
        v=run(row,k,p.lifetime,p.n,rend=p.rend,wall=p.wall);out.append(v)
        Path(p.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        print('FINAL',json.dumps({k:v for k,v in v.items() if k not in ['history','edges','initial_A','final_A']}),flush=True)
