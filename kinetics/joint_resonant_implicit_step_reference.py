"""Implicit entropy-variable collision step on a strictly populated test box."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import expit
from joint_resonant_kinetic_linear_modes_reference import configuration,metric
from sm_dark_entropy_exchange_reference import StandardModel

class Model:
    def __init__(self,ker,H):
        self.k=ker;self.H=H
        self.n=len(ker.c[0].vol);self.d=4*self.n
        self.us=np.clip((ker.k1-ker.c[0].low[ker.I])/ker.c[0].dk[ker.I],0.,1.)
        self.uc=np.clip((ker.k2-ker.c[1].low[ker.J])/ker.c[1].dk[ker.J],0.,1.)

    def arrays(self,q):
        return [q[2*j*self.n:2*(j+1)*self.n].reshape(self.n,2).T for j in (0,1)]

    def moment(self,q):
        ys=[];S=0.;metric=[]
        for ar,c in zip(self.arrays(q),self.k.c):
            z=ar[0,:,None]+ar[1,:,None]*c.u
            if c.g==1 and np.min(z)<=0.:raise ValueError("normal BE domain")
            f=np.exp(-z)/(-np.expm1(-z)) if c.g==1 else expit(-z)
            w=c.w*c.vol[:,None]
            n=(w*f).sum(axis=1)
            u=(w*f*c.u).sum(axis=1)
            ys.append(np.stack((n,u),axis=1).ravel())
            S+=np.sum(w*(z*f-(np.log(-np.expm1(-z)) if c.g==1 else -np.logaddexp(0.,-z))))
            wf=w*f*(1+(1 if c.g==1 else -1)*f)
            metric.append(np.stack([(wf*c.u**u).sum(axis=1) for u in (0,1,2)],axis=1))
        return np.concatenate(ys),S,np.concatenate(metric,axis=0)

    def proj(self,x):
        k=self.k
        return np.concatenate([np.stack((np.bincount(j,weights=x,minlength=self.n),
                              np.bincount(j,weights=u*x,minlength=self.n)),axis=1).ravel()
                           for j,u in ((k.I,self.us),(k.J,self.uc))])

    def collision(self,q):
        k=self.k
        th=self.arrays(q)
        f=k.c[0].f(th[0],k.ph,k.I)
        g=k.c[1].f(th[1],k.ch,k.J)
        a=k.pre*f*g
        b=k.pre*(1+f)*(1-g)
        def sumparent(z):return np.bincount(k.y,weights=z,minlength=k.n)[k.y]
        A=sumparent(a);B=sumparent(b);tot=A+B
        rho=np.divide(k.fac,tot*self.H,out=np.zeros_like(tot),where=tot>0)
        val=b*A-a*B
        flux=self.proj(rho*val)
        def jvp(w):
            d1,d2=self.arrays(w)
            df=-f*(1+f)*(d1[0,k.I]+d1[1,k.I]*self.us)
            dg=-g*(1-g)*(d2[0,k.J]+d2[1,k.J]*self.uc)
            da=k.pre*(df*g+f*dg)
            db=k.pre*(df*(1-g)-(1+f)*dg)
            X=sumparent(da);Y=sumparent(db)
            z=db*A+b*X-da*B-a*Y
            return self.proj(rho*(z-val*np.divide(X+Y,tot,out=np.zeros_like(tot),where=tot>0)))
        return flux,jvp

def step(model,q,h):
    M0,S0,_=model.moment(q)
    rms=np.repeat(M0[::2],2)
    it=0
    best=np.inf
    for _ in range(18):
        M,S,L=model.moment(q)
        F,jvp=model.collision(q)
        g=(M-M0-h*F)/rms
        measure=np.max(abs(g))
        if measure<8e-12:
            assert S>=S0
            return q,dict(newton_iterations=it,residual=measure,entropy_gain_relative=float(S/S0-1))
        def j(w):
            a=w.reshape(-1,2)
            dl=-np.stack((L[:,0]*a[:,0]+L[:,1]*a[:,1],L[:,1]*a[:,0]+L[:,2]*a[:,1]),axis=1).ravel()
            return (dl-h*jvp(w))/rms
        # Small reference grid: a full analytic Jacobian avoids both
        # finite-difference cancellation and an underresolved inner solve.
        jac=np.stack([j(w) for w in np.eye(model.d)],axis=1)
        v=np.linalg.solve(jac,-g)
        best=np.linalg.norm(g)
        for i in range(12):
            cand=q+(2.**-i)*v
            try:A,SS,_=model.moment(cand)
            except ValueError:continue
            F1,_=model.collision(cand)
            err=np.linalg.norm((A-M0-h*F1)/rms)
            if err<best*(1.-(2.**-i)*1e-3):break
        else:raise RuntimeError("implicit residual not reduced")
        q=cand;it+=1
    raise RuntimeError("implicit iteration budget")

def explicit_reference(model,q,extent=3e-8,frames=42):
    M,_,_=model.moment(q)
    x=q.copy()
    h=extent/frames
    def shape(Y):
        out=[]
        for i,c in enumerate(model.k.c):
            old=Y[2*i*model.n:2*(i+1)*model.n].reshape(model.n,2)
            K=old[:,0]*c.low+old[:,1]*c.dk
            ans,err,_=c.reconstruct((old[:,0],K))
            assert err<2e-7
            out.append(np.asarray(ans).T.ravel())
        return np.concatenate(out)
    for _ in range(frames):
        d=model.collision(x)[0]
        trial=shape(M+h*d)
        M=M+.5*h*(d+model.collision(trial)[0])
        x=shape(M)
    return M


def test(r,eos):
    T,H,k,N,t,Q,theta=configuration(r,eos,7,7.0)
    mod=Model(k,H)
    n=mod.n
    u=np.concatenate([s.T.ravel() for s in theta])
    d=np.array([.1*np.sin(j*.33)+.03*np.cos(j*.74) for j in range(len(u))])
    q=u+d
    N0,S0,_=mod.moment(q)
    F,J=mod.collision(q)
    ds=np.sin(np.arange(len(q))*.76)
    eps=4e-6
    check=np.linalg.norm((mod.collision(q+eps*ds)[0]-mod.collision(q-eps*ds)[0])/(2*eps)-J(ds))/np.linalg.norm(J(ds))
    assert check<1e-6
    _,_,invs=metric(k,theta)
    store=[]
    for h in (5e-9,5e-8,2.5e-7):
        z,a=step(mod,q.copy(),h)
        M,S,_=mod.moment(z)
        residual=max(abs((invs.T@M)/(invs.T@N0)-1))
        assert residual<5e-11
        delta=M-N0
        heat=invs[:2*n,2]@delta[:2*n]
        store.append(dict(step_H=h,conserved_moment_error=residual,
            fractional_scalar_heat_of_initial_dark_energy=float(heat/(invs[:,2]@N0)),**a))
    # Same method, same integrated interval but finer steps; this is a
    # consistency rather than equivalence test for large backward steps.
    z=q.copy()
    for _ in range(6):z,_=step(mod,z,5e-9)
    fine=mod.moment(z)[0]
    z,_=step(mod,q.copy(),3e-8)
    coarse=mod.moment(z)[0]
    control=explicit_reference(mod,q.copy())
    be=np.linalg.norm(fine-control,ord=1)/np.linalg.norm(control,ord=1)
    z=q.copy()
    for _ in range(12):z,_=step(mod,z,2.5e-9)
    half=np.linalg.norm(mod.moment(z)[0]-control,ord=1)/np.linalg.norm(control,ord=1)
    assert half<be<.008
    return dict(bins=n,analytic_jacobian_relative_check=check,steps=store,
        BE_vs_explicit_heun_relative_l1=float(be),refined_BE_vs_explicit_heun_relative_l1=float(half),
        reduced_steps_terminal_state_relative_l1=float(np.linalg.norm(fine-coarse,ord=1)/np.linalg.norm(fine,ord=1)),
        scope="Strictly populated, manufactured near-thermal state on a finite support; an implicit collision-step prototype, not the fixed-background full evolved spectrum.")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument('eos');p.add_argument('spectra')
    a=p.parse_args();eos=StandardModel(a.eos)
    r=next(r for r in json.loads(Path(a.spectra).read_text())['rows'] if r['TR']==.05 and abs(r['TF']-.045)<1e-8)
    print(json.dumps(test(r,eos),indent=2,allow_nan=False))
