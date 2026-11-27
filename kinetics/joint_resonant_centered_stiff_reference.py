"""Common-affinity-centered implicit prompt collision step for stiff limits."""
import argparse,json
from pathlib import Path
import numpy as np
from joint_prompt_two_moment_evolution_reference import maxent
from joint_resonant_kinetic_linear_modes_reference import configuration,metric
from joint_resonant_implicit_step_reference import Model,step as previous
from sm_dark_entropy_exchange_reference import StandardModel

class Centered(Model):
    def reference(self,q,guess):
        M,_,_=self.moment(q)
        c=self.k.c; n=self.n
        N=M.reshape(2,n,2)[:,:,0].sum(axis=1)
        E=sum(M[2*j*n:2*(j+1)*n:2]@(s.mass+s.low)+M[2*j*n+1:2*(j+1)*n:2]@s.dk
              for j,s in enumerate(c))
        t,r=maxent(c,N,E,guess)
        qref=np.concatenate([np.array([(s.mass+s.low)/t-r[j][3],s.dk/t]).T.ravel()
                             for j,s in enumerate(c)])
        R0=np.exp(r[0][3]+r[1][3]-self.k.E[self.k.y]/t)
        return qref,R0,t

    def centered(self,z,qref,R0):
        k=self.k
        a0,a1=self.arrays(qref+z)
        f=k.c[0].f(a0,k.ph,k.I)
        g=k.c[1].f(a1,k.ch,k.J)
        b=k.pre*(1+f)*(1-g)
        s,t=self.arrays(z)
        shift=s[0,k.I]+s[1,k.I]*self.us+t[0,k.J]+t[1,k.J]*self.uc
        md=np.expm1(-shift)
        def accum(v):return np.bincount(k.y,weights=v,minlength=k.n)[k.y]
        B=accum(b);J=accum(b*md)
        D=B+R0*(B+J)
        omega=np.divide(k.fac*R0,D*self.H,out=np.zeros_like(D),where=D>0)
        v=J-B*md
        flux=self.proj(omega*b*v)
        def derivative(d):
            p,q=self.arrays(d)
            h=p[0,k.I]+p[1,k.I]*self.us
            l=q[0,k.J]+q[1,k.J]*self.uc
            db=b*(-f*h+g*l)
            mm=-(1+md)*(h+l)
            w=accum(db);jj=accum(db*md+b*mm)
            change=jj-w*md-B*mm
            dD=w+R0*(w+jj)
            return self.proj(omega*(v*db+b*change-b*v*np.divide(dD,D,out=np.zeros_like(D),where=D>0)))
        return flux,derivative

def advance(model,q,h,guess,reference_guess=False,vacuum_scaled=False,
            initial_theta=None,observer=None,centered_diagnostic=None,balanced_merit=False,target_moments=None,
            entropy_guard=None,residual_tolerance=8e-12):
    assert 0.<residual_tolerance<=8e-12
    U0,S0,_=model.moment(q)
    target=U0 if target_moments is None else np.asarray(target_moments)
    if target_moments is not None:
        assert vacuum_scaled and target.shape==U0.shape and np.all(np.isfinite(target))
    qref,R0,t=model.reference(q,guess)
    z=np.zeros_like(q) if reference_guess else q-qref
    if initial_theta is not None:
        assert not reference_guess and np.shape(initial_theta)==np.shape(q)
        z=np.asarray(initial_theta)-qref
    inactive=np.zeros(model.d//2,dtype=bool)
    if vacuum_scaled:
        touched=np.zeros(model.d//2,dtype=bool)
        touched[model.k.I]=True
        touched[model.n+model.k.J]=True
        missing=np.flatnonzero(~touched)
        if len(missing):
            for panel in getattr(model,"panels",()):
                touched[missing]|=np.any(panel[3][:,2*missing]>0.,axis=0)
        inactive=~touched
    active=np.repeat(~inactive,2)
    z[~active]=q[~active]-qref[~active]
    basis=np.eye(model.d)[active]
    w=np.repeat(U0[::2],2)
    for it in range(100 if vacuum_scaled else 22):
        U,S,C=model.moment(qref+z)
        F,jv=model.centered(z,qref,R0)
        if vacuum_scaled:
            # A newly populated bin need not solve to a tolerance relative
            # to its nearly zero *old* count. Normalize both moment defects
            # by the larger of that bin's old and returned count.
            w=np.repeat(np.maximum(U0[::2],U[::2]),2)
            w[~active]=1.
        r=((U-target-h*F)/w)[active]
        err=np.max(abs(r))
        if observer is not None: observer(it,float(err),int(np.flatnonzero(active)[np.argmax(abs(r))]))
        if err<residual_tolerance:
            if entropy_guard is None:
                assert S>=S0
            else:
                assert entropy_guard(qref+z,S,S0)
            if np.any(inactive):
                assert np.all(F[~active]==0.)
                assert np.all(abs(U[~active]-U0[~active])<=2e-12*np.abs(U0[~active]))
            if centered_diagnostic is not None:
                centered_diagnostic.update(reference=qref.copy(),centered=z.copy(),
                    parent_ratio=None if R0 is None else R0.copy())
            return qref+z,dict(newton_iterations=it,residual=float(err),entropy_gain_relative=float(S/S0-1),
                               maximum_exponent_departure_from_eq=float(np.max(abs(z[active]))),
                               structurally_inactive_cells=int(inactive.sum()),
                               structurally_inactive_vacuum_cells=int(np.sum(inactive & (U0[::2]==0.))))
        def der(d):
            a=d.reshape(-1,2)
            dn=-np.stack((C[:,0]*a[:,0]+C[:,1]*a[:,1],C[:,1]*a[:,0]+C[:,2]*a[:,1]),axis=1).ravel()
            return ((dn-h*jv(d))/w)[active]
        M=np.stack([der(d) for d in basis],axis=1)
        da=np.zeros(model.d)
        if balanced_merit:
            rs=1./np.maximum(np.max(abs(M),axis=1),1e-300)
            A=rs[:,None]*M
            cs=1./np.maximum(np.max(abs(A),axis=0),1e-300)
            da[active]=cs*np.linalg.solve(A*cs, -rs*r)
        else:
            da[active]=np.linalg.solve(M,-r)
        norm=np.linalg.norm(r)
        for j in range(28 if vacuum_scaled else 18):
            s=z+2.**-j*da
            try:UU,SS,_=model.moment(qref+s)
            except ValueError:continue
            trial_weights=w
            if balanced_merit and vacuum_scaled:
                trial_weights=np.repeat(np.maximum(U0[::2],UU[::2]),2)
                trial_weights[~active]=1.
            rr=((UU-target-h*model.centered(s,qref,R0)[0])/trial_weights)[active]
            if np.linalg.norm(rr)<norm*(1.-.001*2.**-j):break
        else:raise RuntimeError(f"centered residual not reduced: {err:.3g}, stepH={h:g}")
        z=s
    raise RuntimeError(f"centered iteration budget: {err:.3g}, stepH={h:g}")

def test(r,eos):
    T,H,k,N,t,Q,th=configuration(r,eos,7,7.)
    m=Centered(k,H)
    initial=np.concatenate([x.T.ravel() for x in th])
    initial+=.10*np.sin(.33*np.arange(len(initial)))+.03*np.cos(.74*np.arange(len(initial)))
    n0,S0,_=m.moment(initial)
    ref,R0,t=m.reference(initial,T)
    z=initial-ref
    f,j=m.centered(z,ref,R0)
    old,_=m.collision(initial)
    check0=np.linalg.norm(f-old)/np.linalg.norm(old)
    assert check0<1e-12 and np.linalg.norm(m.centered(z*0,ref,R0)[0])==0.
    dz=np.cos(np.arange(len(z))*.41)
    eps=4e-6
    x=[m.centered(z+i*eps*dz,ref,R0)[0] for i in (-1,1)]
    fd=np.linalg.norm((x[1]-x[0])/(2*eps)-j(dz))/np.linalg.norm(j(dz))
    assert fd<1e-6
    _,_,inv=metric(k,th)
    out=[]
    for h in (5e-9,2.5e-7,5e-5,.002,.02):
        q,log=advance(m,initial.copy(),h,T)
        ans,S,_=m.moment(q)
        e=max(abs((inv.T@ans)/(inv.T@n0)-1))
        assert e<1e-10
        vs=None
        if h<1e-6:
            a,_=previous(m,initial.copy(),h)
            vs=np.max(abs((ans-m.moment(a)[0])/np.repeat(n0[::2],2)))
            assert vs<3e-10
        out.append(dict(step_over_Hubble_time=h,invariant_defect=e,vs_previous_small_step=vs,**log))
    return dict(direct_flux_relative_check=check0,centered_JVP_finite_difference_error=fd,steps=out,
      scope="Fixed-background strictly populated conservative two-moment collision test, not a physical long-term cosmological result.")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("spectra")
    a=p.parse_args();eos=StandardModel(a.eos)
    r=next(r for r in json.loads(Path(a.spectra).read_text())["rows"] if r["TR"]==.05 and abs(r["TF"]-.045)<1e-8)
    print(json.dumps(test(r,eos),indent=2,allow_nan=False))
