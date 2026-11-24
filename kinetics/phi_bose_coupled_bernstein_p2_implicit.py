"""Fully implicit simultaneous physical-capacity Yukawa quantum collision moments.

Reuses the *same positive shared on-shell event weights* and physical
g*p²dp Bernstein-P2 moments for χ1, χ2, φ as the coupled explicit prototype.
The implicit moment equation M(ηnew)-M(ηold)=dt*C(ηnew), if solved inside
positive Bose/Pauli domain, obeys entropy concavity:
Snew-Sold >= ηnew·(Mnew-Mold) = dt*Hnew >= 0.

Numerical Newton convergence and positive root existence are conditional,
not unconditional theorems or physical late cosmological DM predictions.
"""
import math
import numpy as np
from scipy.sparse import block_diag,hstack
from scipy.sparse.linalg import spsolve
from scipy.special import expit
from kinetics.phi_bose_coupled_bernstein_p2_finite import Coupled


def asflat(vecs):
    return np.concatenate(vecs)


def unflat(model,eta):
    return [np.asarray(eta[a:b]) for a,b in zip(model.offsets[:-1],model.offsets[1:])]


def physical_implicit_residual(model,eta,old,dt,with_jac=True):
    """F=Mnew-Mold-dt*Cnew, with analytic sparse Jacobian of all quantum J."""
    es=unflat(model,eta)
    os=unflat(model,old)
    changes=[s.moment_difference(e,o) for s,e,o in zip(model.species,es,os)]
    e1,e2,ep=[A@e for A,e in zip(model.A,es)]
    if np.min(ep)<=0:raise ValueError('Bose eta must remain strictly positive')
    f=expit(-e1);g=expit(-e2);b=1/np.expm1(ep)
    F=f*(1-g)*(1+b)
    R=(1-f)*g*b
    J=F-R
    S=hstack([-model.A[0],model.A[1],model.A[2]],format='csr')
    weighted=model.W*J
    C=np.asarray(S.T@weighted).ravel()
    residual=asflat(changes)-dt*C
    H=float(np.dot(weighted,S@eta))
    if not with_jac:
        return residual,dict(H=H,source=C)
    # dJ/deta1 = -(1-g+b)*f*(1-f)
    # dJ/deta2 =  (f+b)*g*(1-g)
    # dJ/detaphi = (g-f)*b*(1+b)
    dj1=-f*(1-f)*(1-g+b)
    dj2=g*(1-g)*(f+b)
    djp=b*(1+b)*(g-f)
    dJ=hstack([A.multiply((model.W*dj)[:,None])
               for A,dj in zip(model.A,(dj1,dj2,djp))],format='csr')
    dC=S.T@dJ
    G=block_diag([mat for s,e in zip(model.species,es)
                       for mat in s.evaluate(e,with_G=True)['G']],format='csr')
    jac=-G-dt*dC
    return residual,jac,dict(H=H,source=C,weights_positive=bool(min(model.W)>0))


def implicit_step(model,old_etas,dt,maxiter=24,rtol=5e-8):
    """Damped full-coupled sparse Newton solving all physical moment equations.

    A *valid* root satisfies finite-step entropy inequality by concavity,
    but no unconditional numerical positive root existence is claimed.
    """
    assert dt>0
    old=asflat(old_etas)
    initial=model.integrate_spectrum(old_etas)
    # At exact common T, mu1=mu2+muphi, each event has zero affinity.
    # Ignore only event-affinity roundoff at machine scale (not a physical
    # nonzero chemical affinity), to prevent false equilibrium Newton drift.
    affinity=-model.A[0]@old_etas[0]+model.A[1]@old_etas[1]+model.A[2]@old_etas[2]
    machine_affinity=64*np.finfo(float).eps*max(1.,max(abs(model.A[0]@old_etas[0])))
    if max(abs(affinity))<machine_affinity:
        return dict(etas=old_etas,accepted_dt=dt,iterations=0,
                    scaled_residual=0.,stalled_at_precision=True,
                    entropy_increase=0.,entropy_certificate=0.,
                    number_residuals=(0.,0.),energy_defect=0.,
                    min_bose_eta=float(min(old_etas[2])))
    originalM=[s.evaluate(e,with_G=False)['M'] for s,e in zip(model.species,old_etas)]
    baseline=asflat(originalM)
    oldsource=asflat(model.source(old_etas)['C'])
    scale=np.maximum.reduce([np.abs(dt*oldsource),1e-14*np.abs(baseline),
                             np.full(len(baseline),1e-20)])
    if np.max(np.abs(dt*oldsource)/scale)<1e-9:
        return dict(etas=old_etas,accepted_dt=dt,iterations=0,
                    scaled_residual=0.,stalled_at_precision=True,
                    entropy_increase=0.,entropy_certificate=0.,
                    number_residuals=(0.,0.),energy_defect=0.,
                    min_bose_eta=float(min(old_etas[2])))
    new=old.copy()
    errs=[]
    for i in range(maxiter):
        res,jac,info=physical_implicit_residual(model,new,old,dt)
        err=max(abs(res/scale))
        errs.append(err)
        if err<rtol:break
        delta=spsolve(jac,-res)
        if not np.all(np.isfinite(delta)):raise ArithmeticError('Non-finite Newton update')
        alpha=1.
        accepted=False
        for j in range(28):
            candidate=new+alpha*delta
            if min(candidate[model.offsets[2]:])<=0:
                alpha/=2;continue
            try:rr,_=physical_implicit_residual(model,candidate,old,dt,with_jac=False)
            except (OverflowError,ValueError,FloatingPointError):alpha/=2;continue
            enew=max(abs(rr/scale))
            if enew<=err*(1-0.0001*alpha) or enew<rtol:
                accepted=True;new=candidate;break
            alpha/=2
        if not accepted:
            raise ArithmeticError(('Implicit Newton line search failed; no nonexistence theorem',dt,i,err))
    else:
        raise ArithmeticError(('Implicit Newton not converged',dt,errs[-1]))
    newetas=unflat(model,new)
    after=model.integrate_spectrum(newetas)
    newsource=model.source(newetas)
    energy_defect=sum(a-b for a,b in zip(after['Q'],initial['Q']))
    number_residuals=((after['N'][0]-initial['N'][0])+(after['N'][1]-initial['N'][1]),
                      (after['N'][0]-initial['N'][0])+(after['N'][2]-initial['N'][2]))
    sdiff=after['S']-initial['S']
    cert=sdiff-dt*newsource['H']
    return dict(etas=newetas,accepted_dt=dt,iterations=i,
                scaled_residual=errs[-1],residual_history=errs,
                entropy_increase=sdiff,entropy_certificate=cert,
                microscopic_H_new=newsource['H'],
                number_residuals=number_residuals,energy_defect=energy_defect,
                min_bose_eta=float(min(newetas[2])))


def implicit_adaptive_step(model,etas,dt,max_halvings=14):
    """Numerically reject failed roots, conservation, and entropy violations."""
    attempts=[]
    for _ in range(max_halvings+1):
        try:
            r=implicit_step(model,etas,dt)
            s=max(1e-25,dt*abs(model.source(etas)['Rphi']))
            if r['entropy_increase']<-1e-10*s:
                raise ArithmeticError('Physical entropy decreased')
            if r['entropy_certificate']<-1e-9*s:
                raise ArithmeticError('Entropy concavity certificate negative')
            if max(abs(r['energy_defect']),*(abs(v) for v in r['number_residuals']))>1e-8*s:
                raise ArithmeticError('Conserved charges/energy defect')
            if r['min_bose_eta']<=0:raise ArithmeticError('Bose eta nonpositive')
            r['rejections']=attempts
            return r
        except (ArithmeticError,ValueError,np.linalg.LinAlgError,OverflowError) as e:
            attempts.append((dt,str(e)))
            dt/=2
    raise ArithmeticError(('All implicit steps rejected',attempts))


def integrate_to(model,initial_etas,end_time,dt_max):
    assert end_time>0 and dt_max>0
    eta=initial_etas
    t=0.;steps=[]
    while t<end_time-1e-9:
        res=implicit_adaptive_step(model,eta,min(dt_max,end_time-t))
        eta=res['etas'];t+=res['accepted_dt'];steps.append(res)
        if len(steps)>1000:raise RuntimeError('Too many implicit steps')
    return dict(etas=eta,elapsed=t,steps=len(steps),
                rejected=sum(len(s['rejections']) for s in steps),history=steps)


def audit():
    for P in (2.,8.,16.):
        m=Coupled(P=P,Nphi=16,nrad=8,nangle=12)
        e=m.initial_etas()
        for dt in (1e5,3e5):
            v=implicit_adaptive_step(m,e,dt)
            print('P',P,'requested dt',dt,'accepted dt',v['accepted_dt'],
                  'iterations',v['iterations'],'relative residual',v['scaled_residual'],
                  'entropy gain',v['entropy_increase'],
                  'entropy certificate',v['entropy_certificate'],
                  'Q defect',v['energy_defect'],'charge defects',v['number_residuals'],
                  'rejected trials',v['rejections'])

if __name__=='__main__':
    audit()
