"""Spin-summed distinguishable elastic current in a frozen isotropic Coulomb gauge."""
import argparse
import json

import numpy as np

from early_pole_frozen_dynamic_vlasov_reference import Snapshot, load, real_kernels


def vacuum_trace(s,t,ms):
    """Invariant two-current amplitude, stripped of g^4."""
    u_minus_ms = ms - s - t
    return 2*((s-ms)**2+u_minus_ms**2+2*t*ms)/(t*t)


def current_contractions(s,t,ea,eb,omega,masses):
    """Timelike, interference and transverse tensors projected along plasma q."""
    m1,m2=masses
    s,t,ea,eb,omega=np.broadcast_arrays(s,t,ea,eb,omega)
    Q2=omega*omega-t
    Q=np.sqrt(Q2)
    assert np.all(t<0.)
    a_par=(ea*omega+t/2)/Q
    b_par=(eb*omega-t/2)/Q
    a2=np.maximum(ea*ea-m1*m1-a_par*a_par,0.)
    b2=np.maximum(eb*eb-m2*m2-b_par*b_par,0.)
    ea1,eb1=ea+omega,eb-omega
    ab=ea*eb-.5*(s-m1*m1-m2*m2)-a_par*b_par
    timelike=(4*ea*ea1+t)*(4*eb*eb1+t)
    mixed=4*(ea+ea1)*(eb+eb1)*ab
    transverse=16*ab*ab-4*t*(a2+b2)+2*t*t
    return Q,timelike,mixed,transverse


def matrix_trace(s,t,ea,eb,omega,masses,polarization=None):
    Q,l,inter,trans=current_contractions(s,t,ea,eb,omega,masses)
    if polarization is None:
        d0=1/(Q*Q)
        dt=-1/t
    else:
        a,b,c,d=np.moveaxis(polarization(abs(omega)/Q),-1,0)
        d0=1/(Q**2+a+1j*c)
        dt=1/(-t+b-1j*d)
    terms=np.array([l*abs(d0)**2,trans*abs(dt)**2,-2*inter*np.real(d0*np.conj(dt))])
    total=np.sum(terms,axis=0)
    # A stable positive amplitude is expected, but cancellation is possible within
    # the floating point envelope for nearly collinear screened currents.
    assert np.all(total>=-2e-12*np.max(abs(terms),axis=0))
    return np.maximum(total,0.)


def htl_massless(ratio_debye_squared):
    def pol(x):
        a,b,inside,z=real_kernels(x,1.)
        return ratio_debye_squared*np.stack([a,b,np.pi/2*np.where(inside,z,0.),
                   np.pi/4*np.where(inside,z*(1-z*z),0.)],axis=-1)
    return pol


def tabulated_polarization(grid,g=.005,smear=.12):
    path=f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
    model=Snapshot(load(path),16,smear)
    # Resolve the large derivative near light cone with a denser nonuniform mesh.
    x=np.unique(np.r_[np.linspace(0.,.99,1000),1.-np.geomspace(1e-8,.01,1200)])
    v=g*g*model.retarded(x)
    return (lambda z: np.stack([np.interp(z,x,w) for w in v.T],axis=-1)),float(v[0,0])


def check():
    rng=np.random.default_rng(128)
    results=[]
    for masses in ((0.,0.),(.05,.003),(.32,.17)):
        m1,m2=masses
        a,b=rng.uniform(.4,8.,(2,3600))
        ea,eb=np.hypot(a,m1),np.hypot(b,m2)
        c=rng.uniform(-.98,.98,len(a))
        ms=m1*m1+m2*m2
        s=ms+2*(ea*eb-a*b*c)
        lam=(s-ms)**2-4*(m1*m2)**2
        boost=np.sqrt((ea+eb)**2-s)
        X=(s+m1*m1-m2*m2)/(2*s)
        Z=np.sqrt(lam)/(2*s)
        A=np.clip((ea-(ea+eb)*X)/(boost*Z),-1,1)
        theta=rng.uniform(-.96,.999,len(a))
        ang=rng.uniform(0.,np.pi,len(a))
        t=-lam/s*(1-theta)/2
        ea1=(ea+eb)*X+boost*Z*(A*theta+np.sqrt((1-A*A)*(1-theta*theta))*np.cos(ang))
        omega=ea1-ea
        v=vacuum_trace(s,t,ms)
        raw=matrix_trace(s,t,ea,eb,omega,masses)
        err=float(np.max(abs(raw/v-1.)))
        assert err<1e-10
        pol=htl_massless(.005**2/3)
        screened=matrix_trace(s,t,ea,eb,omega,masses,pol)
        reverse=matrix_trace(s,t,ea1,eb-omega,-omega,masses,pol)
        rev_err=float(np.max(np.abs(screened/reverse-1.)))
        assert rev_err<1e-10
        result=dict(masses=masses, vacuum_projection_error=err, detailed_reverse_kernel_error=rev_err)
        results.append(result)
    return dict(scope='Coulomb gauge tensor check and screened Born reference; no full kinetic rate.',cases=results)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.parse_args()
    print(json.dumps(check(),indent=2,allow_nan=False))
