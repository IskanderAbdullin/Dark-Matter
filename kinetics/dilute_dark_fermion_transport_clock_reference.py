"""Same-charge vector transport clock for an injected dark Dirac spectrum."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from three_species_decay_onsager_reference import m

G=.1
MV=10.

def amplitude(s,c):
    t=-.5*(s-4*m*m)*(1-c)
    u=-.5*(s-4*m*m)*(1+c)
    A=s-2*m*m
    d=t-MV*MV
    e=u-MV*MV
    return 2*G**4*((A*A+(u-2*m*m)**2+4*m*m*t)/d**2+
             (A*A+(t-2*m*m)**2+4*m*m*u)/e**2+
             2*A*(s-6*m*m)/(d*e))

def direct_spinors(s,c):
    I=np.eye(2)
    O=np.zeros((2,2))
    S=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]])
    gamma=np.array([np.block([[I,O],[O,-I]])]+[np.block([[O,a],[-a,O]]) for a in S])
    sig=np.array([1.,-1.,-1.,-1.])
    E=np.sqrt(s)/2
    p=np.sqrt((E-m)*(E+m))
    z=np.array([0.,0.,p])
    w=p*np.array([np.sqrt(1-c*c),0.,c])
    def U(vec):
        return np.sqrt(E+m)*np.vstack([I,np.tensordot(vec,S,axes=(0,0))/(E+m)])
    mats=list(map(U,(z,-z,w,-w)))
    def current(i,j):
        return np.array([mats[i].conj().T@gamma[0]@g@mats[j] for g in gamma])
    t=-2*p*p*(1-c)
    u=-2*p*p*(1+c)
    A=np.einsum('j,jac,jbd->abcd',sig,current(2,0),current(3,1))/(t-MV*MV)
    B=np.einsum('j,jbc,jad->abcd',sig,current(3,0),current(2,1))/(u-MV*MV)
    return np.sum(abs(A-B)**2)*G**4/4

def tailavg(s,deg):
    z,w=roots_legendre(deg)
    # P2 weight is invariant under exchange of the two identical finals.
    sig2=sum(weight*1.5*(1-u*u)*amplitude(s,u) for u,weight in zip(z,w))/(64*np.pi*s)
    return sig2

def estimate(hist,edges,T,angle=24):
    # The first cell extends to p=0; a finite representative suffices for
    # this dilute integral; do not use this to infer the occupancy at p=0.
    edges=np.asarray(edges)
    p=np.sqrt(edges[:-1]*edges[1:])*T
    p[0]=edges[1]*T*.5
    e=np.hypot(p,m)
    w=np.array(hist)
    vol=T**3*np.diff(edges**3)/(6*np.pi**2)
    x,wg=roots_legendre(angle)
    invariant=e[:,None,None]*e[None,:,None]-p[:,None,None]*p[None,:,None]*x
    s=2*(m*m+invariant)
    v=np.sqrt(np.maximum(0.,(invariant-m*m)*(invariant+m*m)))/(e[:,None,None]*e[None,:,None])
    Q=tailavg(s,20)
    return np.sum(w[:,None]*w[None,:]*np.sum(v*Q*wg/2,axis=2))/w.sum()**2,vol

def run(file,source):
    data=json.loads(Path(source).read_text())
    eos=StandardModel(file)
    errors=[]
    for S in (.00160001,.003,.03,.4,4.):
        for z in (-.9,-.6,0.,.36,.89):
            a=amplitude(S,z)
            b=direct_spinors(S,z)
            errors.append(abs(a/b-1))
    assert max(errors)<1e-10,max(errors)
    out=[]
    for row in data["rows"]:
        T=row["TF"]
        rs,_,ss,cs=eos.state(T)
        H=np.sqrt(rs/3)/MPL
        Y=sum(row["stable_yields"])
        n=Y*ss
        hist=row["stable_normalized_counts"]
        edges=row["scalar_differential_bins"]
        rate,vol=estimate(hist,edges,T)
        fine,_=estimate(hist,edges,T,38)
        er=abs(rate/fine-1)
        assert er<2e-6,er
        f=n*np.array(hist)/(4*vol)
        out.append(dict(TR=row["TR"],TF=T,identical_charge_l2_clock_over_H=(n*rate/(2*H) if f.max()<.08 else None),
                        dilute_bin_occupancy_gate_passed=bool(f.max()<.08),
                        shell_average_max_FD_occupation=f.max(),
                        maximum_reportable_relative_momentum_bin_step=max(np.diff(np.log(edges[1:]))),
                        angular_quadrature_relative_error=er))
    return dict(rows=out,spinor_crosscheck=max(errors),
                scope="Dilute same-charge vector elastic CM quadrupole transport proxy, not the thermalization eigenvalue of a bath with actual nonthermal blocking.")

if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument('eos');ap.add_argument('spectra')
    a=ap.parse_args()
    print(json.dumps(run(a.eos,a.spectra),indent=2,allow_nan=False))
