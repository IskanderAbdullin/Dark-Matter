"""A number-conserved elastic control contrasting full chemical equilibrium."""
import argparse
import json
from pathlib import Path

import numpy as np
from charged_lepton_portal_local_longitudinal_reference import model
from positive_fd_portal_quadrature_reference import PositivePortal


def project(T,xi,p,r,proxy):
    def m(k):
        return model(T,xi,p,r,k,proxy,elastic_only=True)
    A,H,_,_,_,locked,_=m(0.)
    hs,V=np.linalg.eigh(H)
    W=(V*np.sqrt(hs))@V.T
    J=(V/np.sqrt(hs))@V.T
    C=(-(W@A@J)).real
    B=(-(W@(m(1.)[0]-A)@J).imag)
    rates,U=np.linalg.eigh(.5*(C+C.T))
    if np.max(abs(rates[:3]))>2e-13 or min(rates[3:])<=1e-10:
        raise ValueError("Frozen pair number, energy, or momentum was not conserved.")
    N,F=U[:,:3],U[:,3:]
    speeds,V=np.linalg.eigh(N.T@B@N)
    if abs(speeds[1])>2e-13 or abs(speeds[-1]**2-locked)>2e-13:
        raise ValueError("Frozen-number common sound did not match thermodynamic slope.")
    out=[]
    for i,label in ((1,"nonpropagating"),(2,"sound")):
        v=N@V[:,i]
        b=B@v
        g=(F.T@b)/np.sqrt(rates[3:])
        a=float(g@g)
        u=b.copy(); u[:3]=0.
        z=(F.T@u)/np.sqrt(rates[3:])
        bulk=float(z@z)
        k=.01*min(rates[3:])
        w,ve=np.linalg.eig(W@m(k)[0]@J)
        j=int(np.argmin(abs(w+1j*k*speeds[i])))
        z=ve[:,j]
        loss=float(rates[3:]@abs(F.T@z)**2/np.vdot(z,z).real/k**2)
        if abs(loss/a-1)>4e-3 or abs((-w[j].imag)/k-speeds[i])>2e-4:
            raise ValueError(f"Direct small-k frozen coefficient did not match projection: {T,xi,label,loss,a}.")
        out.append(dict(branch=label,speed_squared=float(speeds[i]**2),damping=a,
                       bulk_term=bulk,heat_term=a-bulk,smaller_k=k,direct_damping=loss))
    return out


def run(receipt):
    p=PositivePortal((13,36,26),tail=65.)
    out=[]
    for r in receipt['rows']:
        for xi in (1.,10.):
            T=r['T']
            a,b=[project(T,xi,p,r,s) for s in (None,"heat_as_momentum")]
            expected=r['spacelike_bulk_over_heat']
            for x,y in zip(a,b,strict=True):
                if (abs(y['bulk_term']/x['bulk_term']-expected)>2e-10
                        or abs(x['heat_term']-y['heat_term'])>1e-11):
                    raise ValueError("Replacing drag failed to isolate the relative-velocity sector.")
            out.append(dict(T=T,visible_to_dark_enthalpy_ratio=xi,true_drag=a,heat_proxy=b))
    return dict(scope="Static ideal thermal limit with both directions of pair reactions removed; not a cosmological transition.",
                rows=out)


if __name__=='__main__':
    z=argparse.ArgumentParser(description=__doc__)
    z.add_argument('combined_portal_receipt')
    a=z.parse_args()
    print(json.dumps(run(json.loads(Path(a.combined_portal_receipt).read_text())),indent=2,allow_nan=False))
