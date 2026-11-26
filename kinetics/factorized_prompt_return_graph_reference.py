"""Parent-energy factorization for positive on-shell tagged scalar return."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre, expit
from scipy.linalg import eigvalsh
from decay_scalar_ir_relaxation_reference import limits
from prompt_resonant_scalar_redistribution_reference import matrix as angular_matrix
from three_species_decay_onsager_reference import M, m, mp, P, Es, width
from tagged_scalar_nonthermal_bath_reference import Bath
from sm_dark_entropy_exchange_reference import StandardModel, MPL

C = 4*width*M*M/(4*np.pi**2*P)

def shell_weights(edges,T=None):
    edges=np.asarray(edges)
    if T is None:
        return np.diff(edges**3)/(6*np.pi**2)
    z,w=roots_legendre(36)
    a=edges[:-1,None]; b=edges[1:,None]
    q=(a*(1-z)+b*(1+z))/2
    V=(b-a)*w*q*q*np.exp(-np.hypot(q,mp)/T)/(4*np.pi**2)
    return V.sum(axis=1)

def graph(edges, *, T, eta=0., bath=None, order=15, theta=16, max_parent=None, exact_FD=False,occupation=None,trial_T=None,components=False,scalar_BE=None):
    """Return G with rows per initial cell; last column is the unresolved exterior."""
    edges=np.asarray(edges)
    assert edges[0]==0 and np.all(np.diff(edges)>0)
    e=np.hypot(edges,mp)
    if max_parent is None:
        max_parent=max(M+70*T,e[-1]+m+48*T) if bath is None else e[-1]+bath.q[-1]
    assert max_parent>M
    _,opens,spans=map(np.asarray,zip(*(limits(p) for p in edges)))
    cut=np.unique(np.clip(np.r_[M,max_parent,opens,opens+spans,
                              M+T*np.geomspace(.00001,(max_parent-M)/T,41)],
                         M,max_parent))
    z,w=roots_legendre(order)
    lo=np.log(cut[:-1]); hi=np.log(cut[1:])
    E=np.exp((lo[:,None]*(1-z)+hi[:,None]*(1+z))/2).ravel()
    dE=((hi-lo)[:,None]*w/2).ravel()*E
    pp=np.sqrt(np.maximum(0.,(E-M)*(E+M)))
    ep=(E*Es+pp*P)/M
    em=(Es*Es+(pp*mp/M)**2)/ep
    low=np.maximum(em[:,None],e)
    upp=np.minimum(ep[:,None],np.r_[e[1:],np.inf])
    span=np.maximum(upp-low,0.)
    good=span.sum(axis=1)>0
    low,span,E,dE=(a[good] for a in (low,span,E,dE))
    u,v=roots_legendre(theta)
    ef=low[:,:,None]+span[:,:,None]*(1+u)/2
    Ec=E[:,None,None]-ef
    if bath is None:
        ap=expit(eta-Ec/T)
        om=expit(Ec/T-eta)
    else:
        ap=bath(Ec) if occupation is None else occupation(Ec)
        om=1-ap
    weights=span[:,:,None]*v/2
    released=np.sum(weights*om,axis=-1)
    denom=released.sum(axis=1)
    assert np.all(denom>0)
    trial=1. if trial_T is None else np.exp(-ef/trial_T)
    captured_flux=np.sum(weights*ap*trial,axis=-1)
    if components:
        if scalar_BE is None:
            return C*dE,captured_flux,released
        t,mu=scalar_BE
        assert mu<mp/t
        f=1/np.expm1(ef/t-mu)
        return C*dE,np.sum(weights*ap*f,axis=-1),np.sum(weights*om*(1+f),axis=-1),E
    raw=(captured_flux[:,:-1].T*(C*dE/denom))@released
    if exact_FD:
        assert bath is None
        # Analytic KMS reduction; comparison against direct two Fermi factors
        # remains separate and is not fixed by averaging measured entries.
        tangent=released.T*(C*dE*np.exp(eta-E/T)/denom)
        K=tangent@released
        actual=np.sum(weights*ap*np.exp(-ef/T),axis=-1)
        cmp=(actual.T*(C*dE/denom))@released
        return raw,K,cmp
    return raw

def tests(eos,sources):
    out=[]
    for T in (.018,.045,.09):
        ed=np.r_[0.,np.geomspace(.5*mp,4*T,23)]
        raw,G,B=graph(ed,T=T,eta=-3.,exact_FD=True)
        V=shell_weights(ed,T)
        A=B[:-1,:-1]
        res=np.linalg.norm(A-A.T)/np.linalg.norm(G[:-1,:-1])
        eq=np.linalg.norm(G-B)/np.linalg.norm(G)
        L=np.diag(A.sum(axis=1))-A
        S=L/np.sqrt(V[:,None]*V)
        lam=eigvalsh(S)
        evec=np.sqrt(V)
        assert min(lam)>-1e-10*max(lam) and np.linalg.norm(S@evec)<1e-13*np.linalg.norm(S)*np.linalg.norm(evec)
        assert max(res,eq)<2e-12
        # Uniform occupation inside each bin is not an exact equilibrium
        # closure. Keep its measured error separate from the KMS quadrature.
        L0=np.diag(raw[:,:-1].sum(axis=1))-raw[:,:-1].T
        bias=np.linalg.norm(L0@(V/shell_weights(ed)))/np.linalg.norm(raw[:,:-1]@(V/shell_weights(ed)))
        out.append(dict(T=T,reciprocity=res,factorization_kms_error=eq,
                   flat_trial_stationary_defect=bias))
    comparison=[]
    for r in sources["rows"]:
        if r["TR"] not in (.02,.05,.1) or abs(r["TF"]/r["TR"]-.9)>1e-6: continue
        T=r["TF"]
        b=Bath(r,eos)
        ed=np.r_[0.,np.geomspace(.5*mp,2.4*T,17)]
        g=graph(ed,T=T,bath=b)
        f=graph(ed,T=T,bath=b,order=20,theta=24)
        a,counts=angular_matrix(b,ed,npin=18,n=10)
        V=shell_weights(ed)
        active=f.sum(axis=1)>max(f.sum(axis=1))*1e-7
        delta=np.sum(abs(g-f),axis=1)[active]/f.sum(axis=1)[active]
        diff0=np.linalg.norm(f/V[:,None]-a,ord=1)/np.linalg.norm(a,ord=1)
        num=float(np.linalg.norm(f.sum(axis=1)/V-counts)/np.linalg.norm(counts))
        assert np.max(delta)<.035 and num<.12
        comparison.append(dict(TR=r["TR"],T=T,refinement_max_significant_row_l1=float(max(delta)),
             old_angular_vs_new_matrix_l1=float(diff0),total_rate_difference=num,
             H=float(np.sqrt(eos.state(T)[0]/3)/MPL)))
    return dict(thermal=out,nonthermal=comparison)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("source")
    a=p.parse_args()
    print(json.dumps(tests(StandardModel(a.eos),json.loads(Path(a.source).read_text())),
                     indent=2,allow_nan=False))
