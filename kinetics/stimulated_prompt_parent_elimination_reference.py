"""Conservative Bose-stimulated prompt-parent elimination at fixed fermion bath."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import expit
from factorized_prompt_return_graph_reference import graph
from tagged_scalar_nonthermal_bath_reference import Bath
from three_species_decay_onsager_reference import mp
from sm_dark_entropy_exchange_reference import StandardModel

def collision(F,weights,a,b, jac=False):
    assert np.all(F>=0)
    U=a*F
    V=b*(1+F)
    A=U.sum(axis=1)
    B=V.sum(axis=1)
    S=A+B
    Num=V*A[:,None]-U*B[:,None]
    source=(weights/S)@Num
    if not jac:return source,np.max(A/S)
    wd=(weights/S)[:,None]
    J=(wd*V).T@a -(wd*U).T@b - (wd*Num/S[:,None]).T@(a+b)
    J+=np.diag((weights/S)@(b*A[:,None]-a*B[:,None]))
    return source,J

def test(eos,source):
    rows=[]
    for T in (.018,.045,.09):
        edges=np.r_[0.,np.geomspace(.5*mp,5*T,21)]
        for eta_phi in (-.3,-1e-5,0.):
            eta_2=-2.8
            w,A,B,E=graph(edges,T=T,eta=eta_2,components=True,scalar_BE=(T,eta_phi))
            a=A.sum(axis=1);b=B.sum(axis=1)
            h=a/(a+b)
            analytic=expit(eta_phi+eta_2-E/T)
            nonzero=analytic>1e-200
            rel=np.max(abs(h[nonzero]/analytic[nonzero]-1))
            forward=A*(1-h[:,None])
            inverse=B*h[:,None]
            err=np.linalg.norm(forward-inverse)/np.linalg.norm(forward)
            assert max(err,rel)<2e-12
            rows.append(dict(T=T,eta_phi=eta_phi,
                   exact_thermal_kinetic_equilibrium_relative_defect=err,FD_parent_relative_error=rel))
    synthetic=[]
    for r in source["rows"]:
        if r["TR"] not in (.02,.05,.1) or abs(r["TF"]/r["TR"]-.9)>1e-6:continue
        T=r["TF"]
        edges=np.r_[0.,np.geomspace(.5*mp,9*T,23)]
        b0=Bath(r,eos)
        w,a,b=graph(edges,T=T,bath=b0,components=True)
        for occ in (.0001,.003,.1,1.,10.):
            n=a.shape[1]
            F=occ*(.08+np.sin(np.linspace(0.,4.,n))**2)/1.08
            C,pop=collision(F,w,a,b)
            C2,J=collision(F,w,a,b,True)
            lim=np.max(b0(b0.q))
            cap=occ*lim/(1-lim+occ*lim)
            assert pop<=cap*(1+2e-14)
            assert abs(C.sum())<2e-13*np.linalg.norm(C)
            assert np.linalg.norm(J.sum(axis=0))<5e-13*np.linalg.norm(J)
            rng=np.random.default_rng(10392)
            d=rng.normal(size=n);d*=.01*occ/max(abs(d))
            eps=1e-4
            plus,minus=(collision(F+i*eps*d,w,a,b)[0] for i in (1,-1))
            err=np.linalg.norm(J@d-(plus-minus)/(2*eps))/np.linalg.norm(J@d)
            assert err<2e-7,err
            X=(w/b.sum(axis=1))[:,None]
            G=(X*a).T@b
            linear=G.T@F-G.sum(axis=1)*F
            z=np.linalg.norm(C-linear)/np.linalg.norm(linear)
            if occ<=.003:assert z<.012
            cutoff=F.copy()
            idx=np.arange(2,n,4)
            cutoff[idx]=0.
            cut,_=collision(cutoff,w,a,b)
            assert np.min(cut[idx])>=0.
            synthetic.append(dict(T=T,controlled_scalar_f_peak=occ,nonthermal_parent_maximum=pop,
                     maximum_parent_occupation_analytic_bound=cap,
                     relative_linear_tag_defect=z,jacobian_refinement=err))
    return dict(thermal=rows,nonthermal_synthetic=synthetic,
                scope="Instantaneous prompt-parent elimination in a prescribed bath; no full stimulated time history or condition for a condensate.")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("spectra")
    a=p.parse_args()
    print(json.dumps(test(StandardModel(a.eos),json.loads(Path(a.spectra).read_text())),
            indent=2,allow_nan=False))
