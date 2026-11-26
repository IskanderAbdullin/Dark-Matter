"""Fixed-T variational H-theorem check below the BE capacity."""
import json
import numpy as np
from scipy.special import expit, xlogy
from true_event_fv_evolution_reference import energy_mesh, Reconstruction


def check(n,M):
    rng=np.random.default_rng(8220944+n)
    q,m,b,(i,j,eb,ef,w)=energy_mesh(n,M,(.08,.18),8)
    B,F=Reconstruction(q,m,True),Reconstruction(q,b,False)
    mb=-rng.uniform(.003,.12,n)
    mf=rng.uniform(-.14,.14,n)
    nb=B.number*B.average(mb-m)
    nf=F.number*F.average(mf-b)
    tb,sb=B.solve(nb)
    tf,sf=F.solve(nf)
    assert np.all(tb<B.cap)
    phi=B.at(tb,sb,eb,i,"sat-high")
    ferm=F.at(tf,sf,ef,j,"sat-high")
    parent=expit(-eb-ef)
    g=w*parent*(1-ferm)*(1+phi)
    loss=w*(1-parent)*ferm*phi
    flux=g-loss
    # Form the logarithm from the physical populations, not the net-rate sign.
    logratio=np.log(g)-np.log(loss)
    energy_relative=float(np.max(abs(logratio+(mb[i]+mf[j]))))
    assert energy_relative<2e-10
    s=float(flux@logratio)
    assert s>0 and np.min(flux*logratio)>-1e-25
    dnb=np.bincount(i,weights=flux,minlength=n)
    dnf=np.bincount(j,weights=flux,minlength=n)
    derivative=float(mb@dnb+mf@dnf)
    # Thermodynamic free-energy derivatives from independent +/- perturbations.
    def parts(model,N):
        av,shift=model.solve(N)
        if model.bose and np.any(av>=model.cap):
            raise ValueError("thermal shape capacity exceeded")
        E=model.k+model.mass
        if model.bose:
            f=1/np.expm1(model.k-shift[:,None])
            entropy=np.log1p(f)+xlogy(f,1+1/f)
        else:
            f=expit(shift[:,None]-model.k)
            entropy=-xlogy(f,f)-xlogy(1-f,1-f)
        return model.number*(model.w*(E*f-entropy)).sum(1)
    h=1e-5
    dp=parts(B,nb+h*dnb)-parts(B,nb-h*dnb)
    df=parts(F,nf+h*dnf)-parts(F,nf-h*dnf)
    fd=float((dp.sum()+df.sum())/(2*h))
    assert abs(derivative+s)/s<1e-9 and abs(fd+s)/s<1e-5
    return dict(M=M,n=n,gain_events=int(np.count_nonzero(flux>0)),
                loss_events=int(np.count_nonzero(flux<0)),
                entropy_production=s,free_energy_derivative=derivative,
                finite_difference_free_energy_derivative=fd,
                max_log_affinity_identity_error=energy_relative,
                relative_finite_difference_residual=abs(fd+s)/s)

if __name__=="__main__":
    print(json.dumps({"cases":[check(n,m) for n in (24,48) for m in (.1,.005)],
                      "status":"fixed-T thermal domain only; no expanding-background H theorem"},
                     indent=2,allow_nan=False))
