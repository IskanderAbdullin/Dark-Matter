"""Independent quantum relative-entropy certificate of physical implicit moments.

For each *physical DOS quadrature point*, the entropy-dual gap is the
positive Bernoulli (FD) or geometric (BE) Kullback-Leibler divergence
D(old||new). At an approximate backward-Euler joint collision root, the
exact discrete identity is:

  S(new)-S(old)-dt H(new)
       = sum_species int DOS D(old||new) + eta(new) . residual.

No physical SM/FRW or scalar zero-mode is included. This is a discrete,
conditional entropy-certificate theorem, not unconditional positive-root
existence.
"""
import numpy as np
from scipy.special import expit

from kinetics.phi_bose_coupled_bernstein_p2_implicit import (
    asflat, physical_implicit_residual,
)
from kinetics.phi_bose_coupled_bernstein_p2_finite import Coupled


def scalar_quantum_kl(old_eta, new_eta, kind):
    """Relative entropy between Bernoulli or geometric quantum occupations.

    Canonical log-partition form is numerically safer than subtracting large
    mode-by-mode entropy densities, especially near a BE infrared endpoint.
    """
    old=np.asarray(old_eta,dtype=float)
    new=np.asarray(new_eta,dtype=float)
    if kind=="fermion":
        psi_old=np.logaddexp(0.,-old)
        psi_new=np.logaddexp(0.,-new)
        f_old=expit(-old)
    elif kind=="boson":
        if np.min(old)<=0 or np.min(new)<=0:
            raise ValueError("Bose eta <=0 in quantum relative entropy")
        psi_old=-np.log(-np.expm1(-old))
        psi_new=-np.log(-np.expm1(-new))
        f_old=1/np.expm1(old)
    else:
        raise ValueError("Unsupported quantum statistics")
    return (psi_new-psi_old)+(new-old)*f_old


def physical_relative_entropy(model,old_etas,new_etas):
    """Physical per-species integrated D(old||new), independently of events."""
    result=[]
    for species,old,new in zip(model.species,old_etas,new_etas):
        o=np.einsum("cqa,ca->cq",species.H,old.reshape((-1,3)))
        n=np.einsum("cqa,ca->cq",species.H,new.reshape((-1,3)))
        local=scalar_quantum_kl(o,n,species.kind)
        result.append(float(np.sum(species.weights*local)))
    return result


def entropy_residual_audit(model,old_etas,new_etas,dt):
    """Independent decomposition of entropy certificate, including root error."""
    old=model.integrate_spectrum(old_etas)
    new=model.integrate_spectrum(new_etas)
    source=model.source(new_etas)
    residual,meta=physical_implicit_residual(
        model,asflat(new_etas),asflat(old_etas),dt,with_jac=False)
    kl=physical_relative_entropy(model,old_etas,new_etas)
    lhs=(new["S"]-old["S"])-dt*source["H"]
    rhs=sum(kl)+float(np.dot(asflat(new_etas),residual))
    return dict(entropy_certificate=lhs,
                integrated_quantum_kl=kl,
                residual_work=float(np.dot(asflat(new_etas),residual)),
                reconstructed_certificate=rhs,
                decomposition_defect=lhs-rhs,
                actual_entropy_increase=new["S"]-old["S"],
                new_microscopic_H=source["H"])


def audit():
    from kinetics.phi_bose_coupled_bernstein_p2_implicit import implicit_step
    out=[]
    for P in (2.,8.,16.):
        m=Coupled(P=P,Nphi=12,nrad=6,nangle=12)
        e=m.initial_etas()
        s=implicit_step(m,e,1e5)
        a=entropy_residual_audit(m,e,s["etas"],1e5)
        out.append(dict(P=P,
                        iterations=s["iterations"],
                        scaled_residual=s["scaled_residual"],
                        decomposition=a,
                        step_min_bose_eta=s["min_bose_eta"]))
    return out


if __name__=="__main__":
    import json
    print(json.dumps(audit(),indent=2))
