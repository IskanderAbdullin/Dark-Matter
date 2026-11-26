"""Entropy-consistent *instantaneous* χℓ 2->2 binned Pauli weak-to-projected rates.

The previous linear barycentric energy deposition conserves χ/ℓ N
and energy but fails the H theorem on nonlinear Fermi log odds if
the event f(E) is evaluated analytically and source entropic logs at
grid nodes. Repair the INCONSISTENCY, not the negative result:

  h_chi(E) = sum_i alpha_i(E) logit(f_chi_i)
  f_event(E) = expit(h_chi(E))     (same for lepton)

Use the *same* nonnegative alpha_i in event evaluation and source
deposition. Then Sdot_bins=-sum(dN_i*logit(f_i)) is IDENTICALLY
the positive on-shell quantum Sdot_event=sum W(F-R)ln(F/R).
C/L/energy remain exactly paired at this quadrature.

LIMIT: this is only a finite-grid, prescribed-occupation instantaneous
projected 5D collision. It is NOT a proven Pauli-bound-preserving time
solver, full χ1/χ2-SM evolution or continuum converged physical Γel.
Coarse logit grids severely alias oscillatory Fermi distributions.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expit,logit
from kinetics.paired_chi_lepton_barycentric_5d import projected_pair


def logit_interpolant(energy_nodes,occupation_nodes):
    nodes=np.asarray(energy_nodes,float)
    occ=np.asarray(occupation_nodes,float)
    if (nodes.ndim!=1 or occ.shape!=nodes.shape or len(nodes)<4
        or not np.all(np.isfinite(nodes)) or
        not np.all(np.diff(nodes)>0) or
        not np.all(np.isfinite(occ)) or
        np.any(occ<=0) or np.any(occ>=1)):
        raise ValueError("strict interior Fermi occupation on ordered energy grid")
    odds=logit(occ)
    def f(energy):
        values=np.asarray(energy,float)
        if (not np.all(np.isfinite(values)) or
            np.any(values<nodes[0]-1e-10) or
            np.any(values>nodes[-1]+1e-10)):
            raise ValueError("energy outside physical logit interpolation grid")
        return expit(np.interp(values,nodes,odds))
    return f


def entropy_consistent_projection(fchi,flepton,*,M=.1,m=0.,Te=1.,
                                  cut_kinetic_T=12.,bins=20,
                                  G=.001,mV=.1,
                                  nr=4,nmu=4,nphi=12,nang=12):
    """Evaluate same microscopic F/R and deposition with logit f.

    Original fchi/flepton are arbitrary *prescribed* legal isotropic
    functions used to sample grid occupations. The new continuous
    functions are reconstructed by affine interpolation of node LOGITS.
    """
    if not (isinstance(bins,int) and bins>=4 and M>0 and m>=0 and
            Te>0 and cut_kinetic_T>0):
        raise ValueError
    upper=M+m+2*cut_kinetic_T*Te
    chi_nodes=np.linspace(M,upper,bins+1)
    lep_nodes=np.linspace(m,upper,bins+1)
    chi_occ=np.asarray(fchi(chi_nodes),float)
    lep_occ=np.asarray(flepton(lep_nodes),float)
    recon_chi=logit_interpolant(chi_nodes,chi_occ)
    recon_lep=logit_interpolant(lep_nodes,lep_occ)
    result=projected_pair(recon_chi,recon_lep,M=M,m=m,Te=Te,G=G,mV=mV,
                          nr=nr,nmu=nmu,nphi=nphi,nang=nang,
                          cut_kinetic_T=cut_kinetic_T,bins=bins)
    result["logit_consistent_with_deposition"]=True
    result["entropy_identity_abs_residual"]=abs(
        result["interpolated_bin_entropy"]-result["direct_event_entropy"])
    result["not_a_Pauli_time_integrator_or_physical_gamma_el"]=True
    return result


if __name__=="__main__":
    import json
    def pack(r):
        return {k:v for k,v in r.items() if not isinstance(v,np.ndarray)}
    oscchi=lambda E:expit(-np.asarray(E)+.3*np.sin(4*np.asarray(E)))
    osclep=lambda E:expit(-np.asarray(E)+.3*np.cos(4*np.asarray(E)))
    rows={}
    for n in (10,20,40,80,160,320):
        rows[str(n)]=pack(entropy_consistent_projection(oscchi,osclep,bins=n))
    print(json.dumps({"oscillatory_nonthermal":rows,"class_ready":False},indent=2))
