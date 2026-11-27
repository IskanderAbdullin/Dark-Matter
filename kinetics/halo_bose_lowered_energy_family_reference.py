"""Virial energy and instantaneous constraints on a homologous lowered Bose family."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from halo_bose_finite_lowered_mixing_reference import solve
from halo_bose_finite_evaporation_kernel_reference import occupied


def moments(F,W,nx=260,nv=96):
    d,_,n,N,s,lr=solve(F,W,return_profile=True)
    x,w=roots_legendre(nx)
    y,wv=roots_legendre(nv)
    l=-15+(x+1)*(lr+15)/2
    r=np.exp(l)
    dxi=r**3*w*(lr+15)/2
    Y=np.clip(s(l)[0],0.,W)
    U2=(W-Y)[:,None]*((1+y)/2)**2
    f=occupied(Y[:,None],U2,d,W)
    n_local=np.sqrt(W-Y)*((U2*f)@wv)
    kin_local=np.sqrt(W-Y)*((U2**2*f)@wv)
    k=dxi@kin_local
    number=dxi@n_local
    assert abs(number/N-1)<7e-5
    Phi_rel0=Y-W-s(lr)[1]
    potential=.5*float(dxi@(Phi_rel0*n_local))
    assert abs(potential/(-2*k)-1)<6e-5
    return np.array([N/n**1.5,k/N]),float(abs(potential/(2*k)+1))


def run(ledger):
    inputs=json.loads(ledger.read_text())['cases']
    out=[]
    for a in inputs:
        F=a['peak_occupation'];W=a['escape_depth']
        ref,vd=moments(F,W)
        _,kap=ref
        J=a['infinity_flux_fraction_per_Gyr']
        E=a['infinity_mean_excess_energy_over_mV02']
        rhs=[-J,J*E/kap]
        dk=[]
        for fixed in ('F','W'):
            ds=[]
            ret=lambda h:moments(F,W+h)[0] if fixed=='F' else moments(F*np.exp(h),W)[0]
            for h in ((.06,.03,.015) if fixed=='F' else (.4,.2,.1)):
                g=(np.log(ret(h))-np.log(ret(-h)))/(2*h)
                A=np.array([[1.5,g[0]],[3.5,g.sum()]])
                ds.append((g,np.linalg.solve(A,rhs),np.linalg.cond(A)))
            final=ds[-1]
            dif=float(np.max(abs(ds[-2][1]/final[1]-1)))
            if fixed=='F':
                assert dif<.035,(W,fixed,ds)
            dk.append(dict(parameter_held_fixed=fixed,structural_log_derivatives=final[0].tolist(),
                           one_step_tangent_dlog_V0_and_dW_or_dlogF_per_Gyr=final[1].tolist(),
                           matrix_condition_number=final[2],
                           relative_difference_of_two_stencils=dif,
                           exploratory_ill_conditioned=bool(final[2]>1.e4 or dif>.035)))
        out.append(dict(peak_occupation=F,escape_depth=W,virial_K_over_MV02=kap,
                        independent_virial_relative_defect=vd,
                        physical_dlog_bound_mass_per_Gyr=-J,
                        physical_dlog_binding_magnitude_per_Gyr=J*E/kap,
                        fixed_shape_would_instead_dlog_binding_per_Gyr=-7*J/3,
                        conditional_tangents=dk))
    return dict(scope="Frozen isolated positive-energy escape and exact virial checks on a family of "
                      "stationary lowered halos. Tangents are algebraic hypothetical constraints on "
                      "shape changes, not isolated dynamical trajectories or a stability proof.",outputs=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('escape_ledger',type=Path)
    print(json.dumps(run(p.parse_args().escape_ledger),indent=2,allow_nan=False))
