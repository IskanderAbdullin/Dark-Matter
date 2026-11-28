"""Imposed chemistry clock sweep at fixed thermodynamic background."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp, quad
from charged_lepton_portal_local_longitudinal_reference import model
from positive_fd_portal_quadrature_reference import PositivePortal


def run(r):
    p=PositivePortal((13,36,26),tail=65.)
    T=r['T']
    xi=10.
    L,H,*_=model(T,xi,p,r,0.)
    E=model(T,xi,p,r,0.,elastic_only=True)[0]
    L1=model(T,xi,p,r,1.,elastic_only=True)[0]
    a,V=np.linalg.eigh(H)
    W=(V*np.sqrt(a))@V.T
    J=(V/np.sqrt(a))@V.T
    norm=lambda X:W@X@J
    C=-norm(E).real
    B=-norm(L1-E).imag
    P=-(1+r['timelike_heat']/r['spacelike_heat'])*norm(L).real-C
    R,U=np.linalg.eigh((C+C.T)/2)
    S,F=U[:,:3],U[:,3:]
    lam=np.linalg.eigvalsh(S.T@P@S)[-1]
    k=0.00025
    start,end=100*k/lam,.01*k/lam
    def A(eps):
        return -C-eps*P-1j*k*B
    def reduce(M):
        U=S-F@np.linalg.solve(F.T@M@F,F.T@M@S)
        return S.T@M@U,U
    def sound(M):
        eig,vec=np.linalg.eig(M)
        idx=int(np.argmin(eig.imag))
        return eig[idx],vec[:,idx]/np.linalg.norm(vec[:,idx])
    _,initial=sound(A(start))
    lam_final,xf=np.linalg.eig(A(end))
    j=int(np.argmin(lam_final.imag))
    rows=[]
    for sweep in (.03,.3,3.,30.):
        finish=np.log(start/end)*sweep
        eps=lambda x:start*np.exp(-x/sweep)
        M=lambda x:A(eps(x))/k
        shifted=lambda x:reduce(M(x))[0]
        vals=np.linspace(0.,finish,131)
        opts=dict(t_eval=vals,rtol=2e-12,atol=2e-14,method='BDF',dense_output=True)
        u=solve_ivp(lambda t,z:M(t)@z,(0.,finish),initial,jac=lambda t,z:M(t),**opts)
        v=solve_ivp(lambda t,z:shifted(t)@z,(0.,finish),S.T@initial,jac=lambda t,z:shifted(t),**opts)
        if not(u.success and v.success):
            raise ValueError("Sweep integration did not complete.")
        n=np.sum(abs(u.y)**2,axis=0)
        if np.max(np.diff(n))>2e-8:
            raise ValueError("Positive collision generator created perturbation free energy.")
        def loss(x):
            z=u.sol(x)
            return float((2*np.vdot(z,(C+eps(x)*P)@z).real/k))
        integral=quad(loss,0.,finish,epsabs=1e-10,epsrel=1e-8,limit=130)[0]
        free_energy_error=abs((n[0]-n[-1])-integral)
        if free_energy_error>2e-7:
            raise ValueError(f"Integrated entropy-loss and endpoint norm did not agree: {T}, {sweep}, {free_energy_error}.")
        exp_damping=quad(lambda x:sound(M(x))[0].real,0.,finish,epsabs=1e-11,limit=120)[0]
        def realblock(x):
            a=M(x)
            return np.block([[a.real,-a.imag],[a.imag,a.real]])
        ref=solve_ivp(lambda x,z:realblock(x)@z,(0.,finish),np.r_[initial.real,initial.imag],
                      method='Radau',jac=lambda x,z:realblock(x),rtol=2e-11,atol=2e-13)
        if not ref.success:
            raise ValueError("Independent real-block integration failed.")
        expected=ref.y[:5,-1]+1j*ref.y[5:,-1]
        if np.linalg.norm(u.y[:,-1]-expected)>4e-7:
            raise ValueError("Two stiff time integrators disagree.")
        final=u.y[:,-1]
        corrected=reduce(A(end))[1]@v.y[:,-1]
        decomp=np.linalg.solve(xf,final)
        rows.append(dict(sweep_duration_in_phase_units=sweep,projected_chemical_rate_to_wave_number_ratio=[100,.01],
                         endpoint_free_energy_balance_error=free_energy_error,
                         radau_bdf_terminal_error=float(np.linalg.norm(u.y[:,-1]-expected)),
                         norm_full=float(np.linalg.norm(final)),norm_reduced=float(np.linalg.norm(corrected)),
                         common_fwd_modal_amplitude=float(abs(decomp[j])),
                         integrated_instantaneous_acoustic_amplitude=float(np.exp(exp_damping)),
                         reconstructed_vector_relative_error=float(np.linalg.norm(corrected-final)/np.linalg.norm(final))))
    return dict(T=T,k_over_elastic_heat_scale=k,sweep=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('receipt')
    a=p.parse_args()
    z=json.loads(Path(a.receipt).read_text())
    print(json.dumps(dict(scope="Imposed chemistry clock with held temperature; not a cosmological transfer model.",
                          states=[run(r) for r in z['rows'][:2]]),indent=2,allow_nan=False))
