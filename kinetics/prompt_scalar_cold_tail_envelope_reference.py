"""Analytic thermal injection upper bound for a prompt-decay soft tail."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp,gammaln
from decay_scalar_ir_relaxation_reference import limits
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from three_species_decay_onsager_reference import M,P
from timelike_portal_onsager_coverage_reference import V,width,gprod

def log_cap(T,p):
    # At p<P (rest-frame decay momentum), the minimum parent energy
    # decreases monotonically as p increases; hence use p_max in the cut.
    assert 0<p<P
    _,minE,_=limits(p)
    x=(minE+M)/T # partner in the same pair costs at least its rest mass
    b=V*width
    spec=(np.hypot(V*V,b)+V*V)/(2*b*b) # max_s s / BW(s)
    max_sigma_v=3*gprod*gprod*spec/(8*np.pi)
    j=np.arange(6)
    gamma6=gammaln(6)-x+logsumexp(j*np.log(x)-gammaln(j+1))
    # Three l+ l- species, 2 spin each, d^3p/(2pi)^3≤E²dE/(2pi²).
    # Convolution in the total energy is E_tot^5 /30.
    return np.log(12*max_sigma_v/(120*np.pi**4))+6*np.log(T)+gamma6,x

def main(eos_file,spectra):
    data=json.loads(Path(spectra).read_text())
    eos=StandardModel(eos_file)
    out=[]
    for r in data['rows']:
        tr,tf=r['TR'],r['TF']
        sR=eos.state(tr)[2]
        rho,_,sF,_=eos.state(tf)
        H=np.sqrt(rho/3)/MPL
        lna=np.log(sR/sF)/3
        for a in (.005,.01,.025):
            p=a*tf*np.exp(lna)
            g,x=log_cap(tr,p)
            bound=np.log(2*lna/(H*sF))+g
            y0=r["stable_yields"][1]
            ratio=(bound-np.log(y0))/np.log(10.)
            out.append(dict(TR=tr,TF=tf,p_over_TF=a,
                              minimum_initial_pair_energy_over_TR=x,
                              log10_soft_yield_upper=bound/np.log(10.),
                              log10_upper_divided_by_numerical_total_scalar_yield=ratio))
    return dict(rows=out,scope="Analytic vacuum direct-lepton injection envelope only; prompt two-body decay, adiabatic fixed SM cooling, no later scalar or fermion interactions. Conditional on the specified widths and masses.")
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument('eos');p.add_argument('spectra')
    a=p.parse_args();print(json.dumps(main(a.eos,a.spectra),indent=2,allow_nan=False))
