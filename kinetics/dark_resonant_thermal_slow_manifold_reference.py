"""Independent fixed-number common-temperature slow manifold for resonant heat."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from conservative_dark_resonant_feedback_reference import fit,thermo
from resonant_on_shell_energy_exchange_reference import thermal_heat,quadr,C
from three_species_decay_onsager_reference import M,mp,Es,P
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from scipy.special import expit

def conductance(T,eta_chi,eta_phi,order=(74,46)):
    z,w=quadr(order[0]);x,v=quadr(order[1])
    y=(z+1)*np.sqrt(68.)/2
    E=M+T*y*y
    dE=w*T*np.sqrt(68.)*y
    p=np.sqrt((E-M)*(E+M))
    u=(E*Es+p*P)/M
    l=(Es*Es+(p*mp/M)**2)/u
    q=np.log(u/l)
    e=l[:,None]*np.exp(q[:,None]*(1+x)/2)
    b=q[:,None]*v*e/2*expit((E[:,None]-e)/T-eta_chi)
    b0=b.sum(axis=1)
    mean=(b*e).sum(axis=1)/b0
    variance=np.sum(b*(e-mean[:,None])**2,axis=1)
    return C*np.sum(dE*np.exp(eta_chi+eta_phi-E/T)*variance)/T**2

def run(row,eos):
    TR=row["TR"];TS=row["TF"];TF=max(.01,.65*TR)
    ss=eos.state(TS)[2]
    e0=eos.state(TF)[2]
    loga=np.log(ss/e0)/3
    n,Td=fit(row,eos)
    rho=np.array([thermo(t,x,i==0)[0] for i,(t,x) in enumerate(zip(Td,n))])
    def common(t,num):
        return [thermo(t,x,i==0) for i,x in enumerate(num)]
    Teq=brentq(lambda t:sum(j[0] for j in common(t,n))-rho.sum(),min(Td),max(Td),xtol=1e-14)
    init=common(Teq,n)
    ent=sum(j[-1] for j in init)
    f=n*np.exp(-3*loga)
    T1=brentq(lambda t:sum(j[-1] for j in common(t,f))-ent*np.exp(-3*loga),.25*Teq,Teq,xtol=1e-14)
    final=common(T1,f)
    qout=[]
    for t,z0,ts in ((Teq,init,TS),(T1,final,TF)):
        K=conductance(t,z0[0][4],z0[1][4])
        eps=2e-6
        Knew=(thermal_heat(t*(1-eps),t,z0[1][4],z0[0][4])-
              thermal_heat(t*(1+eps),t,z0[1][4],z0[0][4]))/(2*eps*t)
        check=abs(Knew/K-1)
        assert check<2e-7,check
        H=np.sqrt((eos.state(ts)[0]+sum(y[0] for y in z0))/3)/MPL
        qout.append(dict(T_visible=ts,T_common=t,conductance=K,
                         relative_derivative_check=check,
                         thermal_mode_per_H=K*sum(1/y[2] for y in z0)/H))
    def temperature(l):
        den=n*np.exp(-3*l)
        t=brentq(lambda v:sum(j[-1] for j in common(v,den))-ent*np.exp(-3*l),
                  .55*T1,1.25*Teq,xtol=1e-14)
        return t,common(t,den)
    def trace(i):
        def f(l):
            v=temperature(l)[1][i]
            return np.exp(4*l)*(v[0]-3*v[1])/rho.sum()
        return quad(f,0.,loga,epsabs=2e-9,epsrel=2e-9)[0]
    trace_v=np.array([trace(0),trace(1)])
    late=np.exp(4*loga)*np.array([j[0] for j in final])
    # The initial fast thermal sharing plus continuous redistribution needs
    # a nonzero exact energy ledger even as the temperature gap tends to zero.
    heat=(late-rho)/rho.sum()-trace_v
    assert abs(sum(heat))<1e-7
    return dict(TR=TR,visible_start=TS,visible_end=TF,
                initial_fitted_temperatures_chi_phi=Td.tolist(),
                common_after_initial_energy_sharing=Teq,terminal_common=T1,
                cumulative_opposite_heat_chi_phi=heat.tolist(),
                relative_total_energy_residual=abs(sum(heat)),
                conductance_checks=qout)

if __name__=="__main__":
    q=argparse.ArgumentParser();q.add_argument("eos");q.add_argument("spectra")
    a=q.parse_args()
    eos=StandardModel(a.eos)
    rows=[run(r,eos) for r in json.loads(Path(a.spectra).read_text())["rows"]
          if r['TR'] in (.02,.05,.1) and abs(r['TF']/r['TR']-.9)<1e-6]
    print(json.dumps(dict(rows=rows,scope="Strong-resonant-thermal-exchange moment control, not a kinetic spectrum."),indent=2,allow_nan=False))
