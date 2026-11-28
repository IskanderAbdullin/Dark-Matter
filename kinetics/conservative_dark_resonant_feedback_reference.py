"""Closed two-temperature number/energy closure with resonant heat feedback."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from sm_dark_entropy_exchange_reference import C,MPL,StandardModel
from three_species_decay_onsager_reference import mp,m
from fd_number_heat_state_reference import state as fd
from mb_number_heat_state_reference import state as mb
from resonant_on_shell_energy_exchange_reference import thermal_heat

def thermo(t,n,fermion):
    D=fd(t,n,m) if fermion else mb(t,n,mp)
    # The generic helper also carries an independent massless radiation
    # component; this is a two-species test, so remove that component exactly.
    rho=D['rho']-C*t**4
    p=D['pressure']-C*t**4/3
    c=D['cv_fixed_number']-4*C*t**3
    E=D['chemical_energy']
    eta=D['eta']+(0. if fermion else np.log(4.)) # g=1 for real scalar
    S=(rho+p)/t-eta*n
    assert rho>0 and c>0 and S>0
    return rho,p,c,E,eta,S

def fit(row,eos):
    t=row['TF']
    vis=eos.state(t)[2]
    Y=np.array(row['stable_yields'])
    n=vis*np.array([Y.sum(),Y[1]])
    e=np.array([Y@np.array(row['redshifted_means'][:2])/Y.sum(),row['redshifted_means'][2]])
    T=[brentq(lambda v:thermo(v,d,i==0)[0]/d-z,t*.2,t*8,xtol=1e-14)
       for i,(d,z) in enumerate(zip(n,e))]
    return n,np.array(T)

def run(row,eos,scale=1.,rtol=1e-8,order=(62,38),solver="Radau"):
    TR=row['TR']
    T0=row['TF']
    Tend=max(.01,.65*TR)
    ss=eos.state(T0)[2]
    ee=eos.state(Tend)[2]
    xend=np.log(ss/ee)/3
    node=np.geomspace(Tend,T0,1500)
    v=[eos.state(T) for T in node]
    X=np.log(ss/np.array([j[2] for j in v]))/3
    table=PchipInterpolator(X[::-1],np.array([j[0] for j in v])[::-1])
    N,T=fit(row,eos)
    a=[thermo(t,n,i==0) for i,(t,n) in enumerate(zip(T,N))]
    rho0=np.array([j[0] for j in a])
    ent0=sum(j[-1] for j in a)
    normal=rho0.sum()
    def rhs(l,y):
        t=T*np.exp(y[:2]); n=N*np.exp(-3*l)
        (r1,p1,c1,e1,h1,s1),(r2,p2,c2,e2,h2,s2)=[thermo(k,d,i==0) for i,(k,d) in enumerate(zip(t,n))]
        H=np.sqrt((table(l)+r1+r2)/3)/MPL
        q=thermal_heat(t[1],t[0],h2,h1,order=order)*scale/H if scale else 0.
        # Projected heat transfer has no number source.
        v=np.array([(-3*(r1+p1-n[0]*e1)-q)/(t[0]*c1),
                    (-3*(r2+p2-n[1]*e2)+q)/(t[1]*c2)])
        return np.r_[v,np.exp(4*l)*np.array([r1-3*p1,r2-3*p2,q])/normal,
                    np.exp(3*l)*q*(t[0]-t[1])/(t[0]*t[1]*ent0)]
    y0=np.zeros(6)
    s=solve_ivp(rhs,(0.,xend),y0,method=solver,rtol=rtol,atol=rtol*.01,max_step=.04,first_step=min(.001,.02/max(1.,np.max(abs(rhs(0.,y0)[:2])))))
    assert s.success,s.message
    t=T*np.exp(s.y[:2,-1]);n=N*np.exp(-3*xend)
    final=[thermo(k,d,i==0) for i,(k,d) in enumerate(zip(t,n))]
    x=s.y[:,-1]
    actual=np.exp(xend*4)*np.array([v[0] for v in final])
    expect=rho0+normal*(x[2:4]+x[4]*np.array([-1.,1.]))
    entropy=np.exp(xend*3)*sum(v[-1] for v in final)
    err=np.max(abs(actual-expect))/normal
    es=abs(entropy/ent0-1-x[-1])
    assert err<2e-6 and es<2e-6 and x[-1]>=-1e-12,(err,es)
    return dict(TR=TR,visible_start=T0,visible_end=Tend,initial_temperature_chi_phi=T.tolist(),
                final_temperature_chi_phi=t.tolist(),final_relative_chi_phi_temperature_gap=float(abs(t[1]/t[0]-1)),
                final_fractional_comoving_scalar_heat_vs_initial_dark_energy=x[4],
                entropy_production_relative_initial=x[-1],restriction_scale=scale,
                energy_residual=float(err),entropy_residual=float(es),nfev=s.nfev)

def main(ef,datafile):
    eos=StandardModel(ef)
    data=json.loads(Path(datafile).read_text())
    x=[]
    for row in data['rows']:
        if row['TR'] not in (.02,.05,.1) or abs(row['TF']/row['TR']-.9)>1e-6:continue
        a=run(row,eos)
        b=run(row,eos,scale=0.)
        c=run(row,eos,rtol=3e-9,order=(86,48),solver='BDF')
        dev=np.max(abs(np.array(c['final_temperature_chi_phi'])/a['final_temperature_chi_phi']-1.))
        assert dev<2e-5,dev
        x.append(dict(exchange=a,zero=b,independent_stiff_solver_and_quadrature=c,relative_terminal_refine=dev))
    return dict(rows=x,scope="Two-temperature FD/MB moment closure after a single prescribed injection snapshot; number fixed, resonant energy only. Not actual kinetic thermalization or a full abundance.")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos");p.add_argument("spectra")
    a=p.parse_args();print(json.dumps(main(a.eos,a.spectra),indent=2,allow_nan=False))
