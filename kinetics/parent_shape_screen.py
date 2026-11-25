"""Momentum-resolved parent *bath-only* closure diagnostic.

Uses saved cosmological theta(x), H(x) histories as an external background.
Compares momentum-dependent effective SM relaxation with the number-only
approximation from YukawaRateBath. Dark 1<->2 reactions are OFF in both
calculations: this is not a new scalar-production solution, nor a validation
of mother kinetic equilibrium during stimulated depletion.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre,expit
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator


def run(row,nq=180,steps_per_log=1000):
    hist=row['history'];x0=row['x0'];g=row['g_bath'];rad0=np.pi**2*g/30
    xlist=np.r_[x0,[h['x'] for h in hist]];zlist=np.log(xlist)
    theta=PchipInterpolator(zlist,np.log(np.r_[1.,[h['theta'] for h in hist]]))
    rho=PchipInterpolator(zlist,np.log(np.r_[row['initial_total_rho_Tref4'],[h['total_rho_Tref4'] for h in hist]]))
    p,w=roots_legendre(nq);p=40*(p+1);w=40*w*4*p*p/(2*np.pi**2)
    k1=row['kappa_sm1'];gamma=row['gamma'] if row['scatter'] else 0.
    def coefficients(z):
        x=np.exp(z);T=np.exp(theta(z));Hf=np.sqrt(rad0/np.exp(rho(z)))
        E=np.hypot(p,x);eq=expit(-E/T);ne=w@eq
        dec=k1*x*x*Hf*x/E
        sc=8*np.pi*gamma*k1*x*T*Hf/(1+(x/T)**2)
        rate=dec+sc
        avg=(w@(rate*eq))/max(ne,1e-300)
        return E,eq,rate,avg,T
    E,eq,_,_,T=coefficients(zlist[0])
    f=eq.copy() if row.get('initial_parent','thermal')=='thermal' else np.zeros_like(eq)
    N=w@f
    out=[]
    for a,b in zip(zlist[:-1],zlist[1:]):
        n=max(1,int(np.ceil((b-a)*steps_per_log)));dt=(b-a)/n
        for j in range(n):
            E,eq,rate,avg,T=coefficients(a+(j+.5)*dt)
            f=eq+(f-eq)*np.exp(-rate*dt)
            ne=w@eq
            N=ne+(N-ne)*np.exp(-avg*dt)
        x=np.exp(b)
        if x>1:break
        E,eq,rate,avg,T=coefficients(b);num=w@f
        xi=brentq(lambda u:w@expit(u-E/T)-num,-750,80)
        ff=expit(xi-E/T)
        # Direct spontaneous scalar decay seed, no daughter blocking or stimulation.
        r=row['fermion_ratio'];m=x*row['scalar_ratio']
        estar=x*(1-r*r+row['scalar_ratio']**2)/2
        pstar=np.sqrt((estar-m)*(estar+m))
        elo=((x*x/(E+p))*estar+p*m*m/(estar+pstar))/x
        ehi=(E*estar+p*pstar)/x
        def seed(f0,cut):
            ecut=np.hypot(m,cut*T)
            frac=np.clip((ecut-elo)/(ehi-elo),0,1)
            return float(w@(f0*x/E*frac))
        seeds={}
        for q in [1e-5,1e-4,1e-3,.1]:
            ref=seed(ff,q);actual=seed(f,q)
            seeds[str(q)]=actual/ref if ref>1e-80 else None
        out.append(dict(x=x,number_resolved=float(num),number_closure=float(N),
            abundance_ratio=float(num/N),normalized_shape_TV=float(.5*(w@abs(f-ff))/num),
            energy_ratio_same_number=float((w@(f*E))/(w@(ff*E))),
            high_p_fraction_actual=float((w[p>5*T]@f[p>5*T])/num),
            high_p_fraction_FD=float((w[p>5*T]@ff[p>5*T])/num),
            spontaneous_low_p_seed_ratio_same_number=seeds))
    return dict(x0=x0,kappa_sm1=k1,initial_parent=row.get('initial_parent','thermal'),
        states=out,scope='Bath only; prescribed background; excludes dark depletion, no scalar abundance prediction')


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',nargs='+',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();results=[]
    for file in a.input:
        for row in json.loads(Path(file).read_text()):
            s=run(row);t=run(row,240,2000)
            drift=max(abs(a['abundance_ratio']/b['abundance_ratio']-1) for a,b in zip(s['states'],t['states']))
            assert drift<5e-5,drift
            s['abundance_numerical_drift']=drift;results.append(s)
            print(json.dumps(dict(x0=s['x0'],kappa_sm1=s['kappa_sm1'],
                max_abundance_closure_error=max(abs(h['abundance_ratio']-1) for h in s['states']),
                max_shape_TV=max(h['normalized_shape_TV'] for h in s['states']),
                max_energy_error=max(abs(h['energy_ratio_same_number']-1) for h in s['states']),
                abundance_numerical_drift=drift)),flush=True)
    Path(a.output).write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
