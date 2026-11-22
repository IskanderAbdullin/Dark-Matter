"""Independent relativistic finite-mediator Pauli 2->2 WEAK moments.

Generalize a same-model χℓ↔χℓ on-shell CM 5D angular event to
ARBITRARY isotropic nonthermal χ(E), ℓ(ω) occupation callbacks:
full F=fχ fℓ (1-fχ')(1-fℓ'), R=fχ' fℓ' (1-fχ)(1-fℓ).
Physical momentum/energy invariant per event and entropy ≥0. Spin
averaged finite-vector amplitude, moving frame and exact lepton mass.

Only incoming-integrated WEAK heat/entropy moments with constant
prescribed f callbacks: no strong Cχ(p), Cℓ(k), no finite bath time,
actual χ1/χ2 coupled abundance, SM weak chemistry or FRW.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss


def gauss(a,b,n):
    z,w=leggauss(n)
    return a+(b-a)*(z+1)/2,w*(b-a)/2


def event(M,m,p,k,mu,G,mV,c,phi,fchi,flepton):
    """One incoming χ+lepton pair, scattering CM polar c and azimuth φ.

    fchi,flepton accept lab-frame energy (scalar or NumPy array) and
    return occupations strictly between zero and one (no Pauli clipping).
    """
    if not(M>0 and m>=0 and p>0 and k>0 and -1<mu<1 and
           G>0 and mV>0 and -1<=c<=1 and
           all(math.isfinite(v) for v in (M,m,p,k,mu,G,mV,c))):
        raise ValueError("physical on-shell 2to2 scattering input required")
    E=math.hypot(p,M)
    w=math.hypot(k,m)
    Et=E+w
    P=math.sqrt(p*p+k*k+2*p*k*mu)
    s=Et*Et-P*P
    if P==0 or s<=(M+m)**2:
        raise ValueError("closed/degenerate CM phase space")
    root=math.sqrt(s)
    ps2=(s-(M+m)**2)*(s-(M-m)**2)/(4*s)
    if ps2<=0:raise ValueError("zero outgoing CM momentum")
    ps=math.sqrt(ps2)
    beta=P/Et
    gamma=Et/root
    parallel=p*(p+k*mu)/P
    ca=max(-1.,min(1.,gamma*(parallel-beta*E)/ps))
    Estar=(s+M*M-m*m)/(2*root)
    phi=np.asarray(phi,float)
    cp=ca*c+math.sqrt(max(0.,1-ca*ca))*math.sqrt(max(0.,1-c*c))*np.cos(phi)
    Ep=gamma*(Estar+beta*ps*cp)
    wp=Et-Ep
    t=-2*ps2*(1-c)
    D=s-M*M-m*m
    U=-D-t
    amplitude=2*G*G*(D*D+U*U+2*(M*M+m*m)*t)*(mV*mV/(mV*mV-t))**2
    f=float(fchi(E));g=float(flepton(w))
    fp=np.asarray(fchi(Ep),float)
    gp=np.asarray(flepton(wp),float)
    if (not 0<f<1 or not 0<g<1 or
        fp.shape!=np.shape(Ep) or gp.shape!=np.shape(wp) or
        np.any(fp<=0) or np.any(fp>=1) or
        np.any(gp<=0) or np.any(gp>=1) or
        not np.all(np.isfinite(fp)) or not np.all(np.isfinite(gp)) or
        not math.isfinite(amplitude) or amplitude<0):
        raise ValueError("finite physical FD interior and nonnegative amplitude needed")
    forward=f*g*(1-fp)*(1-gp)
    reverse=fp*gp*(1-f)*(1-g)
    flux=forward-reverse
    entropy=flux*np.log(forward/reverse)
    dot=E*w-p*k*mu
    vm=math.sqrt(max(0.,dot*dot-M*M*m*m))/(E*w)
    return {"E":E,"omega":w,"Ep":Ep,"omega_p":wp,"dEchi":Ep-E,
            "dElepton":wp-w,"flux":flux,"forward":forward,
            "reverse":reverse,"entropy":entropy,"amplitude":amplitude,
            "s":s,"vm":vm,
            "exact_local_energy_balance":float(np.max(abs((Ep-E)+(wp-w))))}


def weak_moments(fchi,flepton,*,M=.1,m=0.,Te=1.,G=.001,mV=.1,
                 nr=4,nmu=4,nphi=12,nang=12,cut_kinetic_T=12.):
    """Full incoming 5D isotropic χℓ scattering WEAK heat+entropy.

    This uses the historical spin/degeneracy convention gχ=gℓ=4 as a
    numerical reference, not species-complete SM χ/ν flavor counting.
    """
    if not (M>0 and m>=0 and Te>0 and G>0 and mV>0 and
            nr>=4 and nmu>=4 and nphi>=8 and nang>=4 and
            cut_kinetic_T>0):
        raise ValueError
    pmax=math.sqrt((M+cut_kinetic_T*Te)**2-M*M)
    kmax=math.sqrt((m+cut_kinetic_T*Te)**2-m*m)
    ps,wp=gauss(0.,pmax,nr)
    ks,wk=gauss(0.,kmax,nr)
    mus,wmu=gauss(-1.,1.,nmu)
    cs,wc=leggauss(nang)
    phi=np.arange(nphi,dtype=float)*2*math.pi/nphi
    heat=entropy=0.
    max_pair_energy_defect=0.
    for p,ap in zip(ps,wp):
        for k,ak in zip(ks,wk):
            for u,au in zip(mus,wmu):
                integ_q=integ_s=0.
                for c,ac in zip(cs,wc):
                    ev=event(M,m,float(p),float(k),float(u),G,mV,
                             float(c),phi,fchi,flepton)
                    integ_q+=ac*ev["amplitude"]*float(np.mean(ev["dEchi"]*ev["flux"]))
                    integ_s+=ac*ev["amplitude"]*float(np.mean(ev["entropy"]))
                    max_pair_energy_defect=max(max_pair_energy_defect,
                                               ev["exact_local_energy_balance"])
                weight=ap*p*p*ak*k*k*au*ev["vm"]/ev["s"]
                heat+=weight*integ_q
                entropy+=weight*integ_s
    pref=4*4/(512*math.pi**5)
    q=heat*pref
    return {"Qchi":q,"Qlepton":-q,"total_energy_exchange_residual":0.,
            "entropy_production":entropy*pref,
            "max_local_pair_energy_defect":max_pair_energy_defect,
            "nonthermal_isotropic_weak_5D_only":True,
            "not_strong_collision_or_physical_gamma_el":True,
            "not_86eV_DM_prediction":True}


if __name__=="__main__":
    import json
    from scipy.special import expit
    cases={}
    for mass in (0.,.5):
        fchi=lambda E:expit(-np.asarray(E)/.9)
        flep=lambda E:expit(-np.asarray(E))
        cases[str(mass)]={"cold":weak_moments(fchi,flep,m=mass),
           "equilibrium":weak_moments(lambda E:expit(-np.asarray(E)),flep,m=mass)}
    fchi=lambda E:expit(-np.asarray(E)/.9+.15*np.exp(-np.asarray(E)/.9))
    flep=lambda E:expit(-np.asarray(E)+.05*np.exp(-np.asarray(E)))
    cases["nonthermal"]=weak_moments(fchi,flep)
    print(json.dumps(cases,indent=2))
