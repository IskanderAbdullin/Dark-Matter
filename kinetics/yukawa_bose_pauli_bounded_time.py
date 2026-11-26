"""Energy-matched finite-volume χ1 <-> χ2 + real φ Pauli/Bose quantum pilot.

Every forward and inverse moving on-shell event carries one common positive
phase-space weight and is deposited barycentrically in three distinct physical
energy grids. Fermi log-odds and Bose inverse log-fugacity are interpolated
with the SAME weights as deposition. This enforces exact discrete stoichiometric
charges, first energy moment, and an eventwise H-theorem for interior f,b.

Physical densities use integrated on-shell radial phase-space shell volumes.
An entropy/positivity-gated Euler step-doubling scheme recomputes all three
species' sources; it is FIRST ORDER, not a certified global integrator.
No finite SM/lepton thermal bath, FRW, actual 86-eV benchmark or relic claim.
"""
from __future__ import annotations
import math
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit, logit, xlogy
from kinetics.moving_yukawa_bose_pauli_weak import vacuum_width, daughter_energies


def density_cells(nodes, mass, spin):
    e=np.asarray(nodes, float)
    if not (len(e)>=4 and mass>=0 and spin>0 and abs(e[0]-mass)<1e-12
            and np.all(np.diff(e)>0) and
            np.allclose(np.diff(e),e[1]-e[0],rtol=2e-12)):
        raise ValueError('uniform physical E nodes start at species rest energy')
    d=e[1]-e[0]
    lo=np.maximum(e[0],e-d/2)
    hi=np.minimum(e[-1],e+d/2)
    return spin/(6*math.pi**2)*(np.maximum(0,hi*hi-mass*mass)**1.5-
                                   np.maximum(0,lo*lo-mass*mass)**1.5)


def barycentric(nodes, energy):
    e=np.asarray(energy,float)
    if not np.all(np.isfinite(e)) or np.any(e<nodes[0]-2e-12) or np.any(e>nodes[-1]+2e-12):
        raise ValueError('on-shell event outside finite energy grid')
    h=nodes[1]-nodes[0]
    t=np.clip((e-nodes[0])/h,0,len(nodes)-1-1e-13)
    i=np.floor(t).astype(int)
    a=t-i
    return i,a


def interpolation(nodes, values, energy):
    i,a=barycentric(nodes,energy)
    return (1-a)*values[i]+a*values[i+1]


def deposit(dst, nodes, energy, amount):
    e=np.atleast_1d(np.asarray(energy,float))
    r=np.broadcast_to(amount,e.shape)
    i,a=barycentric(nodes,e)
    np.add.at(dst,i,(1-a)*r)
    np.add.at(dst,i+1,a*r)


def fermion_entropy(f,w):
    if np.any(f<0) or np.any(f>1):raise ValueError('fermion bounds')
    return -float(np.dot(w,xlogy(f,f)+xlogy(1-f,1-f)))


def boson_entropy(b,w):
    if np.any(b<0):raise ValueError('boson bounds')
    return float(np.dot(w,xlogy(1+b,1+b)-xlogy(b,b)))


@dataclass
class Config:
    m1: float=.1
    m2: float=.02
    mphi: float=.001
    y: float=.001
    pmax: float=12.
    bins: int=64
    nr: int=16
    nangle: int=16
    g1: int=4
    g2: int=4
    gphi: int=1

    def __post_init__(self):
        if not(vacuum_width(self.m1,self.m2,self.mphi,self.y)>0 and
               self.pmax>0 and self.bins>=8 and self.nr>=2 and self.nangle>=2):
            raise ValueError('open Yukawa channel and valid numerical grid required')
        upper=math.hypot(self.m1,self.pmax)+1e-10
        self.e1=np.linspace(self.m1,upper,self.bins+1)
        self.e2=np.linspace(self.m2,upper,self.bins+1)
        self.ephi=np.linspace(self.mphi,upper,self.bins+1)
        self.w1=density_cells(self.e1,self.m1,self.g1)
        self.w2=density_cells(self.e2,self.m2,self.g2)
        self.wphi=density_cells(self.ephi,self.mphi,self.gphi)
        self.width=vacuum_width(self.m1,self.m2,self.mphi,self.y)


def thermal_initial(cfg, T1=1.1,T2=1.,Tphi=1.,mu1=0.,mu2=0.,muphi=0.):
    if not(T1>0 and T2>0 and Tphi>0 and muphi<cfg.mphi):
        raise ValueError('thermal Bose chemical potential must below mphi')
    return (expit((mu1-cfg.e1)/T1),
            expit((mu2-cfg.e2)/T2),
            1/np.expm1((cfg.ephi-muphi)/Tphi))


def ledger(cfg,f1,f2,b):
    n1=float(np.dot(cfg.w1,f1));n2=float(np.dot(cfg.w2,f2));nphi=float(np.dot(cfg.wphi,b))
    return {'N1_plus_N2':n1+n2,'N1_plus_Nphi':n1+nphi,
            'energy':float(np.dot(cfg.w1*cfg.e1,f1)+np.dot(cfg.w2*cfg.e2,f2)+
                           np.dot(cfg.wphi*cfg.ephi,b)),
            'entropy':fermion_entropy(f1,cfg.w1)+fermion_entropy(f2,cfg.w2)+boson_entropy(b,cfg.wphi),
            'Nchi1':n1,'Nchi2':n2,'Nphi':nphi}


def rhs(cfg,f1,f2,b):
    if (np.any(f1<=0) or np.any(f1>=1) or np.any(f2<=0) or np.any(f2>=1)
            or np.any(b<=0) or not np.all(np.isfinite(f1)) or not np.all(np.isfinite(f2))
            or not np.all(np.isfinite(b))):
        raise ValueError('strict interior Pauli/Bose state required')
    h1,h2=logit(f1),logit(f2)
    eta=np.log1p(1/b)
    dp1=np.zeros_like(f1);dp2=np.zeros_like(f2);dpb=np.zeros_like(b)
    ps,ws=leggauss(cfg.nr)
    ps=(ps+1)*cfg.pmax/2;ws=ws*cfg.pmax/2
    cc,wc=leggauss(cfg.nangle)
    ent=heat1=heat2=heatb=rate=0.
    maxdefect=0.
    for p,wp in zip(ps,ws):
        E1=math.hypot(cfg.m1,float(p))
        E2,Eb=daughter_energies(float(p),cc,cfg.m1,cfg.m2,cfg.mphi)
        fs1=expit(interpolation(cfg.e1,h1,E1))
        fs2=expit(interpolation(cfg.e2,h2,E2))
        ephi=interpolation(cfg.ephi,eta,Eb)
        bs=1/np.expm1(ephi)
        forward=fs1*(1-fs2)*(1+bs)
        reverse=(1-fs1)*fs2*bs
        j=forward-reverse
        micro=cfg.g1/(2*math.pi**2)*wp*p*p*cfg.width*cfg.m1/E1/2
        amounts=micro*wc*j
        tot=float(sum(amounts))
        deposit(dp1,cfg.e1,E1,-tot)
        deposit(dp2,cfg.e2,E2,amounts)
        deposit(dpb,cfg.ephi,Eb,amounts)
        rate+=tot
        heat1-=E1*tot
        heat2+=float(np.dot(E2,amounts))
        heatb+=float(np.dot(Eb,amounts))
        ent+=float(micro*np.dot(wc,j*np.log(forward/reverse)))
        maxdefect=max(maxdefect,float(np.max(abs(E1-E2-Eb))))
    sbin= -float(np.dot(dp1,h1)+np.dot(dp2,h2))+float(np.dot(dpb,eta))
    return (dp1/cfg.w1,dp2/cfg.w2,dpb/cfg.wphi,
            {'Rphi':rate,'Rchi1':-rate,'Rchi2':rate,
             'Qchi1':heat1,'Qchi2':heat2,'Qphi':heatb,
             'heat_residual':heat1+heat2+heatb,
             'number_residual_chi':float(sum(dp1)+sum(dp2)),
             'number_residual_phi':float(sum(dp1)+sum(dpb)),
             'entropy_event':ent,'entropy_nodes':sbin,
             'entropy_identity_abs':abs(ent-sbin),
             'max_event_energy_defect':maxdefect})


def euler_guard(cfg,state,dt):
    d1,d2,db,info=rhs(cfg,*state)
    new=(state[0]+dt*d1,state[1]+dt*d2,state[2]+dt*db)
    if any(np.any(x<=0) for x in new) or any(np.any(x>=1) for x in new[:2]):
        return None
    if ledger(cfg,*new)['entropy']+1e-12<ledger(cfg,*state)['entropy']:
        return None
    return new,info


def evolve(cfg,initial, *, end_time=5.e7,initial_dt=2.e7,atol=1.e-6,max_attempts=2500):
    if not(end_time>0 and initial_dt>0 and atol>0):raise ValueError
    state=tuple(np.array(x,float,copy=True) for x in initial)
    first=ledger(cfg,*state)
    t=0.;dt=initial_dt;accepted=0;rejected=0;max_err=0.;min_S=math.inf
    for _ in range(max_attempts):
        if t>=end_time*(1-1e-14):break
        dt=min(dt,end_time-t)
        full=euler_guard(cfg,state,dt)
        half=euler_guard(cfg,state,dt/2)
        double=None if half is None else euler_guard(cfg,half[0],dt/2)
        if full is not None and double is not None:
            err=max(float(np.max(abs(full[0][0]-double[0][0]))),
                    float(np.max(abs(full[0][1]-double[0][1]))),
                    float(np.max(abs(np.log1p(full[0][2])-np.log1p(double[0][2])))))
            if err<=atol:
                old_S=ledger(cfg,*state)['entropy']
                state=tuple(x.copy() for x in double[0])
                min_S=min(min_S,ledger(cfg,*state)['entropy']-old_S)
                t+=dt;accepted+=1;max_err=max(max_err,err)
                dt*=min(2.,max(1.1,.9*math.sqrt(atol/max(err,1e-30))))
                continue
        rejected+=1;dt*=.5
        if dt<1e-16:raise ArithmeticError('bounded Pauli/Bose step fails')
    else:raise RuntimeError('too many adaptive steps')
    final=ledger(cfg,*state)
    return {'accepted':accepted,'rejected':rejected,'time':t,
            'first':first,'last':final,
            'd_N1_N2':final['N1_plus_N2']-first['N1_plus_N2'],
            'd_N1_Nphi':final['N1_plus_Nphi']-first['N1_plus_Nphi'],
            'd_energy':final['energy']-first['energy'],
            'entropy_gain':final['entropy']-first['entropy'],
            'min_step_entropy_gain':min_S,
            'max_local_full_vs_2half_spectral_error':max_err,
            'min_fermion':min(float(np.min(state[0])),float(np.min(state[1]))),
            'max_fermion':max(float(np.max(state[0])),float(np.max(state[1]))),
            'min_bose':float(np.min(state[2])),
            'max_bose':float(np.max(state[2])),
            'chi2_number_created':final['Nchi2']-first['Nchi2'],
            'phi_number_created':final['Nphi']-first['Nphi'],
            'final_f1_low':float(state[0][0]),'final_phi_low':float(state[2][0]),
            'first_order_step_doubling_only':True,
            'not_a_cosmological_86_eV_relic':True}


if __name__=='__main__':
    import json
    cfg=Config(bins=32,nr=12,nangle=12)
    state=thermal_initial(cfg)
    _,_,_,instant=rhs(cfg,*state)
    print(json.dumps({'instantaneous':instant,
                      'evolution':evolve(cfg,state)},indent=2))
