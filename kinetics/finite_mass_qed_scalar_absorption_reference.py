"""Crossed P+scalar <-> daughter+photon on a fixed cold scalar mode."""
import gzip
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, roots_legendre

from early_pole_decay_mode_clock_reference import Q, intervals
from early_pole_hard_yukawa_photoproduction_reference import dot, amplitudes, invariant_form


@lru_cache(None)
def gl(n):
    return roots_legendre(n)


def coefficients(M,m,ms,J):
    D=(M*M+m*m-ms*ms)/2
    num=4*(2*D*D+2*D*J+2*D*M*m+J*J+2*J*M*m)/J
    off=(-8*D*J-4*D*m*m-8*J*J-8*J*M*m-4*M*m**3)/J**2
    return -4*M*M*(D+M*m),num,4/J,off


def azimuth(M,m,ms,J,P0,Ppar,Pperp,K,u):
    x=K*(P0-Ppar*u)
    y=K*Pperp*np.sqrt((1-u)*(1+u))
    d=np.sqrt((x-y)*(x+y))
    a,b,c,z=coefficients(M,m,ms,J)
    return a*x/d**3+b/d+c*x+z


def frame(E,v,M,m,ms):
    es=np.hypot(ms,Q)
    p=np.sqrt((E-M)*(E+M))
    W2=M*M+ms*ms+2*(E*es-p*Q*v)
    W=np.sqrt(W2)
    J=(W2-m*m)/2
    K=J/W
    et=E+es
    P=np.sqrt((p+Q*v)**2+Q*Q*(1-v*v))
    Ppar=((E+es)*p*(p+Q*v)/P-P*E)/W
    Pperp=p*Q*np.sqrt((1-v)*(1+v))/P
    P0=(W2+M*M-ms*ms)/(2*W)
    l=J/(et+P)
    h=J*(et+P)/W2
    return es,p,W,J,K,et,P,P0,Ppar,Pperp,l,h


def test_geometry(M,m,ms):
    ex,ey=np.eye(4)[1:3]
    ang,wa=gl(90)
    phi=np.pi*(ang+1)
    lower=(1-m*m/(M+ms)**2)*(np.hypot(Q,ms)-Q)/2
    observed=1e5
    max_fact=max_direct=max_gauge=max_proj=max_balance=0.
    for E in (M, M+.0002, .2,1.,5.,32.,80.):
        for v in (-.9,.21,.93):
            es,p,W,J,K,et,P,P0,Ppar,Pperp,l,h=frame(E,v,M,m,ms)
            if l<lower or not h>l:
                raise ArithmeticError('Global soft-exclusion bound failed.')
            observed=min(observed,l)
            max_proj=max(max_proj,abs(P0*P0-Ppar*Ppar-Pperp*Pperp-M*M))
            for u in (-.99,-.41,.23,.81,.999):
                spec=K*(P0-Ppar*u-Pperp*np.sqrt(1-u*u)*np.cos(phi))
                d=(M*M+m*m-ms*ms)/2-spec+J
                raw=invariant_form(-spec,-J,d,M,m)
                calc=azimuth(M,m,ms,J,P0,Ppar,Pperp,K,u)
                if min(raw.min(),calc)<=0.:
                    raise ArithmeticError('Nonpositive crossed matrix element')
                max_fact=max(max_fact,abs(calc/(wa@raw/2)-1.))
                photon=np.array([K,0.,0.,K])
                daughter=np.array([W-K,0.,0.,-K])
                r=Ppar*u+Pperp*np.sqrt(1-u*u)*.67
                trans=np.sqrt(max(0.,(Ppar*Ppar+Pperp*Pperp)-r*r))
                parent=np.array([P0,trans,0.,r])
                total=sum(float(amplitudes(parent,-photon,daughter,M,m,e)[0]) for e in (ex,ey))
                changed=sum(float(amplitudes(parent,-photon,daughter,M,m,e)[0])
                            for e in (ex+.7*photon,ey-.2*photon))
                closed=invariant_form(-dot(parent,photon),-J,dot(parent,daughter),M,m)
                max_direct=max(max_direct,abs(total/closed-1.))
                max_gauge=max(max_gauge,abs(changed/total-1.),abs(float(amplitudes(parent,-photon,daughter,M,m,photon)[0])/total))
                k=K/W*(et+P*u)
                for T in (.4,1.):
                    g=expit(-(et-k)/T)*expit(E/T)/np.expm1(k/T)
                    o=expit(-E/T)*expit((et-k)/T)/(-np.expm1(-k/T))
                    max_balance=max(max_balance,abs(g/(o*np.exp(-es/T))-1.))
    metrics=dict(exact_azimuth_relative_defect=max_fact,independent_amplitude_relative_defect=max_direct,
                 gauge_relative_defect=max_gauge,boost_mass_shell_abs_error=max_proj,
                 pointwise_equilibrium_relative_defect=max_balance,
                 rigorous_lab_photon_lower_bound=lower,minimum_sampled_lab_k=observed)
    if max(max_fact,max_direct,max_gauge,max_proj,max_balance)>1e-7:
        raise ArithmeticError(metrics)
    return metrics


def audit(d,T,cap,rad,polar,inner,equilibrium=False,overlap=False,energy_moments=False):
    M=d['x_end'];m=M*d['mass_ratio'];ms=d['checkpoints'][-1]['scalar_pole_over_T']
    es=np.hypot(ms,Q)
    p_edges,b_edges=[np.hypot(a,d[f'endpoint_{n}_momentum_edges'])
                     for a,n in ((M,'parent'),(m,'fermion'))]
    fp,fb=[np.r_[0.,d[f'endpoint_{k}_occupations'],0.] for k in ('parent','fermion')]
    fs=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],Q)-1]
    if equilibrium: fs=1/np.expm1(es/T)
    za,wa=gl(rad);zm,wm=gl(polar);zk,wk=gl(inner)
    total=np.zeros(2)
    heat=np.zeros(4) if energy_moments else None
    for v,wv in zip(zm,wm,strict=True):
        emin=M if not overlap else np.hypot(.45,M)
        if overlap:
            cap=np.hypot(5.,M)
            bmin,bmax=np.hypot([.45,5.],m)
            event_edges=np.r_[b_edges,bmin,bmax]
            common_radial_boundaries=np.r_[event_edges+.5-es,event_edges+4.-es]
        else:
            event_edges=b_edges
            common_radial_boundaries=[]
        base=np.unique(np.clip(np.r_[emin,cap,np.geomspace(emin,cap,50),common_radial_boundaries,
                                    [] if equilibrium else p_edges],emin,cap))
        if not equilibrium:
            # These are the lab extremal daughter energies at fixed parent momentum.
            # The low daughter endpoint can turn; scan adjacent radial segments.
            extra=[]
            for which in (10,11):
                def ed(E):
                    r=frame(E,v,M,m,ms)
                    return r[5]-r[which]
                xp=np.array([ed(x) for x in base])
                for i in range(len(base)-1):
                    aa,bb=sorted(xp[i:i+2])
                    for boundary in event_edges[(event_edges>aa)&(event_edges<bb)]:
                        root=brentq(lambda a:ed(a)-boundary,base[i],base[i+1],xtol=8e-14)
                        extra.append(root)
            base=np.unique(np.r_[base,extra])
        f=lambda E:expit(-E/T) if equilibrium else fp[np.searchsorted(p_edges,E)]
        e=(base[:-1,None]+base[1:,None])/2+(base[1:,None]-base[:-1,None])*za/2
        ew=((base[1:,None]-base[:-1,None])/2*wa).ravel()
        e=e.ravel()
        for start in range(0,len(e),64):
            E=e[start:start+64,None]
            es,p,W,J,K,et,P,P0,Ppar,Pperp,l,h=frame(E,v,M,m,ms)
            ll,hh=l,h
            if overlap:
                ll=np.maximum.reduce([l,np.full_like(l,.5),et-bmax])
                hh=np.maximum(ll,np.minimum.reduce([h,np.full_like(h,4.),et-bmin]))
            if equilibrium:
                edges=np.concatenate([ll,hh],axis=1)
            else:
                edges=np.sort(np.clip(np.concatenate([ll,hh,et-b_edges],axis=1),ll,hh),axis=1)
            lo,hi=edges[:,:-1],edges[:,1:]
            loglo=np.log(lo)
            dl=np.log1p((hi-lo)/lo)
            k=lo[:,:,None]*np.exp(dl[:,:,None]*(zk+1)/2)
            wk0=dl[:,:,None]/2*wk*k
            uv=np.clip((2*k-l[:,:,None]-h[:,:,None])/(h[:,:,None]-l[:,:,None]),-1.,1.)
            az=azimuth(M,m,ms,J[:,:,None],P0[:,:,None],Ppar[:,:,None],Pperp[:,:,None],
                        K[:,:,None],uv)
            if min(az.min(),p.min())<0.:
                raise ArithmeticError('Crossed scalar absorption has negative measure.')
            b=(expit(-(et[:,:,None]-k)/T) if equilibrium else
               fb[np.searchsorted(b_edges,et-(lo+hi)/2)][:,:,None])
            f0=f(E)[:,:,None]
            n=np.exp(-k/T)/(-np.expm1(-k/T))
            gain=(wk0*az*n*b*(1-f0)).sum(axis=(1,2))
            loss=(wk0*az*(1+n)*f0*(1-b)).sum(axis=(1,2))
            # K/W du = dk/P at fixed E,v. Incoming parent measure
            # p*dE*dv/(8pi²), COM pair measure K/W*du/(8pi).
            clock=(p/(128*np.pi**3*es*P))[:,0]
            norm=wv*(ew[start:start+64]*clock)
            total+=wv*(ew[start:start+64]*clock)@np.column_stack([gain,loss])
            if energy_moments:
                gd=wk0*az*n*b*(1-f0)
                ld=wk0*az*(1+n)*f0*(1-b)
                # For a scalar-producing event in this channel the photon
                # enters. Relativistic charged energy changes by k-Es.
                ks=((gd*k).sum(axis=(1,2)),(ld*k).sum(axis=(1,2)))
                dk=E[:,:,None]-(et[:,:,None]-k)
                ch=((gd*dk).sum(axis=(1,2)),(ld*dk).sum(axis=(1,2)))
                heat+=norm@np.column_stack([*ks,*ch])
    G,L=total
    N=G+fs*(G-L)
    out=dict(gain=float(G),loss=float(L),net=float(N),net_over_positive_turnover=float(N/((1+fs)*G+fs*L)))
    if energy_moments:
        kn=heat[0]+fs*(heat[0]-heat[1])
        ch=heat[2]+fs*(heat[2]-heat[3])
        o=es*N
        defect=o-kn+ch
        hnormal=abs(o)+abs(kn)+abs(ch)+abs(heat[0])*(1+fs)+abs(heat[1])*fs
        if abs(defect)>5e-13*hnormal:
            raise ArithmeticError(('Photon-bath heat ledger did not close',defect,hnormal))
        out.update(photon_bath_gain=float(-kn),scalar_energy_gain=float(o),
                   charged_energy_gain=float(ch),signed_heat_ledger_scaled_defect=float(defect/hnormal),
                   photon_weighted_scalar_gain=float(heat[0]),photon_weighted_scalar_loss=float(heat[1]),
                   photon_heat_over_unsigned_turnover=float(-kn/(heat[0]*(1+fs)+heat[1]*fs)))
    return out


def born(d):
    M=d['x_end'];m=M*d['mass_ratio'];ms=d['checkpoints'][-1]['scalar_pole_over_T']
    es,lo,hi,*_=intervals(M,ms,m,Q)
    a=np.hypot(M,d['endpoint_parent_momentum_edges'])
    b=np.hypot(m,d['endpoint_fermion_momentum_edges'])+es
    cuts=np.unique(np.clip(np.r_[lo,hi,a,b],lo,hi))
    y=(cuts[1:]+cuts[:-1])/2
    p=np.r_[0.,d['endpoint_parent_occupations'],0.][np.searchsorted(a,y)]
    f=np.r_[0.,d['endpoint_fermion_occupations'],0.][np.searchsorted(b,y)]
    fs=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],Q)-1]
    return float(2*((M+m)**2-ms*ms)/(16*np.pi*Q*es)*(np.diff(cuts)@(p*(1-f)+fs*(p-f))))


def gain_tail_bound(d,cap):
    M=d['x_end'];m=M*d['mass_ratio'];ms=d['checkpoints'][-1]['scalar_pole_over_T']
    es=np.hypot(ms,Q)
    emax=np.hypot(m,d['endpoint_fermion_momentum_edges'][-1])
    pmax=np.hypot(M,d['endpoint_parent_momentum_edges'][-1])
    if cap<=max(pmax,emax-es):
        raise ValueError('Tail bound requires neither the parent shell nor Bose enhancement at threshold')
    W0=M+ms
    Jmin=(W0*W0-m*m)/2
    D=(M*M+m*m-ms*ms)/2
    Amin=(1-m*m/(W0*W0))*M*M/2
    A0=4+4/Amin*(2*D*(D+M*m)/Jmin+2*D+2*M*m+D+ms*ms-m*m)
    A1=4/Amin*(es+Q)
    thermal=np.exp(-(cap+es-emax))/(-np.expm1(-(cap+es-emax)))
    value=thermal*(A0*(cap+1)+A1*(cap*cap+2*cap+2))/(64*np.pi**3*es)
    return dict(parent_energy_cap=cap,absolute_gain_tail_upper_bound=float(value),
                loss_tail_exactly_zero=True)


def run():
    fixtures=[]
    root=Path(__file__).parent/'results'
    with gzip.open(root/'early-pole-frozen-yukawa-real-photon-crossings-2026-11-27.json.gz') as f:
        archival=json.load(f)['grid']
    for g in ('fine','third'):
        path=Path(__file__).parent/'results'/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz'
        with gzip.open(path) as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];ms=d['checkpoints'][-1]['scalar_pole_over_T']
        controls=test_geometry(M,m,ms)
        rows=[audit(d,1.,56.,*mesh) for mesh in ((3,4,12),(6,6,20),(8,8,32))]
        emax=max(abs(v/rows[-2][k]-1.) for k,v in rows[-1].items())
        tail=audit(d,1.,42.,6,6,20)
        bound=gain_tail_bound(d,56.)
        if bound['absolute_gain_tail_upper_bound']>1e-9:
            raise ArithmeticError(bound)
        therm=[audit(d,T,38.,6,6,32,True) for T in (.4,1.)]
        if max(abs(x['net_over_positive_turnover']) for x in therm)>2e-12:
            raise ArithmeticError('Fixed-mode crossed scalar absorption violated equilibrium')
        if emax>8e-4: raise ArithmeticError(('Radial ladder failed',emax))
        narrow=audit(d,1.,56.,10,8,24,overlap=True)
        historical=np.array(archival[g]['scattering'][-1]['daughter_plus_photon'][:2])/2
        now=np.array([narrow[k] for k in ('gain','loss')])
        overlap_error=float(np.max(abs(now/historical-1.)))
        if overlap_error>.006:
            raise ArithmeticError(('Previous restricted crossing did not reproduce',historical,now,overlap_error))
        B=born(d)
        fixtures.append(dict(grid=g,Q=Q,geometry=controls,meshes=rows,shorter_cap42=tail,rigorous_tail_bound=bound,
                             equilibrium=therm,relative_last_mesh_drift=float(emax),
                             old_restricted_domain_current_one_charge=narrow,
                             old_restricted_domain_original_one_charge=historical.tolist(),
                             max_restricted_legacy_relative_difference=overlap_error,
                             restricted_net_fraction_of_full=narrow['net']/rows[-1]['net'],
                             one_charge_Born_net=B,relative_to_Born_at_g_point005=.005**2*rows[-1]['net']/B))
    return dict(scope='One fixed cold scalar mode, frozen charged and scalar archive, photon bath T=1, vacuum scalar absorption and its inverse. No thermal widths or virtual matching.',data=fixtures)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
