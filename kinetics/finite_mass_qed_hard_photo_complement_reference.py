"""Frozen P+photon <-> b+scalar real scattering above a declared lab photon cut."""
import gzip,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit
from finite_mass_qed_scalar_absorption_reference import Q,azimuth,frame,gl,born

def audit(d,T,cap,rad,polar,inner,kmin,equilibrium=False,overlap=False):
    m=d['x_end'];M=m*d['mass_ratio'];ms=d['checkpoints'][-1]['scalar_pole_over_T']
    es=np.hypot(ms,Q)
    p_edges,b_edges=[np.hypot(a,d[f'endpoint_{n}_momentum_edges'])
                     for a,n in ((M,'fermion'),(m,'parent'))]
    fp,fb=[np.r_[0.,d[f'endpoint_{k}_occupations'],0.] for k in ('fermion','parent')]
    fs=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],Q)-1]
    if equilibrium: fs=1/np.expm1(es/T)
    za,wa=gl(rad);zm,wm=gl(polar);zk,wk=gl(inner)
    # For a fixed incoming daughter energy and lab photon energy, solve
    # the angular boundary exactly from (W^2-m^2-2k E_total)^2=4 k^2 P^2.
    polars=[-1.,1.]
    def angle_edges(E,k):
        q=Q
        P=np.sqrt(np.maximum(0.,(E-M)*(E+M)))
        Et=E+es
        pb2=(Et-k-m)*(Et-k+m)
        w0=M*M+ms*ms+2*es*E-m*m
        ok=(k>0.)&(pb2>0.)&(P>0.)
        v0=(w0[ok]/2-k[ok]*(Et[ok]-k[ok]))/(P[ok]*q)
        h=k[ok]*np.sqrt(pb2[ok])/(P[ok]*q)
        for roots in (v0-h,v0+h):
            polars.extend(roots[(roots>-1.)&(roots<1.)])
    if not equilibrium:
        for kc in ([kmin,4.] if overlap else [kmin]):
            ee=np.r_[p_edges,b_edges+kc-es]
            ee=ee[(ee>M)&(ee<cap)]
            angle_edges(ee,np.full_like(ee,kc))
        a0,b0=np.broadcast_arrays(p_edges[:,None],b_edges[None,:])
        kp=a0+es-b0
        ok=(a0>M)&(a0<cap)&(kp>kmin)
        angle_edges(a0[ok],kp[ok])
    polars=np.unique(polars)
    w0=np.diff(polars)/2
    v_all=((polars[:-1,None]+polars[1:,None])/2+w0[:,None]*zm).ravel()
    w_all=(w0[:,None]*wm).ravel()
    total=np.zeros(2)
    for v,wv in zip(v_all,w_all,strict=True):
        # J=0 is the exact photon opening. Here the incoming daughter is
        # much lighter than the outgoing parent; stay on the physical side.
        entrance=brentq(lambda E:frame(E,v,M,m,ms)[3],M,cap,xtol=1e-14)
        emin=entrance if not overlap else max(entrance,np.hypot(.45,M))
        opening=brentq(lambda E:frame(E,v,M,m,ms)[11]-kmin,entrance,cap,xtol=1e-14)
        emin=max(emin,opening)
        low_peak_root=None
        if frame(cap,v,M,m,ms)[10]>kmin:
            low_peak_root=brentq(lambda E:frame(E,v,M,m,ms)[10]-kmin,entrance,cap,xtol=1e-14)
        if overlap:
            cap=np.hypot(5.,M)
            bmin,bmax=np.hypot([.45,5.],m)
            event_edges=np.r_[b_edges,bmin,bmax]
            common_radial_boundaries=np.r_[event_edges+.5-es,event_edges+4.-es]
        else:
            event_edges=b_edges
            common_radial_boundaries=event_edges+kmin-es
        # k_max opens only O(k_min) above the exact 2->2 invariant
        # threshold. A global E grid never resolves this narrow radial peak.
        scale=max(opening-entrance,1e-12)
        edge_zoom=opening+scale*np.geomspace(.005,max(.01,(cap-opening)/scale),44)
        base=np.unique(np.clip(np.r_[emin,cap,np.geomspace(emin,cap,50),edge_zoom,common_radial_boundaries,
                                    [] if low_peak_root is None else [low_peak_root],
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
            ll=np.maximum(l,kmin)
            hh=np.maximum(ll,h)
            if overlap:
                ll=np.maximum.reduce([l,np.full_like(l,.5),et-bmax])
                hh=np.maximum(ll,np.minimum.reduce([h,np.full_like(h,4.),et-bmin]))
            rstar=np.hypot(Ppar,Pperp)
            rho=M*M/(rstar*(P0+rstar))
            pol=rho*np.array([.04,.15,.5,1.5,4.5,15.,45.,150.,450.,1500.])
            peak=np.where(Ppar>0,h-(h-l)/2*pol,l+(h-l)/2*pol)
            if equilibrium:
                edges=np.sort(np.clip(np.concatenate([ll,hh,peak],axis=1),ll,hh),axis=1)
            else:
                edges=np.sort(np.clip(np.concatenate([ll,hh,et-b_edges,peak],axis=1),ll,hh),axis=1)
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
            total+=wv*(ew[start:start+64]*clock)@np.column_stack([gain,loss])
    G,L=total
    N=G+fs*(G-L)
    return dict(gain=float(G),loss=float(L),net=float(N),net_over_positive_turnover=float(N/((1+fs)*G+fs*L)))


def tail_bound(d,cap):
    m=d['x_end'];M=m*d['mass_ratio'];ms=d['checkpoints'][-1]['scalar_pole_over_T']
    es=np.hypot(ms,Q)
    final_max=np.hypot(m,d['endpoint_parent_momentum_edges'][-1])
    initial_max=np.hypot(M,d['endpoint_fermion_momentum_edges'][-1])
    Smin=M*M+ms*ms+2*cap*(es-Q)
    if cap<=max(initial_max,final_max-es) or Smin<=m*m:
        raise ArithmeticError('tail assumptions failed')
    Amin=(1-m*m/Smin)*M*M/2
    Jmin=(Smin-m*m)/2
    D=(M*M+m*m-ms*ms)/2
    A0=4+4/Amin*(2*D*(D+M*m)/Jmin+2*D+2*M*m+D+ms*ms-m*m)
    A1=4/Amin*(es+Q)
    boltz=np.exp(-(cap+es-final_max))/(-np.expm1(-(cap+es-final_max)))
    return boltz*(A0*(cap+1)+A1*(cap*cap+2*cap+2))/(64*np.pi**3*es)


def run(grids=('fine','third')):
    root=Path(__file__).parent/'results'
    with gzip.open(root/'early-pole-frozen-yukawa-real-photon-crossings-2026-11-27.json.gz') as f:
        old=json.load(f)['grid']
    results=[]
    for g in grids:
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        b=born(d)
        rows=[]
        for cut in (.00004,.00016,.00064,.00256,.01024,.04096,.5):
            x=audit(d,1.,58.,5,6,18,cut)
            y=audit(d,1.,58.,8,8,26,cut)
            drift=max(abs(x[k]/y[k]-1) for k in ('gain','loss','net'))
            rows.append(dict(lab_photon_minimum_energy=cut,fine=y,
                             net_relative_to_Born_at_g_point005=.005**2*y['net']/b,
                             refinement_relative_drift=drift))
        own=audit(d,1.,58.,8,8,26,.5,overlap=True)
        previous=np.array(old[g]['scattering'][-1]['parent_plus_photon'][:2])/2
        diff=max(abs(own[k]/p-1) for k,p in zip(('gain','loss'),previous,strict=True))
        if diff>.009:raise ArithmeticError((g,'legacy overlap',own,previous,diff))
        eq=audit(d,1.,58.,6,6,32,.00016,equilibrium=True)
        if abs(eq['net_over_positive_turnover'])>2e-12:
            raise ArithmeticError((g,'balance',eq))
        results.append(dict(grid=g,rows=rows,unscaled_gain_absolute_tail_bound=float(tail_bound(d,58.)),old_window_legacy_prior=previous.tolist(),old_window_replay=own,
                            old_window_relative_difference=float(diff),equilibrium=eq))
    return dict(scope='Real photon absorption/reciprocal scattering, fixed q, only above a bath-frame k cut; no inclusive subtraction or virtual term.',data=results)


if __name__=='__main__':
    import sys
    print(json.dumps(run(tuple(sys.argv[1:]) or ('fine','third')),indent=2,allow_nan=False))
