"""CM coordinates for a frozen two-flavor elastic heat and rest-frame matching comparison."""
import argparse
import json
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import expit,roots_legendre
from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_frozen_screened_elastic_amplitude_reference import matrix_trace,htl_massless,tabulated_polarization,vacuum_trace
from early_pole_live_joint_hard_pair_packets_reference import gauss_each


def integrate(order, transfer, bounds, resolve_incident, fixtures, masses, screening, smear, occupation_model):
    nr,ni,nc,na=order
    names=('fine','third')
    profiles=[]
    for g in names:
        z=load(f'results/early-pole-spectrum-{g}-m06-2026-11-26.json.gz')
        profiles.append([(np.asarray(z[f'endpoint_{k}_momentum_edges']),
                          np.r_[0.,z[f'endpoint_{k}_occupations'],0.]) for k in ('parent','fermion')])
    low,high=bounds
    curves=[]
    if occupation_model!='steps':
        for pair in profiles:
            one=[]
            for edge, F in pair:
                mid=(edge[:-1]+edge[1:])/2
                if occupation_model=='pchip':
                    one.append(PchipInterpolator(np.r_[edge[0],mid,edge[-1]],
                                  np.r_[F[1],F[1:-1],F[-2]],extrapolate=False))
                else:
                    x=edge**3/3
                    y=np.r_[0.,np.cumsum(np.diff(x)*F[1:-1])]
                    cdf=PchipInterpolator(x,y,extrapolate=False)
                    fn=cdf.derivative()
                    checks=fn(np.linspace(0,1,13)[:,None]*np.diff(x)+x[:-1])
                    assert np.nanmin(checks)>-1e-12 and np.nanmax(checks)<1.
                    one.append(lambda p,fn=fn:np.maximum(fn(p**3/3),0.))
            curves.append(one)
    m1,m2=masses
    if not 0. < low < high or m1+m2 >= 1.:
        raise ValueError("Invalid hard support or masses for s >= 1")
    m12,m22=m1*m1,m2*m2
    ms=m12+m22
    knots=np.r_[.45,.65,.8,1.,1.25,1.32,1.56,1.7,2.2,2.8,3.7,5.,6.5,8.,10.]
    def radial(species):
        breaks=list(knots)
        if resolve_incident:
            for state in profiles:
                edges,occup=state[species]
                # Resolve the sharp initial loss/inversion near the threshold. Full
                # incident smooth-radial convergence still tested by the order sweep.
                rapid=(edges>=1.28)&(edges<=1.59)&(abs(np.diff(occup))>.012)
                breaks.extend(edges if resolve_incident==2 else edges[rapid])
        return np.unique(np.r_[low,high,[i for i in breaks if low<i<high]])
    p,dp=gauss_each(radial(0),nr)
    b=radial(1)
    zp,wp=roots_legendre(nr)
    z,w=roots_legendre(ni)
    # Remove the kink at the exact incident s>=1 boundary for each incident p.
    ea=np.hypot(p,m1)
    K=(1.-ms)/2
    # Stable root of ea*sqrt(q²+m2²)+p*q = K for s_max = 1.
    floor=(K*K-ea*ea*m22)/(ea*np.sqrt(K*K-m12*m22)+p*K)
    left=np.maximum(b[:-1][None,:],floor[:,None])
    half=np.maximum(0.,(b[1:] - left)/2)
    B=left[...,None]+half[...,None]*(1+zp)
    I=(2*dp*p*p/np.pi**2)[:,None,None]
    J=2*half[...,None]*wp*B*B/np.pi**2
    P=np.broadcast_to(p[:,None,None],B.shape)
    EA,EB=np.hypot(P,m1),np.hypot(B,m2)
    smin=np.maximum(1.,ms+2*(EA*EB-P*B))
    ds=np.maximum(0.,ms+2*(EA*EB+P*B)-smin)
    ss=smin[...,None]+.5*ds[...,None]*(1+z)
    lam=np.maximum((ss-ms)**2-4*m12*m22,0.)
    # With the isotropic angle and reciprocal-pair factor, dcos/4=ds/(8 pq).
    rate=np.sqrt(lam)/(2*(EA*EB)[...,None])
    W=(I*J*ds/(16*P*B))[...,None]*rate*w
    take=W>0.
    q1,q2=(np.broadcast_to(p[...,None],W.shape)[take] for p in (P,B))
    e1,e2=(np.hypot(q,m) for q,m in zip((q1,q2),masses,strict=True))
    s=ss[take]
    lam=lam[take]
    W=W[take]
    E=e1+e2
    boost=np.sqrt(np.maximum(E*E-s,0.))
    ecm=(s+m12-m22)/(2*s)
    bcm=np.sqrt(lam)/(2*s)
    A=np.clip(np.divide(e1-E*ecm,boost*bcm,out=np.zeros_like(s),where=boost>0.),-1.,1.)
    phi,pw=roots_legendre(na)
    cz,cw=roots_legendre(nc)
    Edges=np.unique(np.r_[low, high,*[p[0][(p[0]>low)&(p[0]<high)]
                                      for x in profiles for p in x]])
    edge1,edge2=(np.hypot(Edges,m) for m in masses)
    Qcuts=np.array([.005,.01,.02,.04,.08,.16,.32,.64,1.28,2.56,5.12])
    pol=[tabulated_polarization(g,smear=smear) for g in names]
    mD=[x[1] for x in pol]
    functions=[x[0] if screening=='snapshot' else htl_massless(x[1]) for x in pol]
    widths=np.diff(transfer)
    logs,tw=gauss_each(np.log(transfer),nc)
    tau=np.exp(logs)
    tw*=tau
    out=np.zeros((len(widths),4 if fixtures else 2,3))
    Qhist=np.zeros((len(widths),2,len(Qcuts)+1))
    corr=np.zeros_like(Qhist)
    def occupation(a,c):
        return c[1][np.searchsorted(c[0],a)]
    if curves:
        fi0,gi0=([np.nan_to_num(curve[k](q)) for curve in curves]
                 for k,q in enumerate((q1,q2)))
    else:
        fi0,gi0=([occupation(q,species[k]) for species in profiles]
                 for k,q in enumerate((q1,q2)))
    for jt,(ta,wt) in enumerate(zip(tau,tw,strict=True)):
        zone=jt//nc
        # Fixed bands in tau_cm = |t|/(4 p_cm^2). The massless cut is |t|/s.
        ct=1.-2*ta
        t=-lam/s*ta
        u=ms-s-t
        spin=((s-ms)**2+u**2+2*t*ms)/(16*np.pi*s*t*t)
        v=E*ecm+boost*bcm*A*ct
        r=boost*bcm*np.sqrt((1-ct*ct)*(1-A*A))
        for left in range(0,len(E),224):
            h=slice(left,left+224)
            e=E[h];vel=v[h];w0=np.maximum(r[h],1e-30)
            transferE=np.sqrt(np.maximum(Qcuts**2+t[h,None],0.))
            Qpoints=np.concatenate([e1[h,None]-transferE,e1[h,None]+transferE],axis=1)
            angles=np.concatenate([(edge1-vel[:,None])/w0[:,None],
              (Qpoints-vel[:,None])/w0[:,None],
              (e[:,None]-edge2-vel[:,None])/w0[:,None],
              np.broadcast_to([-1.,1.],(len(e),2))],axis=1)
            p=np.sort(np.arccos(np.clip(angles,-1.,1.)),axis=1)
            d=np.diff(p,axis=1)
            i,j=np.nonzero(d>1.e-14)
            z0,z1=p[i,j],p[i,j+1]
            mid=(z0+z1)/2
            ea=vel[i]+w0[i]*np.cos(mid)
            eb=e[i]-ea
            idx=(ea>=edge1[0])&(ea<=edge1[-1])&(eb>=edge2[0])&(eb<=edge2[-1])
            pa=np.sqrt(np.maximum(0.,(ea-m1)*(ea+m1)))
            pb=np.sqrt(np.maximum(0.,(eb-m2)*(eb+m2)))
            i,z0,z1,mid,pa,pb=(a[idx] for a in (i,z0,z1,mid,pa,pb))
            width=z1-z0
            weights=2*wt*(W*spin)[h][i]/np.pi
            # Piecewise-constant FD gain by interval, analytic first azimuthal moment.
            delta=e[i]-vel[i]-e2[h][i]-w0[i]*np.cos(mid)*np.sinc(width/(2*np.pi))
            q=np.sqrt((np.hypot(pa,m1)-e1[h][i])**2-t[h][i])
            ib=np.searchsorted(Qcuts,q)
            soft=ib<=6
            jsoft=i[soft]
            qh=e1[h][jsoft,None]
            rs=s[h][jsoft,None]
            rt=t[h][jsoft,None]
            if curves:
                EN=vel[i,None]+w0[i,None]*np.cos(mid[:,None]+width[:,None]*phi/2)
                OE=e[i,None]-EN
                p_nodes=np.sqrt(np.maximum(0.,(EN-m1)*(EN+m1)))
                b_nodes=np.sqrt(np.maximum(0.,(OE-m2)*(OE+m2)))
                dE=OE-e2[h][i,None]
                omega=EN[soft]-qh
            else:
                omega=vel[jsoft,None]+w0[jsoft,None]*np.cos(mid[soft,None]+width[soft,None]*phi/2)-qh
            projection=vacuum_trace(rs,rt,ms)
            for_check=matrix_trace(rs,rt,qh,e2[h][jsoft,None],omega,masses)
            if for_check.size:
                assert np.max(abs(for_check/projection-1.))<2e-7
            for k,basis in enumerate(profiles):
                fi,gi=fi0[k][h][i],gi0[k][h][i]
                amp=matrix_trace(rs,rt,qh,e2[h][jsoft,None],omega,masses,functions[k])
                if curves:
                    fo=np.nan_to_num(curves[k][0](p_nodes))
                    go=np.nan_to_num(curves[k][1](b_nodes))
                    fi,gi=fi[:,None],gi[:,None]
                    gain=fi*gi*(1-fo)*(1-go)-fo*go*(1-fi)*(1-gi)
                    flux=(weights*width)[:,None]*gain*(pw/2)
                    heat=flux*dE
                    f=heat.sum(axis=1)
                    replace=np.sum(heat[soft]*(amp/projection-1),axis=1)
                    absflux=np.sum(abs(flux))
                else:
                    fo,go=occupation(pa,basis[0]),occupation(pb,basis[1])
                    gain=fi*gi*(1-fo)*(1-go)-fo*go*(1-fi)*(1-gi)
                    flux=weights*width*gain
                    f=flux*delta
                    replace=flux[soft]*np.sum(-omega*(pw/2)*(amp/projection-1),axis=1)
                    absflux=np.sum(abs(flux))
                Qhist[zone,k]+=np.bincount(ib,weights=f,minlength=len(Qcuts)+1)
                corr[zone,k]+=np.bincount(ib[soft],weights=replace,minlength=len(Qcuts)+1)
                out[zone,k]+=np.r_[f.sum(),absflux,np.sum(abs(f))]
            if not fixtures:
                continue
            q=e2[h][i,None]
            pa=vel[i,None]+w0[i,None]*np.cos(mid[:,None]+width[:,None]*phi/2)
            pb=e[i,None]-pa
            wk=weights[:,None]*width[:,None]*pw/2
            for k,T in enumerate((1.,.65),start=2):
                fi,gi=expit(-e1[h][i,None]+.1),expit(-q/T-.4)
                fo,go=expit(-pa+.1),expit(-pb/T-.4)
                gain=wk*(fi*gi*(1-fo)*(1-go)-fo*go*(1-fi)*(1-gi))
                out[zone,k]+=np.r_[np.sum(gain*(pb-q)),np.sum(abs(gain)),
                   np.sum(abs(gain*(pb-q)))]
    np.testing.assert_allclose(Qhist.sum(axis=2),out[:,:2,0],rtol=2e-10,atol=1e-15)
    return out,Qhist,corr,Qcuts,mD


def run(level, bounds=(.45,8.), resolve_incident=1, fixtures=False, massive=False,
        screening='snapshot', smear=.12, lower=1e-7, occupation_model='steps'):
    order=[(2,4,5,3),(3,6,7,3),(4,8,9,3),(1,3,4,2)][level]
    transfer=[lower,*[t for t in (1e-7,1e-6,1e-5,1e-4,.001,.01,.03,.10,.15,.25,.4,.7,1.) if t>lower]]
    masses=(.05,.003) if massive else (0.,0.)
    A,Q,C,edges,mD=integrate(order,transfer,bounds,resolve_incident,fixtures,masses,screening,smear,occupation_model)
    strict=transfer.index(.15)
    H=np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18
    factor=.005**4/H
    if fixtures:
        assert np.max(abs(A[:,2,0]))<5e-16 and np.min(A[:,3,0])>0.
    names=('fine','third','common_T','parent_hotter') if fixtures else ('fine','third')
    return dict(scope='Frozen Born distinct-flavor t-channel with a conditional soft Coulomb-gauge '
                      'propagator substitution and fixed hard state; no complete transport rate.',
                order=order, masses=masses, bounds=bounds, incident_edges=resolve_incident, tau_cm_edges=transfer,
                Q_cuts=edges.tolist(), Q_heat_by_tau_and_archive=(Q*factor).tolist(),
                screening_model=screening, occupation_model=occupation_model, response_smear=smear, mD2=mD,
                Q_screening_shift_by_tau_and_archive=(C*factor).tolist(),
                matched_Q_cuts=edges[:7].tolist(),
                matched_heat_per_reference_log_a=((Q.sum((0,2))[:,None]+C.sum(0).cumsum(1)[:,:7])*factor).tolist(),
                rows={k:dict(heat_to_daughter_per_reference_log_a=float(sum(v[strict:,0])*factor),
                             heat_by_tau_band_per_reference_log_a=list(v[:,0]*factor),
                             absolute_number_flux_above_tau_015=float(sum(v[strict:,1])),
                             absolute_of_interval_net_heat_above_tau_015=float(sum(v[strict:,2])*factor))
                      for k,v in zip(names,A.transpose(1,0,2),strict=True)})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,default=0,choices=range(4))
    p.add_argument('--bounds',nargs=2,type=float,default=[.45,8.])
    p.add_argument('--incident-edges',type=int,choices=range(3),default=1)
    p.add_argument('--without-fixtures',action='store_true')
    p.add_argument('--massive',action='store_true')
    p.add_argument('--screening',choices=('snapshot','massless'),default='snapshot')
    p.add_argument('--smear',type=float,default=.12)
    p.add_argument('--lower',type=float,default=1e-7,choices=(1e-9,1e-8,1e-7,1e-6,1e-5))
    p.add_argument('--occupation-model',choices=('steps','pchip','cell-conservative'),default='steps')
    z=p.parse_args()
    print(json.dumps(run(z.level,z.bounds,z.incident_edges,not z.without_fixtures,z.massive,z.screening,z.smear,z.lower,z.occupation_model),indent=2))
