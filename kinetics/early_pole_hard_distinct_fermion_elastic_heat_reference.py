"""Frozen signed heat in a vacuum hard t-channel between distinct fermions."""
import argparse
import json
import numpy as np
from scipy.special import expit,roots_legendre
from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_live_joint_hard_pair_packets_reference import gauss_each


def integrate(order, transfer, bounds, resolve_incident, fixtures, masses):
    nr,ni,nc,na=order
    names=('fine','third')
    profiles=[]
    for g in names:
        z=load(f'results/early-pole-spectrum-{g}-m06-2026-11-26.json.gz')
        profiles.append([(np.asarray(z[f'endpoint_{k}_momentum_edges']),
                          np.r_[0.,z[f'endpoint_{k}_occupations'],0.]) for k in ('parent','fermion')])
    low,high=bounds
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
    widths=np.diff(transfer)
    tau,tw=gauss_each(np.asarray(transfer),nc)
    out=np.zeros((len(widths),4 if fixtures else 2,3))
    def occupation(a,c):
        return c[1][np.searchsorted(c[0],a)]
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
            angles=np.concatenate([(edge1-vel[:,None])/w0[:,None],
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
            for k,basis in enumerate(profiles):
                fi,gi=fi0[k][h][i],gi0[k][h][i]
                fo,go=occupation(pa,basis[0]),occupation(pb,basis[1])
                gain=fi*gi*(1-fo)*(1-go)-fo*go*(1-fi)*(1-gi)
                flux=weights*width*gain
                out[zone,k]+=np.r_[np.dot(flux,delta),np.sum(abs(flux)),
                   np.dot(abs(flux),abs(delta))]
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
    return out


def run(level, bounds=(.45,8.), resolve_incident=1, fixtures=True, massive=False):
    order=[(2,4,5,3),(3,6,7,3),(4,8,9,3),(1,3,4,2)][level]
    transfer=[.10,.15,.25,.4,.7,1.]
    masses=(.05,.003) if massive else (0.,0.)
    A=integrate(order,transfer,bounds,resolve_incident,fixtures,masses)
    H=np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18
    factor=.005**4/H
    if fixtures:
        assert np.max(abs(A[:,2,0]))<5e-16 and np.min(A[:,3,0])>0.
    names=('fine','third','common_T','parent_hotter') if fixtures else ('fine','third')
    return dict(scope='Vacuum kinematics and spin kernel on reversible vacuum-hard '
                      'windows; no soft transport, no Compton, no physical rate calibration.',
                order=order, masses=masses, bounds=bounds, incident_edges=resolve_incident, tau_cm_edges=transfer,
                rows={k:dict(heat_to_daughter_per_reference_log_a=float(sum(v[1:,0])*factor),
                             heat_by_tau_band_per_reference_log_a=list(v[:,0]*factor),
                             absolute_number_flux_above_tau_015=float(sum(v[1:,1])),
                             absolute_of_interval_net_heat_above_tau_015=float(sum(v[1:,2])*factor))
                      for k,v in zip(names,A.transpose(1,0,2),strict=True)})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,default=0,choices=range(4))
    p.add_argument('--bounds',nargs=2,type=float,default=[.45,8.])
    p.add_argument('--incident-edges',type=int,choices=range(3),default=1)
    p.add_argument('--without-fixtures',action='store_true')
    p.add_argument('--massive',action='store_true')
    z=p.parse_args()
    print(json.dumps(run(z.level,z.bounds,z.incident_edges,not z.without_fixtures,z.massive),indent=2))
