"""Split the actual incoming and outgoing hard support of an archived pair fold."""
import argparse
import json
import numpy as np
from scipy.special import roots_legendre
from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_live_joint_hard_pair_packets_reference import gauss_each


def integrate(snapshots,species,orders):
    nr,ni,nf,na=orders
    # Resolve the incident cut at 2, not with an indicator on an unsplit mesh.
    incident=np.r_[.5,.65,.85,1.,1.25,1.55,2.,2.5,3.3,4.,4.5,6.,8.,8.5,12.,17.,25.,40.]
    p,w=gauss_each(incident,nr)
    V=w*p*p/np.pi**2
    ia,ib=np.triu_indices(len(p))
    x,y=p[ia],p[ib]
    weight=V[ia]*V[ib]*(1+(ia!=ib).astype(float))
    rel,relw=roots_legendre(ni)
    cut=np.minimum(1.,1/(4*x*y))
    a=cut[:,None]+.5*(1-cut[:,None])*(1+rel)
    s=4*x[:,None]*y[:,None]*a
    shape=s.shape
    s=s.ravel()
    W=(weight[:,None]*(1-cut[:,None])*a*relw).ravel()
    xv,yv=[np.broadcast_to(k[:,None],shape).ravel() for k in (x,y)]
    Ec=xv+yv
    keep=W>0.
    s,W,xv,yv,Ec=[k[keep] for k in (s,W,xv,yv,Ec)]
    boost=np.sqrt(np.maximum(Ec*Ec-s,0.))
    rho=np.clip((xv-yv)/boost,-1.,1.)
    cos,cw=roots_legendre(nf)
    za,wa=roots_legendre(na)
    m=.05 if species=='parent' else .003
    beta=np.sqrt(1-4*m*m/s)
    M=np.hypot(np.unique(np.r_[.45,2.2,5.,10.,*[d[f'endpoint_{species}_momentum_edges']
                                        for d in snapshots]]),m)
    align=[]
    for d in snapshots:
        e=np.hypot(d[f'endpoint_{species}_momentum_edges'],m)
        f=np.asarray(d[f'endpoint_{species}_occupations'])
        align.append((e,np.r_[0.,f,0.]))
    g1,g2=1/np.expm1(xv),1/np.expm1(yv)
    # 1/2 identical incident particles times 2 outgoing count.
    result=np.zeros((len(snapshots),7,4))
    # Full hard, old active hard, photons only, outgoing fermions only.
    for z,wz in zip(.8*cos,.8*cw,strict=True):
        r=4*m*m/s
        u=1-(beta*z)**2
        bw=beta*(4*(1+r)/u-2-4*r*r/(u*u))/(16*np.pi*s)
        E0=(Ec+boost*beta*rho*z)/2
        A=boost*beta*np.sqrt((1-rho*rho)*(1-z*z))/2
        for left in range(0,len(E0),128):
            h=slice(left,left+128)
            total=Ec[h]
            b=E0[h]
            a=A[h]
            ph=np.concatenate([(M-b[:,None])/a[:,None],
                               (total[:,None]-M-b[:,None])/a[:,None],
                               np.broadcast_to([-1.,1.],(len(b),2))],axis=1)
            ph=np.sort(np.arccos(np.clip(ph,-1,1)),axis=1)
            wi=np.diff(ph,axis=1)
            ii,jj=np.nonzero(wi>1e-14)
            Z=(ph[ii,jj]+ph[ii,jj+1])/2
            v=wi[ii,jj]/2
            az=Z[:,None]+v[:,None]*za
            e1=b[ii,None]+a[ii,None]*np.cos(az)
            e2=total[ii,None]-e1
            assert min(e1.min(),e2.min())>=m
            isph=(xv[h][ii]<=2.)&(yv[h][ii]<=2.)
            isf=(e1>=np.hypot(.45,m))&(e2>=np.hypot(.45,m))&(
                e1<=np.hypot(2.2,m))&(e2<=np.hypot(2.2,m))
            w=(W[h]*bw[h]*wz)[ii,None]*(v[:,None]*wa)/np.pi
            p1,p2=g1[h][ii,None],g2[h][ii,None]
            mid=((xv[h][ii]<=4.)&(yv[h][ii]<=4.))[:,None]
            mid=mid&(e1>=np.hypot(.45,m))&(e2>=np.hypot(.45,m))&(
                e1<=np.hypot(5.,m))&(e2<=np.hypot(5.,m))
            wide=((xv[h][ii]<=8.)&(yv[h][ii]<=8.))[:,None]
            wide=wide&(e1>=np.hypot(.45,m))&(e2>=np.hypot(.45,m))&(
                e1<=np.hypot(10.,m))&(e2<=np.hypot(10.,m))
            low=(e1<np.hypot(.45,m))|(e2<np.hypot(.45,m))
            for k,(ee,ff) in enumerate(align):
                a1,a2=(ff[np.searchsorted(ee,e)] for e in (e1,e2))
                forward=w*p1*p2*(1-a1)*(1-a2)
                reverse=w*(1+p1)*(1+p2)*a1*a2
                net=forward-reverse
                heat=.5*net*total[ii,None]
                for j,mask in enumerate((1.,isph[:,None]&isf,isph[:,None],isf,mid,wide,low)):
                    result[k,j]+=np.array([np.sum(q*mask) for q in (net,heat,forward,reverse)])
    return result


def run(level):
    orders=[(2,4,10,2),(3,6,16,2),(4,8,22,3),(5,10,28,3)][level]
    grids=['fine','third']
    sn=[load(f'results/early-pole-spectrum-{g}-m06-2026-11-26.json.gz') for g in grids]
    names=['entire_selected_hard','active_photon_fermion','active_photons','active_fermions',
           'expanded_photon_4_fermion_5','expanded_photon_8_fermion_10','outgoing_fermion_below_.45']
    H=np.sqrt(np.pi**2*106.75/90)*(86/.05)/2.435e18
    data=[]
    for species in ('parent','fermion'):
        v=integrate(sn,species,orders)
        for g,R in zip(grids,v,strict=True):
            row={}
            for lab,r in zip(names,R,strict=True):
                row[lab]=dict(vacuum_hard_count_per_reference_log_a=r[0]*.005**4/H,
                              vacuum_hard_heat_per_reference_log_a=r[1]*.005**4/H,
                              two_moment_share_of_hard=(r[:2]/R[0,:2]).tolist(),
                              absolute_signed_input_weights=r.tolist())
            data.append(dict(grid=g,species=species,rows=row))
    return dict(scope='Frozen charged distributions, fixed unit-temperature BE photon bath, '
                      'vacuum central hard subset only. Not the cosmological rate, no soft or medium matching.',
                quadrature=orders,rows=data)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,default=0)
    a=p.parse_args()
    print(json.dumps(run(a.level),indent=2,allow_nan=False))
