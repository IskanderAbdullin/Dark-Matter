"""Instantaneous isotropic 2<->2 quantum elastic moment diagnostic.

L_int=-lambda phi^4/4!, one real scalar, fixed quasiparticle mass. The
angular phase-space measure is exact for a constant matrix element. All
momenta are in Tref units. Integrate in kinetic energy to avoid subtracting
nearly equal rest energies. Occupations are the original piecewise-constant
FV PSD; outgoing integration splits at BOTH daughter cell boundaries.

No time evolution, no 2<->4, no vertex resummation, no condensate field.
Large calculated moment rates diagnose the collisionless model's domain of
validity; they are not a solved interacting abundance or thermalization time.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from quartic_consistency import shell_nodes,lambda_from_sigma,gap,variance
from postproduction_elastic_screen import pair_factor

def kinetic(q,m):
    return q*q/(np.hypot(q,m)+m)

def momentum(K,m):
    return np.sqrt(np.maximum(K*(K+2*m),0.))

def nodes(edges,f,m,order,cut=.1):
    e=kinetic(np.asarray(edges),m);c=kinetic(cut,m)
    e2=np.sort(np.unique(np.r_[e,c])) if e[0]<c<e[-1] else e
    mid=(e2[1:]+e2[:-1])/2
    idx=np.clip(np.searchsorted(e,mid)-1,0,len(f)-1)
    u,w=roots_legendre(order)
    K=mid[:,None]+np.diff(e2)[:,None]*u/2
    wt=np.diff(e2)[:,None]*w/2
    F=np.broadcast_to(np.maximum(f,0)[idx,None],K.shape)
    mask=F>0
    return e,K[mask],wt[mask],F[mask]

def lookup(e,f,K):
    i=np.searchsorted(e,K,side='right')-1
    valid=(i>=0)&(i<len(f))
    return np.where(valid,f[np.clip(i,0,len(f)-1)],0.)

def integrals(edges,f,m,order_in=2,order_out=2,cut=.1,block=128):
    e,K,wt,F=nodes(edges,f,m,order_in,cut)
    f=np.maximum(f,0);p=momentum(K,m);c=kinetic(cut,m)
    # Ordered incoming pairs with symmetry factor 2 off the diagonal.
    i,j=np.triu_indices(len(K))
    result=np.zeros(5);absnet=np.zeros(2);start=time.monotonic()
    u,w=roots_legendre(order_out)
    for k in range(0,len(i),block):
        a,b=i[k:k+block],j[k:k+block]
        S=K[a]+K[b];S2=S/2
        pref=(wt[a]*wt[b]*F[a]*F[b]*np.where(a==b,1.,2.))[:,None]
        # Integrate only K3<=K4, multiply by 2. Include hot-cut boundaries.
        bps=np.r_[e,c]
        edges3=np.concatenate((np.zeros((len(S),1)),S2[:,None],
            np.broadcast_to(bps,(len(S),len(bps))),S[:,None]-bps[None,:]),axis=1)
        edges3=np.sort(np.clip(edges3,0,S2[:,None]),axis=1)
        L,H=edges3[:,:-1],edges3[:,1:];d=H-L;mid=(L+H)*.5
        p1=p[a,None];p2=p[b,None]
        lo12=np.abs(p1-p2);hi12=p1+p2
        incoming_q2=p1*p1+p2*p2
        incoming_hot=(K[a,None]>c).astype(float)+(K[b,None]>c).astype(float)
        for xi,wi in zip(u,w):
            K3=mid+d*xi/2;K4=S[:,None]-K3
            p3=momentum(K3,m);p4=momentum(K4,m)
            angle=np.maximum(np.minimum(hi12,p3+p4)-np.maximum(lo12,np.abs(p3-p4)),0.)
            f3=lookup(e,f,K3);f4=lookup(e,f,K4)
            base=pref*d*wi*angle  # outgoing symmetry 2 cancels GL interval 1/2
            # The f1*f2*f3*f4 term cancels EXACTLY under swapping incoming
            # and outgoing pairs in the moment equation. Remove analytically,
            # not numerically; retain spontaneous + cubic Bose terms.
            bos=1+f3+f4
            delta_q2=p3*p3+p4*p4-incoming_q2
            delta_hot=(K3>c).astype(float)+(K4>c).astype(float)-incoming_hot
            result[0]+=base.sum()
            result[1]+=np.sum(base*bos*delta_q2)
            result[2]+=np.sum(base*bos*delta_hot)
            result[3]+=np.sum(base*delta_q2)
            result[4]+=np.sum(base*delta_hot)
            absnet[0]+=np.sum(base*bos*np.abs(delta_q2))
            absnet[1]+=np.sum(base*bos*np.abs(delta_hot))
    q,wq=shell_nodes(np.asarray(edges),f,8)
    N=float(wq.sum());M2=float(wq@(q*q));V=variance(q,wq,m)
    # Check normalization by the independent dilute Moller-rate integral.
    # n*Gamma = lambda^2/(512*pi^5) * I0 = lambda^2/(64*pi) * pair_factor
    ref=8*np.pi**4*pair_factor(q,wq,m,16)
    return dict(incoming_nodes=len(K),incoming_pairs=len(i),wall_seconds=time.monotonic()-start,
        dilute_integral=float(result[0]),dilute_integral_reference=float(ref),
        dilute_normalization_rel_error=float(result[0]/ref-1),
        number=N,q2_moment=M2,variance=V,
        raw_q2_net=float(result[1]),raw_hot_net=float(result[2]),
        raw_q2_classical=float(result[3]),raw_hot_classical=float(result[4]),
        q2_net_over_absolute=float(result[1]/max(absnet[0],1e-300)),
        hot_net_over_absolute=float(result[2]/max(absnet[1],1e-300)))

def physical_rate(r,row,h,lam):
    x=h['x'];M=row['parent_GeV'];g=row['g_bath']
    Href=1.66*np.sqrt(g)*M*M/1.220890e19
    Hf=np.sqrt(np.pi**2/30*g/h['total_rho_Tref4'])
    factor=lam*lam/(1024*np.pi**5)*(M*x/Href)*Hf
    return dict(relative_q2_change_per_H=factor*r['raw_q2_net']/r['q2_moment'],
        absolute_hot_fraction_change_per_H=factor*r['raw_hot_net']/r['number'],
        relative_q2_classical_change_per_H=factor*r['raw_q2_classical']/r['q2_moment'],
        absolute_hot_fraction_classical_change_per_H=factor*r['raw_hot_classical']/r['number'])

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--sigma',type=float,default=.001);p.add_argument('--x',type=float,default=80.)
    p.add_argument('--order',type=int,default=2);p.add_argument('--massless',action='store_true')
    a=p.parse_args();row=json.loads(Path(a.input).read_text())[0]
    if row.get('scalar_history_psd'):
        j=int(np.argmin([abs(s['x']-a.x) for s in row['scalar_history_psd']]))
        h=row['history'][j];f=np.asarray(row['scalar_history_psd'][j]['f_phi'])
    else:
        h=row['final'];f=np.asarray(row['frozen_cell_spectrum']['f_phi'])
    e=np.asarray(row['frozen_cell_spectrum']['edges']);q,w=shell_nodes(e,f,8)
    lam=float(lambda_from_sigma(a.sigma));mu0=row['scalar_ratio']*h['x']
    mu=0. if a.massless else gap(q,w,mu0,lam)
    out=dict(input=a.input,x=h['x'],n=row['n'],sigma_over_m_cm2_g=a.sigma,lambda_quartic=lam,
        mu_eff_over_Tref=mu,mu_over_peak=mu/h['peak_bin_q'],massless=a.massless,orders=[])
    for oi,oo in [(a.order,a.order),(a.order+1,a.order+1)]:
        r=integrals(e,f,mu,oi,oo);r.update(order_in=oi,order_out=oo,**physical_rate(r,row,h,lam))
        out['orders'].append(r);print(json.dumps(r),flush=True)
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
