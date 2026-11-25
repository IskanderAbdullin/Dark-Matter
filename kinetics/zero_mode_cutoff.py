"""Zero-mode exchange with a regular isotropic gas at fixed positive mass.

Tree-level lambda phi^4/4!; one condensate leg and three gas legs, all
on shell. No formation prescription for a singular distribution.

In kinetic energy K=E-m, net d ln N0/dt = lambda^2/(128 pi^3 m)
times integral_{b,c>=k} [f(b)f(c)-f(b+c)(1+f(b)+f(c))] db dc.
For a regular integrable f and k=0 the quadratic terms cancel exactly
by exchanging the triangle integration order. A finite lower cutoff
leaves an explicitly positive strip, which can dominate at large f.
Evaluate that strip directly, not as a difference of two huge rates.
"""
import argparse,json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from quartic_consistency import shell_nodes,gap,lambda_from_sigma,MPL

def kinetic(p,m):
    return p*p/(np.hypot(p,m)+m)

def exchange(K,f,k):
    K=np.asarray(K,float);f=np.maximum(np.asarray(f,float),0.)
    assert len(K)==len(f)+1 and k>=0 and np.all(np.diff(K)>0)
    L=np.maximum(K[:-1],k);U=K[1:]
    active=(U>L)&(f>0)
    L,U,F=L[active],U[active],f[active]
    A=float(F@(U-L))
    lo=np.maximum(L,2*k); width=np.maximum(U-lo,0.)
    B=float(F@(width*(lo-2*k+width/2)))
    # The interval [a-k,a] intersect [L_j,U_j] has length
    # max(0,min(k,a-L_j,U_j+k-a,U_j-L_j)).
    # Its slope changes only at the four interval endpoints.
    strip=0.
    if k>0:
        for l,u,w in zip(L,U,F):
            keep=(L<u)&(U+k>l)
            ll,uu,ff=L[keep],U[keep],F[keep]
            if not len(ll):continue
            cuts=np.unique(np.clip(np.r_[l,u,ll,uu,ll+k,uu+k],l,u))
            mid=(cuts[1:]+cuts[:-1])/2
            overlap=np.maximum(0.,np.minimum(np.minimum(k,mid[:,None]-ll),
                                np.minimum((uu-mid[:,None])+k,uu-ll)))
            strip+=2*w*float(np.diff(cuts)@(overlap@ff))
    return dict(gain_quadratic=A*A,positive_cutoff_strip=strip,
        spontaneous_loss=B,net=strip-B)

def screen(row,sigma,targets=(.2,.5,80.)):
    e=np.array(row['frozen_cell_spectrum']['edges'])
    hist=row.get('history',[row['final']])
    ps=row.get('scalar_history_psd',[dict(x=row['final']['x'],f_phi=row['frozen_cell_spectrum']['f_phi'])])
    lam=float(lambda_from_sigma(sigma));out=[]
    for x0 in targets:
        i=int(np.argmin([abs(s['x']-x0) for s in ps]))
        s,h=ps[i],hist[i];x=s['x'];f=np.maximum(np.asarray(s['f_phi']),0.)
        q,w=shell_nodes(e,f)
        if w.sum()==0:continue
        m=gap(q,w,1e-12*x,lam);K=kinetic(e,m)
        TH=MPL*x/(1.66*np.sqrt(row['g_bath'])*row['parent_GeV'])*np.sqrt(np.pi**2/30*row['g_bath']/h['total_rho_Tref4'])
        fac=lam**2/(128*np.pi**3*m)*TH
        values=[]
        for pc in np.r_[0.,np.geomspace(1e-11,1e-3,33)]:
            v=exchange(K,f,kinetic(pc,m))
            Ncut=float(f@np.maximum(np.minimum(e[1:],pc)**3-e[:-1]**3,0)/(6*np.pi**2))
            v.update(q_cut=float(pc),removed_fraction=Ncut/w.sum(),net_over_H=fac*v['net'])
            values.append(v)
        out.append(dict(x=x,m_eff=m,lambda_quartic=lam,N=float(w.sum()),fmax=float(max(f)),values=values))
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--sigma',type=float,default=1e-7)
    p.add_argument('--output',required=True);a=p.parse_args()
    r=screen(json.loads(Path(a.input).read_text())[0],a.sigma)
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    for s in r:
        pos=[v for v in s['values'] if v['net']>0.]
        print('x',s['x'],'m',s['m_eff'],'zero cutoff',s['values'][0]['net_over_H'],
              'first positive',pos[0] if pos else None,'max positive',max(v['net_over_H'] for v in s['values']))
