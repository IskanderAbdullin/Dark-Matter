"""Late NR contact switch-on: instantaneous frozen-PSD moment diagnostic.
Physical rates normalized by n*sigma*v_peak; no evolved late PSD.
All energy variables y=(p/p_peak)^2. Full f^4 cancellation analytically.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from quantum_elastic_moments import lookup

def compute(edges,f,p0,order,block=96):
    e=(np.asarray(edges)/p0)**2
    f=np.maximum(f,0.)
    cuts=np.array([2.,6.,20.])**2
    e2=np.unique(np.r_[e,cuts]);e2=e2[(e2>=e[0])&(e2<=e[-1])]
    mid=(e2[1:]+e2[:-1])/2
    idx=np.clip(np.searchsorted(e,mid)-1,0,len(f)-1)
    u,w=roots_legendre(order)
    K=mid[:,None]+np.diff(e2)[:,None]*u/2
    wt=np.diff(e2)[:,None]*w/2
    F=np.broadcast_to(f[idx,None],K.shape)
    sel=F>0
    K,wt,F=K[sel],wt[sel],F[sel]
    N0=float(np.sum(wt*F*np.sqrt(K)))
    moments=np.r_[np.sum(wt*F*np.sqrt(K)*K*K),[np.sum(wt*F*np.sqrt(K)*(K>c)) for c in cuts]]
    i,j=np.triu_indices(len(K))
    net=np.zeros(4);absolute=np.zeros(4);classical=np.zeros(4);normal=0.
    for k in range(0,len(i),block):
        a,b=i[k:k+block],j[k:k+block]
        S=K[a]+K[b];S2=S/2
        pref=(wt[a]*wt[b]*F[a]*F[b]*np.where(a==b,1.,2.))[:,None]
        bounds=np.r_[e,cuts]
        E=np.concatenate([np.zeros((len(S),1)),S2[:,None],
                 np.broadcast_to(bounds,(len(S),len(bounds))),S[:,None]-bounds[None,:]],axis=1)
        E=np.sort(np.clip(E,0,S2[:,None]),axis=1)
        L,H=E[:,:-1],E[:,1:];d=H-L;mid=(L+H)/2
        in1=K[a,None];in2=K[b,None]
        for xi,wi in zip(u,w):
            y3=mid+d*xi/2;y4=S[:,None]-y3
            angle=2*np.sqrt(np.maximum(np.minimum(np.minimum(in1,in2),np.minimum(y3,y4)),0.))
            F3=lookup(e,f,y3);F4=lookup(e,f,y4)
            base=pref*d*wi*angle
            bos=1+F3+F4
            ds=[y3*y3+y4*y4-in1*in1-in2*in2]
            for c in cuts:ds.append((y3>c).astype(float)+(y4>c).astype(float)-(in1>c).astype(float)-(in2>c).astype(float))
            for ix,dy in enumerate(ds):
                net[ix]+=np.sum(base*bos*dy)
                absolute[ix]+=np.sum(base*bos*abs(dy))
                classical[ix]+=np.sum(base*dy)
            normal+=float(base.sum())
    # C=lambda^2/(1024*pi^5); dK^3=(p_peak^2/(2m))^3.
    # n sigma v_peak=lambda^2*p_peak^4*M0/(512*pi^3*m^3).
    rate=net/(4*N0*moments)
    return dict(order=order,number_integral=N0,normalized_moments=(moments/N0).tolist(),
                relative_rate_over_n_sigma_vpeak=rate.tolist(),
                classical_relative_rate_over_n_sigma_vpeak=(classical/(4*N0*moments)).tolist(),
                net_over_absolute=(net/np.maximum(absolute,1e-300)).tolist(),
                classical_pair_norm=normal,labels=['p4','N(q>2peak)','N(q>6peak)','N(q>20peak)'])

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    row=json.loads(Path(a.input).read_text())[0];h=row['final'];p0=h['peak_bin_q']
    e=row['frozen_cell_spectrum']['edges'];f=row['frozen_cell_spectrum']['f_phi']
    rs=[compute(e,f,p0,k) for k in (1,2,3)]
    # Tref today inferred from conserved SM entropy and final bath theta.
    Tref0=2.7255*8.617333262e-5*(3.909/row['g_bath'])**(1/3)/h['theta']
    v0=p0*Tref0/86*2.99792458e10
    h0=.674;Om=.315;Or=9.2e-5;H0=h0*100*1e5/3.085677581e24
    rho0=.12*1.8783416e-29
    late=[]
    for z in (0.,1.,3.,10.,30.,100.):
        H=H0*np.sqrt(Om*(1+z)**3+Or*(1+z)**4+1-Om-Or)
        R=rho0*v0*(1+z)**4/H
        rr=np.asarray(rs[-1]['relative_rate_over_n_sigma_vpeak'])*R
        late.append(dict(z=z,vpeak_cm_s=v0*(1+z),n_sigma_v_over_H_per_cm2_g=R,
                         relative_rates_over_H_per_cm2_g=rr.tolist()))
    result=dict(input=a.input,Tref0_eV=Tref0,peak_q=p0,peak_velocity_today_cm_s=v0,
                rho0_g_cm3=rho0,orders=rs,late=late,scope='instantaneous NR frozen-PSD contact collision moments; turn-on not evolved, tree kinetics may fail')
    Path(a.output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
