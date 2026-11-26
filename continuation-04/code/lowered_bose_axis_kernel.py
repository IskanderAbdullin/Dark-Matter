"""Axisymmetric hard-sphere UU angular kernel; E=q², p=q_z."""
import argparse,json
import numpy as np
from scipy.special import ellipkm1
from scipy.integrate import quad
from scipy.stats import qmc

def kernel(E,p):
    A=E-p*p;valid=np.all(A>0,axis=0);a,b,c,d=np.sqrt(np.maximum(A,0))
    L1=(a-b)**2;U1=(a+b)**2;L2=(c-d)**2;U2=(c+d)**2
    l=np.minimum(L1,L2);L=np.maximum(L1,L2);U=np.minimum(U1,U2);u=np.maximum(U1,U2)
    den=(U-l)*(u-L);ok=valid&(U>L)&(den>0)
    m1=np.divide((L-l)*(u-U),den,out=np.ones_like(den),where=ok)
    ans=np.zeros_like(den);ans[ok]=2*ellipkm1(np.clip(m1[ok],0,1))/(np.pi**3*np.sqrt(den[ok]))
    return ans

def direct(r):
    a,b,c,d=r;lo=sorted([(a-b)**2,(c-d)**2]);hi=sorted([(a+b)**2,(c+d)**2]);l,L=lo;U,u=hi
    if U<=L:return 0.
    return quad(lambda t:2/np.sqrt((L+(U-L)*np.sin(t)**2-l)*(u-L-(U-L)*np.sin(t)**2)),0,np.pi/2,epsabs=1e-12,epsrel=1e-12)[0]/np.pi**3

def expected(E,i,j):
    q=np.sqrt(E);S=E[0]+E[1];a=max(abs(q[0]-q[1]),abs(q[2]-q[3]));b=min(q[0]+q[1],q[2]+q[3]);m2=(a*a+a*b+b*b)/3;m4=(a**4+a**3*b+a*a*b*b+a*b**3+b**4)/5;inv2=1/(a*b) if a else 0.;inv4=m2/(a**3*b**3) if a else 0.;d=E-E[[1,0,3,2]];ap2=(m2+2*d+d*d*inv2)/4
    if i==j:dot=E[i];dot2=E[i]**2
    elif i//2==j//2:dot=(m2-S)/2;dot2=(m4-2*S*m2+S*S)/4
    else:
        x,y=d[i],d[j];dot=(m2+x+y+x*y*inv2)/4;prod=(m4+2*(x+y)*m2+x*x+y*y+4*x*y+2*x*y*(x+y)*inv2+x*x*y*y*inv4)/16;dot2=1.5*prod+.5*(E[i]*E[j]-E[i]*ap2[j]-E[j]*ap2[i])
    return dot/3,(3*dot2-E[i]*E[j])/10

def marginal(E,m,seed):
    E=np.asarray(E,dtype=float);q=np.sqrt(E);x=qmc.Sobol(3,scramble=True,seed=seed).random_base2(m);p3=(2*x-1)*q[:3];p=np.c_[p3,p3[:,0]+p3[:,1]-p3[:,2]].T;I=kernel(E[:,None],p);vol=8*np.prod(q[:3]);normal=vol*I.mean();theory=2*q.min()/np.pi
    h=1.5*p*p-E[:,None]/2;out=dict(E=E.tolist(),m=m,seed=seed,normal=normal,theory=theory,relative=normal/theory-1,singular=int(np.sum(~np.isfinite(I))))
    for i,j in [(0,0),(0,1),(0,2),(0,3),(2,3)]:
        c1,c2=expected(E,i,j);out[f'cov{i}{j}']=[float(np.mean(I*p[i]*p[j])/I.mean()),float(np.mean(I*h[i]*h[j])/I.mean()),c1,c2]
    return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--m',type=int,default=18);p.add_argument('--seed',type=int,default=0);a=p.parse_args()
 for E in [[1,2,1.4,1.6],[.03,3,.9,2.13],[.5,.5,.5,.5],[.1,1.9,.3,1.7]]:print(json.dumps(marginal(E,a.m,a.seed)))
