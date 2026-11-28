"""Frozen isotropic massive Vlasov response and independent noise moments."""
from pathlib import Path
import argparse
import gzip
import json
import numpy as np
from scipy.special import roots_legendre

from early_pole_snapshot_hard_loop_screening_reference import C,raw

THERMAL=2/3


def real_kernels(s,v):
    z=np.divide(s,v,out=np.full(np.broadcast_shapes(np.shape(s),np.shape(v)),1e12),where=v>0)
    z=np.maximum(z,0.)
    inside=z<1.
    zz=np.where(z>30.,2.,z)
    h=.5*(np.log1p(zz)-np.log(np.maximum(abs(1-zz),1e-14)))
    L=1-zz*h
    T=.5*(zz*zz+zz*(1-zz*zz)*h)
    a=np.minimum(1/np.maximum(z,1.),1/30.)**2
    L=np.where(z>30.,-a*(1/3+a*(1/5+a/7)),L)
    T=np.where(z>30.,1/3+a*(1/15+a/35),T)
    return L,T,inside,z


class Snapshot:
    def __init__(self,d,order,smear=0.):
        self.data=[]
        m0=float(d['x_end'])
        y,w=roots_legendre(order)
        for lab,mass in (('parent',m0),('fermion',m0*float(d['mass_ratio']))):
            p=np.asarray(d[f'endpoint_{lab}_momentum_edges'])
            f=np.asarray(d[f'endpoint_{lab}_occupations'])
            h=np.diff(p)
            j=-np.diff(np.r_[0.,f,0.])
            if smear:
                a=smear*np.minimum(np.r_[p[0],h],np.r_[h,h[-1]])
                pr=(p[:,None]+a[:,None]*y).ravel()
                jr=(j[:,None]*(w/2)).ravel()
            else:
                pr,jr=p,j
            er=np.hypot(pr,mass)
            self.data.append((mass,p,f,pr/er,C*pr*er*jr))
        self.smear=smear
        self.y=y
        self.w=w

    def retarded(self,s):
        s=np.asarray(s)[...,None]
        ans=np.zeros(s.shape[:-1]+(4,))
        for _,p,f,v,d in self.data:
            L,T,inside,z=real_kernels(s,v)
            ImL=np.pi/2*np.where(inside,z,0.)
            ImT=np.pi/4*np.where(inside,z*(1-z*z),0.)
            ans+=(np.stack([L@d, T@(d*v*v), ImL@d, ImT@(d*v*v)],axis=-1))
        return ans

    def noise(self,s):
        ans=np.zeros(np.shape(s)+(2,))
        s=np.asarray(s)
        for m,p,f,_,_ in self.data:
            start=m*s/np.sqrt(np.maximum(1-s*s,1e-30))
            lo=np.maximum(p[:-1],start[...,None])
            hi=np.maximum(p[1:],lo)
            h=(hi-lo)[...,None]/2
            k=(hi+lo)[...,None]/2+h*self.y
            E=np.hypot(k,m);v=k/E
            nv=C*k*E*h*self.w*(f*(1-f))[:,None]
            ans[...,0]+=np.pi/2*np.sum(nv,axis=(-1,-2))
            ans[...,1]+=np.pi/4*np.sum(nv*np.maximum(v*v-s[...,None,None]**2,0.),
                                       axis=(-1,-2))
        return ans

    def evaluate(self,s):
        y=self.retarded(s)
        p=self.noise(s)
        return y,p

    def plasma_frequency_sq(self):
        return sum((d*v*v).sum() for _,_,_,v,d in self.data)/3


def load(name):
    return json.loads(gzip.decompress(Path(name).read_bytes()))


def run():
    rows=[]
    point=np.array([0.,.01,.1,.5,.9,.99,.999,1.0001,1.1,2.,100.])
    for grid in ('fine','third'):
        d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
        for f in (0.,.12,.28):
            snap=Snapshot(d,18,f)
            R,N=snap.evaluate(point)
            assert abs(R[-1,1]/snap.plasma_frequency_sq()-1)<5e-5
            a=[raw(p,h,m,0.)['static_mD2_over_g2'] for m,p,h,_,_ in snap.data]
            if f==0.:
                assert abs(R[0,0]-sum(a))<1e-13
            rows.append(dict(grid=grid,smear_fraction=f,
                             response=R.tolist(),noise_weights=N.tolist(),
                             plasma_frequency_sq_over_g2=snap.plasma_frequency_sq()))
    return dict(scope='Frozen isotropic collisionless Abelian snapshot. Noise is independent '
                      'fermion occupation variance, not KMS for the mixed state.',
                phase_velocity=point.tolist(),
                retarded_fields=['Coulomb_real','transverse_real','Coulomb_abs_imag','transverse_abs_imag'],
                statistical_fields=['Coulomb_sym_omega_divided','transverse_sym_omega_divided'],
                rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
