"""Split the outgoing angular continuum at all receiving-node boundaries."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.special import roots_legendre

from early_pole_scalar_number_change_angular_reference import data
from early_pole_scalar_number_change_stencil_reference import NumberChange
from early_pole_scalar_number_change_threshold_reference import LAMBDA, threshold


class IntervalNumberChange(NumberChange):
    def __init__(self,q,edges,m,qlo,amp_scale,power=3,angular_order=2):
        E=np.hypot(q,m)
        vol=np.diff(edges**3)/(6*np.pi**2)
        cold=np.flatnonzero((edges[:-1]>=qlo)&(edges[1:]<=1e-5))
        I,mult,lab,cos=data(tuple(cold),tuple(q),power)
        total=E[I].sum(axis=1)
        shat=total[:,None]**2-lab**2
        beta=np.sqrt(1-4*m*m/shat)
        span=.5*lab*beta
        assert min(span.ravel())>0.
        mid=total/2
        imin=mid[:,None]-span
        imax=mid[:,None]+span
        j0=np.searchsorted(E,imin.min(),side='right')-1
        j1=np.searchsorted(E,imax.max(),side='left')
        egrid=E[j0:j1+1]
        raw=np.sort(np.c_[np.broadcast_to(egrid,(len(I),len(egrid))),total[:,None]-egrid],axis=1)
        lower=np.maximum(imin[:,:,None],raw[:,None,:-1])
        upper=np.minimum(imax[:,:,None],raw[:,None,1:])
        ix,ig,ip=np.nonzero(upper>lower)
        top,bot=upper[ix,ig,ip],lower[ix,ig,ip]
        factor=mult[ix]*np.prod(vol[I[ix]]*m/E[I[ix]],axis=1)
        factor*=beta[ix,ig]*(top-bot)/(len(cos)*np.sqrt(.75)*2*span[ix,ig])
        offset,w=roots_legendre(angular_order)
        oxy=(top[:,None]+bot[:,None])/2+(top-bot)[:,None]*offset/2
        w0=(factor[:,None]*w/2).ravel()
        a=oxy.ravel()
        ix=np.repeat(ix,angular_order)
        b=total[ix]-a
        hi=np.searchsorted(E,a,side='right'); hj=np.searchsorted(E,b,side='right')
        assert min(np.min(hi),np.min(hj))>0 and max(np.max(hi),np.max(hj))<len(E)
        v=(a-E[hi-1])/(E[hi]-E[hi-1])
        u=(b-E[hj-1])/(E[hj]-E[hj-1])
        raw=np.c_[I[ix],hi-1,hi,hj-1,hj]
        signed=np.c_[-np.ones_like(I[ix]),1-v,v,1-u,u]
        self.active=np.unique(raw)
        jj=np.searchsorted(self.active,raw)
        shp=(len(a),len(self.active))
        rr=np.repeat(np.arange(len(a)),8)
        self.S=coo_matrix((signed.ravel(),(rr,jj.ravel())),shape=shp).tocsr()
        self.A=coo_matrix((abs(signed).ravel(),(rr,jj.ravel())),shape=shp).tocsr()
        self.E=E[self.active]
        self.vol=vol[self.active]
        self.cold=np.isin(self.active,cold)
        self.W=w0*amp_scale*LAMBDA**4/m**8*threshold()['event_coefficient_times_lambda4_over_m8']
        self.incoming=jj[:,:4]
        self.outgoing=jj[:,4:]
        self.outw=signed[:,4:]
        assert np.max(abs(self.S@self.E))<2e-14*m
        assert np.max(abs(self.S@np.ones(len(self.active))+2))<5e-15


if __name__=='__main__':
    from early_pole_scalar_angular_chemical_coupled_reference import AngularChemical
    from early_pole_live_parent_energy_reference import state_moments
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spec')
    a=p.parse_args()
    f=Path(a.spec)
    raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
    s=json.loads(raw)
    rows=[]
    for hot in (0,36,100,300):
        d=AngularChemical({(67,112,122):(65,90,72),(97,156,176):(95,125,108)}[tuple(s['cells'])],s['mass_ratio'],0.,3)
        rem,meta=d.setup_hot_grid(s,hot)
        n=np.asarray(rem['endpoint_scalar_occupations'])
        x=s['x_end']
        fs=np.array(rem['endpoint_fermion_occupations'])
        ps=np.array(rem['endpoint_parent_occupations'])
        E,D,P,m=state_moments(d,x,np.r_[n,fs,ps])
        H=np.sqrt(np.pi**2*106.75/90)*(86./x)/2.435e18*np.sqrt((d.radiation+D+P)/d.radiation)
        for gauss in (1,2,4):
            obj=IntervalNumberChange(d.d.p,d.d.ed,m,1e-6,1.,3,gauss)
            fl=obj.flux(n[obj.active])
            e=abs(fl@obj.E)/np.sum(abs(fl)*obj.E)
            assert e<1e-12
            rows.append(dict(outgoing_gauss_points_per_interval=gauss,events=obj.S.shape[0],
                             formal_net_loss_per_log_a=float(-fl.sum()/(d.d.w@n)/H),
                             energy_conservation_relative=e,**meta,static=obj.verify(m)))
    print(json.dumps(rows,indent=2,allow_nan=False))
