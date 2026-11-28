"""Cold angular 4-to-2 comparison with an optional locally refined hot mesh."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from early_pole_live_parent_energy_reference import state_moments
from early_pole_scalar_chemical_sensitivity_reference import Chemical
from early_pole_scalar_number_change_angular_reference import AngularNumberChange
from early_pole_scalar_number_change_interval_reference import IntervalNumberChange


class AngularChemical(Chemical):
    def __init__(self,shape,ratio,strength,power,split=0):
        super().__init__(*shape,ratio,strength)
        self.power=power
        self.split=split

    def collider(self,x,m):
        h=np.sqrt(np.pi**2*106.75/90)*(86./x)/2.435e18
        cl=(IntervalNumberChange(self.d.p,self.d.ed,m,self.qlo,self.strength,
                                 power=self.power,angular_order=self.split)
            if self.split else AngularNumberChange(self.d.p,self.d.ed,m,self.qlo,self.strength,power=self.power))
        return cl,h

    def setup_hot_grid(self,s,hot_cells,bounds=(.0006,.0025)):
        ed=self.d.ed
        assert np.array_equal(ed,s['scalar_momentum_edges'])
        old=np.asarray(s['endpoint_scalar_occupations'])
        x0=s['x_end']
        f=s['endpoint_fermion_occupations'];p=s['endpoint_parent_occupations']
        _,d0,_,m0=state_moments(self,x0,np.r_[old,f,p])
        n0=self.d.w@old
        if hot_cells:
            new=np.unique(np.r_[ed,np.linspace(*bounds,hot_cells+1)])
            points=np.sqrt(new[:-1]*new[1:])
            i=np.searchsorted(ed,points,side='right')-1
            mapped=old[i]
            self.d.ed=new
            self.d.p=points
            self.d.w=np.diff(new**3)/(6*np.pi**2)
            self.d.n=len(points)
            self.ns=self.d.n
            self.iP=self.ns+self.nf
            self.n=self.iP+self.np
            self.w=np.r_[self.d.w,self.d.wf,self.pw]
            s=dict(s,scalar_momentum_edges=new.tolist(),
                   scalar_momentum_nodes=points.tolist(),endpoint_scalar_occupations=mapped.tolist())
        else:
            mapped=old
        _,d1,_,m1=state_moments(self,x0,np.r_[mapped,f,p])
        er=d1/d0-1
        gr=m1/m0-1
        assert abs(self.d.w@mapped/n0-1)<1e-12
        assert abs(er)<8e-7 and abs(gr)<8e-7
        return s,dict(extra_hot_partition_cells=hot_cells,scalar_cells_after_remesh=self.ns,
                      energy_relative_remesh_defect=er,pole_relative_remesh_defect=gr)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectrum')
    p.add_argument('--strength',type=float,default=1e-7)
    p.add_argument('--power',type=int,default=2)
    p.add_argument('--hot-cells',type=int,default=0)
    p.add_argument('--split-outgoing',dest='split',type=int,default=0)
    p.add_argument('--xend',type=float,default=.052)
    a=p.parse_args()
    path=Path(a.spectrum)
    raw=gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()
    s=json.loads(raw)
    shape={(67,112,122):(65,90,72),(97,156,176):(95,125,108)}
    d=AngularChemical(shape[tuple(s['cells'])],s['mass_ratio'],a.strength,a.power,a.split)
    s,meta=d.setup_hot_grid(s,a.hot_cells)
    z=d.from_snapshot(s,a.xend)
    z.update(meta,angular_samples=2**a.power,outgoing_interval_order=a.split)
    print(json.dumps(z,indent=2,allow_nan=False))
