"""Local evolving hard window with archived fermions and spectator scalar gap."""
import argparse
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import logit
from early_pole_live_pauli_gauge_moment_decay_reference import Tool,Cell
from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_snapshot_gauge_ledger_reference import photon,cell_J


class Local(Tool):
    def __init__(self,grid,refine,order):
        d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
        super().__init__(2,.005,order,True)
        x=d['x_end']
        masses=[x,d['mass_ratio']*x]
        tags=['parent','fermion']
        p=[np.asarray(d[f'endpoint_{l}_momentum_edges']) for l in tags]
        f=[np.asarray(d[f'endpoint_{l}_occupations']) for l in tags]
        lo,hi=1.30,1.56
        native=np.unique(np.r_[lo,*[v[(v>lo)&(v<hi)] for v in p],hi])
        # Coalesce numerical copies of the same geometric edge shared by
        # the two archives; a roundoff-width cell has no well-conditioned energy moment.
        self.merged_edges=int(sum(np.diff(native)<1e-10))
        native=native[np.r_[True,np.diff(native)>=1e-10]]
        edges=np.unique(np.concatenate([np.linspace(l,h,refine+1)
                                        for l,h in zip(native[:-1],native[1:],strict=True)]))
        self.n=len(edges)-1
        self.cells=[]
        init=[]
        totals=np.array([cell_J(a,b,c,.5) for a,b,c in zip(p,f,masses,strict=True)])
        for b,u,m in zip(p,f,masses,strict=True):
            for l,h in zip(edges[:-1],edges[1:],strict=True):
                c=Cell(l,h,m,order)
                self.cells.append(c)
                i=np.searchsorted(b,(l+h)/2)-1
                init.append(c.initial(logit(u[i]),0.))
        self.ec=np.array([v.ec for v in self.cells])
        self.de=np.array([v.de for v in self.cells])
        J0=np.array([c.fit(u)[1] for u,c in zip(init,self.cells,strict=True)])
        self.Js=totals-np.r_[sum(J0[:self.n]),sum(J0[self.n:])]
        assert min(self.Js)>=0.
        self.jgamma=photon(1.,0.,.5)
        sx=np.asarray(d['scalar_momentum_edges'])
        i=np.searchsorted(sx,self.q)-1
        y=np.asarray(d['endpoint_scalar_signed_shell_number_over_T3'])
        assert np.sum(abs(np.minimum(y,0.)))<1e-10*y.sum()
        self.ms=d['checkpoints'][-1]['scalar_pole_over_T']
        self.es=np.hypot(self.ms,self.q)
        p0=np.asarray(d['scalar_momentum_nodes'])
        z=np.arange(len(y))!=i
        self.specq,self.specn=p0[z],np.maximum(y[z],0.)
        self.wscalar=y[i]/d['endpoint_scalar_occupations'][i]
        self.clocknorm=d['endpoint_scalar_occupations'][i]
        self.lambda_s=2*self.ms**2/(y[i]/self.es+np.sum(self.specn/np.hypot(self.specq,self.ms)))
        self.init=np.r_[np.array(init).ravel(),y[i]]

    def scalar(self,n):
        def test(ms):
            r=np.hypot(self.q,ms)
            pp=np.hypot(self.specq,ms)
            return ms*ms-.5*self.lambda_s*(n/r+self.specn@np.reciprocal(pp))
        ms=brentq(test,1e-8,.01,xtol=1e-17)
        es=np.hypot(self.q,ms)
        rho=n*es+self.specn@np.hypot(self.specq,ms)-ms**4/(2*self.lambda_s)
        return es,rho


def run(grid,refine,order):
    M=Local(grid,refine,order)
    return dict(scope='Local hard decay window initialized from archived f; only one scalar bin evolves. '
                'Other scalars contribute a live Hartree field at fixed spectator counts. '
                'Exposure divided by initial scalar occupancy, no physical rate or expansion.',
                grid=grid,refinement=refine,cells_per_species=M.n,coalesced_roundoff_edges=M.merged_edges,
                scalar_start_number=M.init[-1],initial_scalar_occupation=M.clocknorm,
                rows=[M.solve(x,.10) for x in ('bare','closed')])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grid',choices=['fine','third'],default='third')
    p.add_argument('--refine',type=int,default=1)
    p.add_argument('--order',type=int,default=12)
    a=p.parse_args()
    print(json.dumps(run(a.grid,a.refine,a.order),indent=2,allow_nan=False))
