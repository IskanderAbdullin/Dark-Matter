"""Compare finite live hard supports at a shared pair reference clock."""
import argparse
import json
import numpy as np
from scipy.special import roots_legendre
from scipy.optimize import brentq

from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_snapshot_gauge_ledger_reference import cell_J,photon
from early_pole_live_joint_hard_pair_packets_reference import Packets
from early_pole_decay_mode_clock_reference import G,KAPPA,intervals


class ArchivePackets(Packets):
    vacuum_clock=False

    def rhs_joint(self,y,closed,decay_scale):
        if self.vacuum_clock and decay_scale:
            E=self.scalar(y[-5])[0]
            ms2=E*E-self.scalarq**2
            decay_scale*=self.Ws*G*KAPPA*.05*((.05+.003)**2-ms2)/(self.scalarq*E)
        return super().rhs_joint(y,closed,decay_scale)

    def scalar(self,n):
        q=self.scalarq
        def gap(m):
            return m*m-.5*self.lam*(n/np.hypot(q,m)+self.specn@np.reciprocal(np.hypot(self.specq,m)))
        m=brentq(gap,1e-8,.02,xtol=1e-17)
        E=np.hypot(q,m)
        rho=n*E+self.specn@np.hypot(self.specq,m)-m**4/(2*self.lam)
        f=n/self.Ws
        S=self.Ws*(np.log1p(f)+f*np.log1p(1/f))
        return E,rho,S,f


def fit_archive_cell(C,edges,f):
    q,sw=roots_legendre(12)
    lo,hi=C.lo,C.hi
    p=np.unique(np.r_[lo,edges[(edges>lo)&(edges<hi)],hi])
    dh=np.diff(p)[:,None]/2
    z=((p[1:]+p[:-1])[:,None]/2+dh*q).ravel()
    w=(dh*sw).ravel()*z*z*2/np.pi**2
    u=f[np.searchsorted(edges,z)-1]
    Y=np.r_[w@u,(w*u)@C.basis(z)]
    return Y


def mount(grid,wide,level,aperture_refined=False):
    d=load(f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz')
    edges=([.5,.92,1.12,1.55,2.0,2.55,3.15,4.],
           [.45,.65,.85,1.05,1.25,1.7,2.2,2.85,3.7,5.],
           [.45,.69,.8,1.03,1.3,1.7,2.2,2.85,3.7,5.]) if wide else None
    if wide==2:
        edges=tuple(z+extra for z,extra in zip(
            edges,([4.7,5.7,6.8,8.],[6.1,7.25,8.55,10.],[6.1,7.25,8.55,10.]),strict=True))
    if aperture_refined:
        if edges is None:
            edges=([.5,.92,1.12,1.55,2.0],
                   [.45,.65,.85,1.05,1.25,1.7,2.2],
                   [.45,.69,.8,1.03,1.3,1.7,2.2])
        add=[1.32,1.56]
        for t in ('parent','fermion'):
            old=np.array(d[f'endpoint_{t}_momentum_edges'])
            add.extend(old[(old>1.32)&(old<1.56)])
        def join(existing):
            a=np.unique(np.r_[existing,add])
            return a[np.r_[True,np.diff(a)>1e-10]].tolist()
        edges=(edges[0],join(edges[1]),join(edges[2]))
        assert len(edges[1])==len(edges[2])
    quad=[(2,3,10,-3),(3,5,14,-3)][level]
    M=ArchivePackets(12,quad,0.,edges=edges)
    total=np.zeros(3)
    active=np.zeros(3)
    original=np.zeros(3)
    # Matching hard background outside active cells is static; reciprocal
    # field energy includes its cross term, never fictitious particle transfers.
    total[0]=photon(1.,0.,.5)
    original[0]=total[0]-photon(1.,0.,M.edge_sets[0][-1])
    for i in range(M.r[0].stop):
        G=M.cells[i]
        M.init[i]=G.moments(1/np.expm1(G.p))
        active[0]+=G.fit(M.init[i])[2]
    for k,name,m in ((1,'parent',.05),(2,'fermion',.003)):
        e=np.asarray(d[f'endpoint_{name}_momentum_edges'])
        f=np.asarray(d[f'endpoint_{name}_occupations'])
        total[k]=cell_J(e,f,m,.45)
        original[k]=total[k]-cell_J(e,f,m,M.edge_sets[k][-1])
        for i in range(M.r[k].start,M.r[k].stop):
            C=M.cells[i]
            M.init[i]=fit_archive_cell(C,e,f)
            active[k]+=C.fit(M.init[i])[1]
    M.background=total-original
    assert np.min(M.background)>=0.
    assert max(abs((M.background+active)/total-1))<.004
    q=M.scalarq
    k=np.searchsorted(d['scalar_momentum_edges'],q)-1
    M.Ns=d['endpoint_scalar_signed_shell_number_over_T3'][k]
    occ=d['endpoint_scalar_occupations'][k]
    M.Ws=M.Ns/occ
    ms=d['checkpoints'][-1]['scalar_pole_over_T']
    spec=np.asarray(d['endpoint_scalar_signed_shell_number_over_T3'])
    assert -np.minimum(spec,0.).sum()<1e-10*spec.sum()
    keep=np.arange(len(spec))!=k
    M.specn=np.maximum(spec[keep],0.)
    M.specq=np.asarray(d['scalar_momentum_nodes'])[keep]
    M.lam=2*ms*ms/(M.Ns/np.hypot(q,ms)+M.specn@np.reciprocal(np.hypot(M.specq,ms)))
    x=.05
    H=np.sqrt(np.pi**2*106.75/90)*(86./x)/2.435e18
    # Spin γγ channel has the reference g⁴/H prefactor; the scalar
    # test multiplier remains imposed. A complete time calibration
    # requires matched physical treatment of both channels and coverage.
    M.norm=H/.005**4
    return M,d,quad,((M.background+active)/total-1).tolist()


def one(grid,wide,level,vacuum_clock=False):
    M,d,quad,relative=mount(grid,wide,level)
    M.vacuum_clock=vacuum_clock
    # Check reduction of the retained scalar coefficient to the existing
    # independent vacuum-width formula in the hard Lorentz-invariant limit.
    E=M.scalar(M.Ns)[0]
    q=M.scalarq
    m=np.sqrt(E*E-q*q)
    _,_,_,_,Ps=intervals(.05,m,.003,q)
    b=G*.05**2/(2*E*q*Ps)*KAPPA*.05**2*(2*Ps/.05)*((.05+.003)**2-m*m)/.05**2
    a=G*KAPPA*.05*((.05+.003)**2-m*m)/(E*q)
    assert abs(a/b-1)<2e-15
    clock=1. if vacuum_clock else 50.*M.Ws/M.Ns
    rows=[]
    for p,q in ((0.,clock),(1.,clock),(1.,0.)):
        M.pair_scale=p
        rows.append(M.go_joint(True,q,duration=8e-5))
    return dict(scope='Archived charged hard distributions, Tref BE photons initially, '
                'fixed-scale live hard and one scalar bin, spectator Hartree moments, reference vacuum or imposed decay clock; not physical evolution.',
                grid=grid,wide=wide,quad=quad,live_edges=M.edge_sets,vacuum_decay_reference_clock=vacuum_clock,
                frozen_background=M.background.tolist(),
                initial_field_moment_reconstruction_relative=relative,rows=rows,
                remap_rejection=M.reject)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grid',choices=['fine','third'],default='third')
    p.add_argument('--level',type=int,default=0)
    p.add_argument('--wide',action='store_true')
    p.add_argument('--wider',action='store_true')
    p.add_argument('--vacuum-clock',action='store_true')
    a=p.parse_args()
    print(json.dumps(one(a.grid,2 if a.wider else int(a.wide),a.level,a.vacuum_clock),indent=2,allow_nan=False))
