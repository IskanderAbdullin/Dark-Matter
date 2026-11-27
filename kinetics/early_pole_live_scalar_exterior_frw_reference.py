"""Short joined exterior evolution with a prescribed log-scale path and held polarization."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
from types import MethodType

import numpy as np
from scipy import __version__ as scipy_version
from scipy.integrate import solve_ivp

from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_cold_screened_split_reference import atomic_write,state_digest
from early_pole_archived_live_compton_reference import ComptonKernel
from early_pole_archived_live_flavor_conversion_reference import ConversionKernel
from early_pole_joined_scalar_exterior_source_reference import band_kernel
from early_pole_joined_screened_frw_reference import exterior,finite_energy_stable,elastic_at_midpoint
from early_pole_live_joint_hard_packets_frw_reference import Expanding
from early_pole_live_scalar_exterior_decay_reference import warm_shape
from early_pole_native_aperture_witness_reference import GRID
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels


def run(bundle,projection,duration,panels,expands,warm_order=None,checkpoint=None,resume=False,
        checkpoint_every=6,stop_after=None):
    if checkpoint_every<1 or (stop_after is not None and (not checkpoint or not 1<=stop_after<=2*panels)):
        raise ValueError("Checkpoint controls must name valid stopped macro steps.")
    if resume and not checkpoint:
        raise ValueError("A resumed run needs a checkpoint.")
    old=None
    if checkpoint and Path(checkpoint).exists():
        if not resume:
            raise FileExistsError("Checkpoint exists; resume explicitly or choose a fresh path.")
        old=json.loads(Path(checkpoint).read_text())
    elif resume:
        raise FileNotFoundError(checkpoint)
    physics=None
    if checkpoint:
        code=[]
        for fn in (run,exterior,finite_energy_stable,elastic_at_midpoint,
                   elastic_at_midpoint.__globals__['step'],Expanding.weights,warm_shape,band_kernel,
                   Screened.__init__,Screened.scalar,ComptonKernel.collision,
                   ConversionKernel.collision,ElasticChannels.collision):
            path=Path(fn.__code__.co_filename)
            code.append(hashlib.sha256(path.read_bytes()).hexdigest())
        physics=dict(format=1,projection=projection,duration=duration,panels=panels,expands=expands,
                     warm_order=warm_order,bundle_sha=hashlib.sha256(Path(bundle).read_bytes()).hexdigest(),
                     code=code,numpy=np.__version__,scipy=scipy_version)
    results=[] if old is None else old["finished_arms"]
    if old is not None and (not isinstance(results,list) or len(results)>1 or old.get('arm')!=len(results)):
        raise ValueError("Checkpoint arm and completed ledger disagree.")
    if old is not None and stop_after is not None and stop_after<old["global_steps_done"]:
        raise ValueError("The requested stopping boundary is before the checkpoint.")
    with gzip.open(bundle,'rt') as fp:
        raw=json.load(fp)['branches'][projection]
    grid=raw['grid']
    if checkpoint:
        arch=GRID/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        physics['archive_sha']=hashlib.sha256(arch.read_bytes()).hexdigest()
        if old is not None and old.get("physics")!=physics:
            raise ValueError("Checkpoint does not match the bundle, archive, controls, or solver.")
    shot=raw['snapshots'][-1]
    Y0=np.array(shot['hard_moments'])
    N0=np.array(shot['cold_bins'])
    cold_n=len(N0)
    C=Screened(grid,1,cold_n,raw['duration'],(1,1,2),.08,1.,raw['aperture_grid'],projection)
    M=C.M
    for x,y in zip(M.edge_sets,raw['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)
    M.weights=MethodType(Expanding.weights,M)
    M.energy_increment=MethodType(finite_energy_stable,M)
    M.spectators=exterior(M,grid)
    M.eta=0.
    assert max(abs(M.spectators(0.)[0]-M.background))<3e-12
    with gzip.open(GRID/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz','rt') as fp:
        archive=json.load(fp)
    e=np.asarray(archive['scalar_momentum_edges'])**3
    selected=(C.spq>=.1)&(C.spq<=5.)&(C.spect>0.)
    cold_q,cold_w=C.q.copy(),C.sw.copy()
    q,w,N=C.spq[selected],np.diff(e)[selected]/(6*np.pi**2),C.spect[selected]
    initial_gap=C.scalar(N0)[-1]
    shape={}
    if warm_order is not None:
        old_rho=C.scalar(N0)[1]
        q,w,N,*_=warm_shape(e[:-1][selected],e[1:][selected],q,w,N,initial_gap,warm_order)
        q,w,N=q.ravel(),w.ravel(),N.ravel()
        tail=C.spect[~selected]
        C.spq=np.r_[C.spq[~selected],q]
        C.spect=np.r_[tail,N]
        selected=np.r_[np.zeros(len(tail),dtype=bool),np.ones(len(q),dtype=bool)]
        new=C.scalar(N0)
        assert abs(new[1]-old_rho)<2e-21
        shape=dict(warm_order=warm_order,marker_pole=initial_gap,reconstructed_pole=new[-1],
                   initial_scalar_energy_shift=float(new[1]-old_rho))
    spectator=C.spect.copy()
    expected_pole=C.scalar(N0)[-1]
    elastic=ElasticChannels(C,grid,1.)
    co,cv=ComptonKernel(M),ConversionKernel(M)
    ef,sf,cof,cvf=C.factor,elastic.soft.fac,co.fac,cv.fac
    nor=M.norm
    dt=duration/panels

    def clock(t):
        eta=t if expands else 0.
        M.eta=eta
        a=np.exp(eta)
        M.norm=nor/a
        C.factor,elastic.soft.fac,co.fac,cv.fac=ef*a,sf*a,cof*a,cvf*a

    def branch(active):
        arm=int(active)
        C.q,C.sw=(np.r_[cold_q,q],np.r_[cold_w,w]) if active else (cold_q,cold_w)
        C.n=len(C.q)
        C.spect=np.where(selected,0.,spectator) if active else spectator.copy()
        Nall=np.r_[N0,N] if active else N0
        g0=C.scalar(Nall)[-1]
        assert abs(g0-expected_pole)<4e-18
        assert abs(Nall.sum()+C.spect.sum()-N0.sum()-spectator.sum())<1e-16
        z0=np.r_[Y0.ravel(),Nall,0.,0.,0.]

        def totals(t,z):
            clock(t)
            y=z[:C.j0].reshape(-1,2)
            _,_,_,rho,s=M.weights(y,True)
            _,sr,ss,_,pole=C.scalar(z[C.j0:-3])
            return rho+sr,s+ss,pole,float(y[:,0].sum())

        r0,s0,p0,c0=totals(0.,z0)
        work_max=0.
        def rhs(t,z):
            nonlocal work_max
            clock(t)
            Y=z[:C.j0].reshape(-1,2)
            E,F,G,*_=M.weights(Y,True)
            trace=M.trace
            se,_,_,f,pole=C.scalar(z[C.j0:-3])
            hard=np.zeros_like(Y)
            ds=np.zeros(C.n)
            entropy=0.
            for i,(p,wi,fi) in enumerate(zip(C.q,C.sw,f,strict=True)):
                v,d,pr=band_kernel(C,E,F,p,pole,fi,8,ledger=True,scale_log=M.eta)
                hard+=wi*d
                ds[i]=wi*v[2]
                entropy+=wi*pr
            dec_work=float(np.sum(G*hard)+se@ds)
            work_max=max(work_max,abs(dec_work))
            assert abs(dec_work)<3e-13 and entropy>=0.
            work=dec_work
            pair=M.pair_rhs(Y,True)
            hard+=pair[:-2].reshape(-1,2)
            entropy+=pair[-2];work+=pair[-1]
            for kernel in (co,cv):
                d,s,w0,*_=kernel.collision(E,F,G)
                hard+=d;entropy+=s;work+=w0
            return np.r_[hard.ravel(),ds,entropy,work,trace if expands else 0.]

        atol=np.r_[np.full(C.j0,3e-15),np.full(cold_n,2e-18),
                   np.full(C.n-cold_n,2e-23),3e-16,3e-20,3e-19]
        def advance(a,b,state):
            sol=solve_ivp(rhs,(a,b),state,method='DOP853',rtol=2e-11,atol=atol)
            assert sol.success,sol.message
            return sol.y[:,-1].copy()
        init_rate=rhs(0.,z0)[C.j0:-3]
        z=z0.copy()
        implicit=[]
        begin=0
        if old is not None and old['arm']==arm:
            z=np.asarray(old['state'])
            begin=old['arm_step']
            implicit=old['implicit']
            work_max=old['maximum_decay_work']
            if (old['initial_digest']!=state_digest(z0) or old['state_digest']!=state_digest(z)
                or z.shape!=z0.shape or not np.all(np.isfinite(z)) or not 0<=begin<=panels
                or old['global_steps_done']!=arm*panels+begin or len(implicit)!=begin):
                raise ValueError("Stored warm FRW state is invalid.")

        def save(k):
            atomic_write(checkpoint,dict(physics=physics,finished_arms=results,arm=arm,arm_step=k,
                global_steps_done=arm*panels+k,initial_digest=state_digest(z0),state_digest=state_digest(z),
                state=z.tolist(),implicit=implicit,maximum_decay_work=work_max))

        for j in range(begin,panels):
            m=(j+.5)*dt
            z=advance(j*dt,m,z)
            clock(m)
            y,d=elastic_at_midpoint(elastic,z[:C.j0].reshape(-1,2),dt)
            z[:C.j0]=y.ravel()
            z[-3]+=d['entropy_gain']
            implicit.append(d)
            z=advance(m,(j+1)*dt,z)
            k=j+1
            end_now=stop_after==arm*panels+k
            if checkpoint and (k%checkpoint_every==0 or k==panels or end_now):
                save(k)
            if end_now:
                return None
        r1,s1,p1,c1=totals(duration,z)
        net=z[C.j0:-3]-Nall
        ledger=r1-r0-z[-1]
        assert abs(ledger)<3e-12 and abs(ledger-z[-2])<3e-12
        assert abs(s1-s0-z[-3])<3e-11 and abs(c1-c0)<3e-12
        return dict(warm_active=active,cold_gain=float(net[:cold_n].sum()),warm_gain=float(net[cold_n:].sum()),
                    initial_count_rates=[float(init_rate[:cold_n].sum()),float(init_rate[cold_n:].sum())],
                    initial_gap=p0,final_gap=p1,hard_count_defect=c1-c0,
                    redshift_trace_integral=z[-1],collision_work=z[-2],energy_change=r1-r0,
                    redshift_subtracted_defect=ledger,entropy_defect=s1-s0-z[-3],
                    maximum_decay_work=work_max,elastic_splits=sum('elastic_substeps' in d for d in implicit))

    for idx in range(len(results),2):
        out=branch(bool(idx))
        if out is None:
            return dict(checkpoint_only=True,steps_done=stop_after,total_steps=2*panels)
        results.append(out)
    a,b=results
    np.testing.assert_allclose(a['initial_count_rates'][0],b['initial_count_rates'][0],rtol=2e-13)
    return dict(scope='Short comoving warm continuation under prescribed radiation clock and held polarization.',
                projection=projection,grid=grid,duration=duration,panels=panels,expands=expands,
                final_log_scale=duration if expands else 0.,initial_warm_nodes=q.tolist(),arms=[a,b],**shape)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    p.add_argument('--duration',type=float,default=2e-6)
    p.add_argument('--panels',type=int,default=6)
    p.add_argument('--expands',action='store_true')
    p.add_argument('--warm-order',type=int)
    p.add_argument('--checkpoint')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--checkpoint-every',type=int,default=6)
    p.add_argument('--stop-after',type=int)
    a=p.parse_args()
    if a.duration<=0. or a.panels<1 or (a.warm_order is not None and a.warm_order<3):
        p.error('Positive duration, positive macro steps and at least three distributed nodes required.')
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
