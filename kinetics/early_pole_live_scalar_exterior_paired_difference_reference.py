"""Evolve baseline and the small difference under a shared time grid."""
import argparse
import gzip
import json

import numpy as np

from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_joined_scalar_exterior_source_reference import band_kernel
from early_pole_native_aperture_witness_reference import GRID


def run(path,projection,duration,steps,order,outer,joined=False,warm_multiplier=1.,elastic_energy_budget=None,attribute_feedback=False,warm_order=None,stable_pole_tangent=False,anchor_pole=False,linearize_cold_pole=False,cold_pole_relative_probe=5e-5):
    if linearize_cold_pole and (not anchor_pole or attribute_feedback):
        raise ValueError("The paired cold tangent requires anchoring and its own rate ledger.")
    if not np.isfinite(cold_pole_relative_probe) or not 1e-6<=cold_pole_relative_probe<=1e-3:
        raise ValueError("Use a finite relative pole probe in [1e-6,1e-3].")
    if cold_pole_relative_probe!=5e-5 and not linearize_cold_pole:
        raise ValueError("The relative pole probe requires the matched linear option.")
    if stable_pole_tangent and not attribute_feedback:
        raise ValueError("The anchored pole tangent requires the feedback attribution ledger.")
    if not np.isfinite(warm_multiplier) or warm_multiplier<0.:
        raise ValueError("Warm collision multiplier must be finite and nonnegative.")
    if elastic_energy_budget is not None and (not joined or not np.isfinite(elastic_energy_budget)
                                               or elastic_energy_budget<=0.):
        raise ValueError("An explicit elastic energy budget needs a joined run and a positive finite value.")
    with gzip.open(path,'rt') as f:raw=json.load(f)['branches'][projection]
    g=raw['grid']
    with gzip.open(GRID/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz','rt') as f:d=json.load(f)
    snap=raw['snapshots'][-1]
    y0=np.asarray(snap['hard_moments'])
    n0=np.asarray(snap['cold_bins'])
    if joined:
        from early_pole_archived_cold_screened_elastic_reference import Screened
        from early_pole_archived_live_compton_reference import ComptonKernel
        from early_pole_archived_live_flavor_conversion_reference import ConversionKernel
        from early_pole_screened_active_response_update_reference import ActiveResponse
        from early_pole_screened_elastic_channel_sum_reference import ElasticChannels
        from early_pole_screened_elastic_discrete_gradient_reference import step
        C=Screened(g,1,len(n0),raw['duration'],(1,1,2),.08,1.,raw['aperture_grid'],projection)
        active=ActiveResponse(C,.12)
        soft=ElasticChannels(C,g,1.)
        add=(ComptonKernel(C.M),ConversionKernel(C.M))
    else:
        C=Cold(g,1,0,len(n0),True,exposure=raw['duration'],reference_aperture_grid=raw['aperture_grid'],
               archive_projection=projection)
    for a,b in zip(C.M.edge_sets,raw['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(a,b)
    c_q,c_w=C.q.copy(),C.sw.copy()
    orig=C.spect.copy()
    warm=(C.spq>=outer[0])&(C.spq<=outer[1])&(orig>0)
    q=C.spq[warm]; edges=np.array(d['scalar_momentum_edges'])
    overlap=np.maximum(0.,np.minimum(edges[1:]**3,1e-15)-np.maximum(edges[:-1]**3,1e-18))
    vol=((np.diff(edges**3)-overlap)/(6*np.pi**2))[warm]
    assert len(vol)>0 and min(vol)>0.
    w0=orig[warm]
    shape_ledger={}
    if warm_order is not None:
        if warm_order<3 or np.any(overlap[warm]>0.):
            raise ValueError("Distributed warm shells require at least three nodes and full archive shells.")
        from early_pole_live_scalar_exterior_decay_reference import warm_shape
        original=C.scalar(n0)
        lo,hi=edges[:-1][warm]**3,edges[1:][warm]**3
        q_new,w_new,n_new,*_=warm_shape(lo,hi,q,vol,w0,original[-1],warm_order)
        q,vol,w0=q_new.ravel(),w_new.ravel(),n_new.ravel()
        tails=orig[~warm]
        C.spq=np.r_[C.spq[~warm],q]
        orig=np.r_[tails,w0]
        warm=np.r_[np.zeros(len(tails),dtype=bool),np.ones(len(w0),dtype=bool)]
        C.spect=orig
        new=C.scalar(n0)
        assert abs(new[1]-original[1])<2e-21
        shape_ledger=dict(warm_quadrature_order=warm_order,original_marker_pole=original[-1],
                          common_initial_pole=new[-1],initial_energy_change=new[1]-original[1])
    specs=(orig,np.where(warm,0.,orig))
    banks=((c_q,c_w),(np.r_[c_q,q],np.r_[c_w,vol]))
    h0,nw=len(n0)+C.j0,len(q)
    initial=np.r_[y0.ravel(),n0]
    warm_slice=slice(2*h0,2*h0+nw)
    extra_start=warm_slice.stop
    n_attribution=6 if stable_pole_tangent else 3
    U=np.zeros(extra_start+n_attribution*len(n0) if attribute_feedback else extra_start)

    def set_bank(arm):
        C.spect=specs[arm]
        C.q,C.sw=banks[arm]
        C.n=len(C.q)

    anchor_audit=[]
    closure_audit=[]
    tangent_audit=[]

    def split(z,arm,anchored=None,linear=None):
        set_bank(arm)
        y=z[:C.j0].reshape(-1,2)
        N=z[C.j0:]
        E,F,G,*_=C.M.weights(y,True)
        En,_,_,f,m=C.scalar(N)
        if anchored is not None:
            gap_error=abs(m-anchored)
            assert gap_error<2e-17
            anchor_audit.append(float(gap_error))
            m=anchored
            En=np.hypot(C.q,m)
        rate=np.zeros((2,h0))
        scalar=np.zeros(C.n)
        correction=np.zeros(h0)
        for i,w in enumerate(C.sw):
            exterior=i>=len(n0)
            local_mass=m if linear is None or exterior else linear[0]
            v,H,_=band_kernel(C,E,F,C.q[i],local_mass,f[i],radial_order=order,ledger=True)
            if linear is not None and not exterior:
                mp,dm=linear
                epsilon=mp*cold_pole_relative_probe
                def pair(mu):
                    out,hard,_=band_kernel(C,E,F,C.q[i],mu,f[i],radial_order=order,ledger=True)
                    return out[2],hard
                sa,ha=pair(mp-epsilon)
                sb,hb=pair(mp+epsilon)
                fac=w*dm/(2*epsilon)
                correction[:C.j0]+=fac*(hb-ha).ravel()
                correction[C.j0+i]=fac*(sb-sa)
            effective_w=w*warm_multiplier if exterior else w
            rate[int(exterior),:C.j0]+=effective_w*H.ravel()
            scalar[i]=effective_w*v[2]
        rate[0,C.j0:]=scalar[:len(n0)]
        work_basic=En@scalar+np.sum(G*rate[:,:C.j0].sum(axis=0).reshape(y.shape))
        work_fix=En[:len(n0)]@correction[C.j0:]+np.sum(G*correction[:C.j0].reshape(y.shape))
        defect=work_basic+work_fix
        if linear is not None:
            closure_audit.append(float(abs(defect)))
            tangent_audit.append(float(correction[C.j0:].sum()))
        assert abs(defect)<3e-13
        return rate,scalar[len(n0):],m,correction

    def other_slow(Y):
        out=[C.M.pair_rhs(Y,True)[:-2]]
        state=C.M.weights(Y,True)[:3]
        for source in add:
            d,*_=source.collision(*state)
            out.append(d.ravel())
        return np.asarray(out)

    def cold_at(Y,N,pole):
        E,F,*_=C.M.weights(Y.reshape(-1,2),True)
        return np.array([w*band_kernel(C,E,F,q,pole,float(x/w),radial_order=order)[2]
                         for q,w,x in zip(c_q,c_w,N,strict=True)])

    def increments(base,diff,warm_delta):
        N=np.r_[base[C.j0:],specs[0]]
        dN=np.zeros_like(N)
        dN[:len(n0)]=diff[C.j0:]
        dN[len(n0)+np.flatnonzero(warm)]=warm_delta
        return np.r_[c_q,C.spq],N,dN

    def gap_delta(base,diff,warm_delta,m0):
        # Exact difference of the two common gap equations; rationalize
        # the difference in reciprocal one-particle energies.
        nodes,N,dN=increments(base,diff,warm_delta)
        A=np.hypot(nodes,m0)
        m=m0
        dm=0.
        for _ in range(8):
            B=np.hypot(nodes,m)
            new=(.5*C.lam*(dN@np.reciprocal(B)))/(
                (m+m0)*(1.+.5*C.lam*(N@np.reciprocal(A*B*(A+B)))))
            if abs(new-dm)<1e-29:
                return float(new)
            dm=new
            m=m0+dm
        raise RuntimeError("The anchored scalar gap correction did not converge.")

    def attribute(base,diff,first,second,m0,m1,warm_delta):
        hard_base=base[:C.j0]
        hard_add=hard_base+diff[:C.j0]
        nb=base[C.j0:]
        # Order is hard state, then scalar pole, then cold occupation.
        at_hard=cold_at(hard_add,nb,m0)
        at_pole=cold_at(hard_add,nb,m1)
        d0=at_hard-first[0,C.j0:]
        d1=at_pole-at_hard
        d2=second[0,C.j0:]-at_pole
        if not stable_pole_tangent:
            return np.r_[d0,d1,d2]
        # At fixed q, d(m² - λJ/2)/dm is analytically positive.
        # Use the evolution variables for δN, before floating subtraction
        # of a short accumulated gain from its larger initial counts.
        nodes,N,dN=increments(base,diff,warm_delta)
        em=np.hypot(nodes,m0)
        dm=(.5*C.lam*(dN@np.reciprocal(em)))/(2*m0+.5*C.lam*m0*(N@em**-3))
        t=[]
        for epsilon in (m0*1e-4,m0*5e-5):
            slope=(cold_at(hard_add,nb,m0+epsilon)-cold_at(hard_add,nb,m0-epsilon))/(2*epsilon)
            t.append(dm*slope)
        # Replacing only the pole piece in the telescoping sum is a
        # perturbative sensitivity, not an alternative evolved trajectory.
        return np.r_[d0,d1,d2,t[0],t[1],d0+t[1]+d2]

    def evaluate(U):
        base=initial+U[:h0]
        diff=U[h0:2*h0]
        first,_,m0,_=split(base,0)
        full=np.r_[base+diff,w0+U[warm_slice]]
        dm=gap_delta(base,diff,U[warm_slice],m0) if anchor_pole else 0.
        if linearize_cold_pole and abs(dm)>1e-4*m0*cold_pole_relative_probe:
            raise ValueError("The requested cold pole displacement is outside the near-null control.")
        anchored=m0+dm if anchor_pole else None
        second,w,m1,correction=split(full,1,anchored,(m0,dm) if linearize_cold_pole else None)
        delta=(second[0]-first[0])+second[1]
        if linearize_cold_pole:
            delta+=correction
        if joined:
            a0=other_slow(base[:C.j0].reshape(-1,2))
            a1=other_slow((base[:C.j0]+diff[:C.j0]).reshape(-1,2))
            first[0,:C.j0]+=a0.sum(axis=0)
            delta[:C.j0]+=(a1-a0).sum(axis=0)
        if attribute_feedback:
            return np.r_[first[0],delta,w,attribute(base,diff,first,second,m0,m1,U[warm_slice])]
        return np.r_[first[0],delta,w]

    def total(state,arm):
        set_bank(arm)
        y=state[:C.j0].reshape(-1,2)
        _,_,_,en,s=C.M.weights(y,True)
        _,r,b,*_=C.scalar(state[C.j0:])
        return np.array([en+r,s+b])

    e0=total(initial,0)
    e1=total(np.r_[initial,w0],1)
    assert abs(e0[0]-e1[0])<3e-16
    h=duration/steps
    rows=[]
    elastic_drifts=[[],[]]
    difference_roundoff=[]
    def rk(v,dt):
        A=evaluate(v)
        B=evaluate(v+dt*A/2)
        Cc=evaluate(v+dt*B/2)
        D=evaluate(v+dt*Cc)
        return dt*(A+2*B+2*Cc+D)/6
    for i in range(steps):
        if joined:
            U+=rk(U,h/2)
            hb=initial[:C.j0]+U[:C.j0]
            hd=U[h0:h0+C.j0].copy()
            ha=hb+hd
            moved=[]
            for arm,v in enumerate((hb,ha)):
                yy=v.reshape(-1,2)
                active.refresh(yy,'full')
                if elastic_energy_budget is None:
                    end,diag=step(soft,yy,h)
                else:
                    end,diag=step(soft,yy,h,energy_budget=elastic_energy_budget)
                elastic_drifts[arm].append(diag['rho_drift'])
                moved.append(end.ravel()-v)
            U[:C.j0]+=moved[0]
            U[h0:h0+C.j0]+=moved[1]-moved[0]
            hb_new=initial[:C.j0]+U[:C.j0]
            residual=(hb_new+U[h0:h0+C.j0])-(ha+moved[1])
            difference_roundoff.append(float(np.max(abs(residual))))
            U+=rk(U,h/2)
        else:
            U+=rk(U,h)
        if (i+1)%max(1,steps//4)==0 or i+1==steps:
            line=dict(exposure=h*(i+1),cold_gain=float(U[C.j0:h0].sum()),
                      outer_feedback_on_cold=float(U[h0+C.j0:2*h0].sum()),
                      warm_gain=float(U[warm_slice].sum()))
            if attribute_feedback:
                contributions=U[extra_start:].reshape(n_attribution,len(n0))
                total_delta=U[h0+C.j0:2*h0]
                defect=float(np.max(abs(contributions[:3].sum(axis=0)-total_delta)))
                assert defect < 1e-28
                line.update(cold_feedback_hard_then_pole_then_population=contributions[:3].sum(axis=1).tolist(),
                            maximum_attribution_sum_defect=defect)
                if stable_pole_tangent:
                    line.update(stable_linear_pole_pieces=contributions[3:5].sum(axis=1).tolist(),
                                estimate_replacing_direct_pole_piece=float(contributions[5].sum()))
            rows.append(line)
    base=initial+U[:h0]
    t0=total(base,0)-e0
    t1=total(np.r_[base+U[h0:2*h0],w0+U[warm_slice]],1)-e1
    extra={}
    if joined:
        hard_count=np.sum(y0[:,0])
        count=[float(np.sum(x.reshape(-1,2)[:,0])-hard_count)
               for x in (base[:C.j0],base[:C.j0]+U[h0:h0+C.j0])]
        assert max(abs(t0[0]),abs(t1[0]),*map(abs,count))<3e-12
        assert min(t0[1],t1[1])>=0. and max(difference_roundoff)<2e-16
        if elastic_energy_budget is not None:
            peaks=[max(map(abs,seq)) for seq in elastic_drifts]
            assert max(peaks)<elastic_energy_budget
            extra.update(elastic_energy_budget=elastic_energy_budget,peak_step_energy_defects=peaks)
        extra.update(joined_macro_steps=steps,slow_rk4_half_steps=2*steps,total_hard_count_defect=count,
                   summed_elastic_energy_defects=[float(sum(x)) for x in elastic_drifts],
                   maximum_elastic_hard_pair_reconstruction_residual=max(difference_roundoff))
    else:
        assert max(abs(t0[0]),abs(t1[0]))<2e-13 and min(t0[1],t1[1])>=0.
    return dict(grid=g,projection=projection,duration=duration,common_rk4_steps=steps,radial_order=order,
                energy_closure=[t0[0],t1[0]],entropy_change=[t0[1],t1[1]],rows=rows,
                **(dict(warm_node_cuts=outer) if outer!=[.1,5.] else {}),
                **(dict(warm_collision_multiplier=warm_multiplier,
                        manufactured_rate_scaling=True) if warm_multiplier!=1. else {}),
                **(dict(uses_anchored_live_gap=True,maximum_pole_difference_from_independent_root=max(anchor_audit))
                   if anchor_pole else {}),
                **(dict(cold_pole_linearized_at_common_state=True,
                        cold_pole_relative_probe=cold_pole_relative_probe,
                        maximum_active_decay_work=max(closure_audit),
                        maximum_cold_rate_tangent=max(map(abs,tangent_audit)))
                   if linearize_cold_pole else {}),
                **shape_ledger,**extra,scope=('Short fixed-scale all-channel matched split on one grid; explicit difference variable.'
                               if joined else 'Same cold hard box and decay-only physics; explicit difference variable on a shared RK4 grid.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('path')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    p.add_argument('--duration',type=float,default=7.2e-5)
    p.add_argument('--steps',type=int,default=24)
    p.add_argument('--order',type=int,default=10)
    p.add_argument('--outer',type=float,nargs=2,default=[.1,5.])
    p.add_argument('--joined',action='store_true')
    p.add_argument('--warm-multiplier',type=float,default=1.,
                   help='Manufactured scaling of reciprocal promoted decay packets, not their occupation.')
    p.add_argument('--elastic-energy-budget',type=float)
    p.add_argument('--attribute-feedback',action='store_true',
                   help='Integrate a telescoping cold response in hard, scalar pole, then cold populations.')
    p.add_argument('--warm-order',type=int)
    p.add_argument('--stable-pole-tangent',action='store_true',
                   help='Independent first-order near-null pole diagnostic of the feedback ledger.')
    p.add_argument('--anchor-pole',action='store_true',
                   help='Solve the active pole from the common baseline gap and the evolved count difference.')
    p.add_argument('--linearize-cold-pole',action='store_true',
                   help='Replace the near-null cold pole response by matched reciprocal tangents.')
    p.add_argument('--cold-pole-relative-probe',type=float,default=5e-5)
    a=p.parse_args()
    if a.duration<=0 or a.steps<1 or a.order<2 or not 1e-5<a.outer[0]<a.outer[1]<=5.:
        p.error('Invalid positive exposure, steps or order.')
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
