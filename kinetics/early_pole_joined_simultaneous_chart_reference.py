"""Unsplit finite early joined reference in admissible charged-logit coordinates."""
import argparse
import gzip
import hashlib
import sys
import scipy
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_archived_shared_cold_window_reference import Cold
from early_pole_archived_live_compton_reference import ComptonKernel
from early_pole_archived_live_flavor_conversion_reference import ConversionKernel
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels


def run(bundle,projection,duration,panels,rtol,base_atol,prepared_input=None,retune_local_daughter=False,compensation_only=False,number_only_local_retune=False,without_elastic=False,elastic_channel='both',record_mode_history=False,local_from_paired_archive=False,retune_donor_support=None,local_daughter_fraction=1.,checkpoint_input=None,emit_checkpoint=False,adaptive=None):
    if not 0.<local_daughter_fraction<=1. or (local_daughter_fraction!=1. and not local_from_paired_archive):
        raise ValueError('A fractional replacement in (0,1] requires the exact local archive option.')
    with gzip.open(bundle,'rt') as f:
        branches=json.load(f)['branches']
    saved=branches[projection]
    shot=saved['snapshots'][0]
    y0=np.asarray(shot['hard_moments'])
    n0=np.asarray(shot['cold_bins'])
    depth=saved.get('aperture_subdivision',1)
    if depth<1:
        raise ValueError('Invalid archived aperture subdivision.')
    if depth!=1 and (prepared_input is not None or retune_local_daughter
                     or compensation_only or number_only_local_retune or retune_donor_support is not None):
        raise ValueError('Prepared control requires its validated original aperture support.')
    C=Screened(saved['grid'],1,len(n0),duration,(1,1,2),.08,1.,saved['aperture_grid'],projection,
               subdivision=depth)
    M=C.M
    np.testing.assert_array_equal(y0,M.init)
    np.testing.assert_array_equal(n0,C.ns)
    for x,y in zip(M.edge_sets,saved['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(x,y)
    swap=None
    if prepared_input is not None:
        src=json.loads(Path(prepared_input).read_text())
        swap=src.get('prepared_charged_swap')
        if not src.get('manufactured_hard_shape_control') or swap not in ('daughter','local-daughter','outside-daughter'):
            raise ValueError("Only accepted whole, aperture-local or outside daughter inputs are supported.")
        if any(src[k]!=saved[k] for k in ('projection','grid','hard_edge_sets','aperture_grid')):
            raise ValueError("Prepared input uses a different hard support.")
        if not abs(src['energy_drift'])<3e-12 or not abs(src['entropy_closure'])<3e-11:
            raise ValueError("Prepared input must come from an accepted solver history.")
        ss=src['snapshots'][0]
        np.testing.assert_array_equal(ss['cold_bins'],n0)
        y0=np.asarray(ss['hard_moments'])
        expected=Cold.rhs(C,np.r_[y0.ravel(),n0,0.,0.],False)[C.j0:-2]
        np.testing.assert_allclose(expected,ss['instantaneous_source_by_q'],rtol=1e-12,atol=2e-14)
    local_archive_meta=None
    if local_from_paired_archive:
        if projection!='steps' or depth not in (1,2,4,8) or any((prepared_input is not None,
                retune_local_daughter, compensation_only, number_only_local_retune, retune_donor_support is not None)):
            raise ValueError('Exact local archive swap requires validated native hard supports.')
        counterpart=branches['cell-conservative']
        if counterpart.get('aperture_subdivision',1)!=depth:
            raise ValueError('Different native and continuous archive depth.')
        other=counterpart['snapshots'][0]
        np.testing.assert_array_equal(n0,other['cold_bins'])
        for e,f in zip(M.edge_sets,counterpart['hard_edge_sets'],strict=True):
            np.testing.assert_array_equal(e,f)
        ed=np.asarray(M.edge_sets[2])
        lo,hi=1.4138888888888888,1.4666666666666666
        loc=np.flatnonzero((ed[:-1]>=lo-2e-14)&(ed[1:]<=hi+2e-14))
        if len(loc)!=2*depth or np.any(np.diff(loc)!=1):
            raise ValueError('Expected nested cells over the original local physical support.')
        np.testing.assert_allclose([ed[loc[0]],ed[loc[-1]+1]],[lo,hi],rtol=0.,atol=2e-14)
        y0=y0.copy()
        ix=M.r[2].start+loc
        from_archive=np.asarray(other['hard_moments'])[ix]
        if local_daughter_fraction==1.:
            y0[ix]=from_archive
        else:
            y0[ix]+=local_daughter_fraction*(from_archive-y0[ix])
        if depth!=1 or local_daughter_fraction!=1.:
            local_archive_meta=dict(indices=ix.tolist(),
                momentum_interval=[float(ed[loc[0]]),float(ed[loc[-1]+1])],
                initialized_cold_rate=float(Cold.rhs(C,np.r_[y0.ravel(),n0,0.,0.],False)[C.j0:-2].sum()))
        swap='local-daughter'
        if local_daughter_fraction!=1.:
            local_archive_meta['replacement_fraction']=local_daughter_fraction
            swap='local-daughter-fractional'
    retune_check=None
    if retune_local_daughter or compensation_only or number_only_local_retune or retune_donor_support is not None:
        if projection!='steps' or prepared_input is not None:
            raise ValueError('Retuned local preparation requires the native graph and no alternate input.')
        from early_pole_daughter_complement_retune_reference import prepare
        independent=branches['cell-conservative']
        for x,y in zip(M.edge_sets,independent['hard_edge_sets'],strict=True):
            np.testing.assert_array_equal(x,y)
        other=independent['snapshots'][0]
        np.testing.assert_array_equal(n0,other['cold_bins'])
        ed=np.asarray(M.edge_sets[2])
        local=np.flatnonzero((ed[:-1]>=1.4138888888888888-2e-14)&(ed[1:]<=1.4666666666666666+2e-14))
        if len(local)!=2:
            raise ValueError('Exactly two initial aperture cells required.')
        fixed=M.r[2].start+local
        native=y0
        match=(True,not number_only_local_retune)
        donor=None
        if retune_donor_support is not None:
            if retune_donor_support not in ('below','no_neighbors') or any((
                    retune_local_daughter, compensation_only, number_only_local_retune)):
                raise ValueError('Require one accepted donor-support control.')
            from early_pole_daughter_equal_moments_donor_support import choose_donors
            span=M.r[2]
            donor=choose_donors(np.arange(span.start,span.stop),fixed)[retune_donor_support]
        _,y0,retune_check=prepare(M,native,np.asarray(other['hard_moments']),fixed,match,donor_ids=donor)
        swap='number-only-retuned-local-daughter' if number_only_local_retune else 'retuned-local-daughter'
        if retune_donor_support is not None:
            swap+='-'+retune_donor_support.replace('_','-')
            retune_check['donor_support_variant']=retune_donor_support
        if compensation_only:
            y0[fixed]=native[fixed]
            swap='local-complement-only'
            retune_check['actual_complement_only_number_and_bare_energy_difference']=[
                float(np.sum((y0-native)[M.r[2],0])),
                float(np.sum(M.shape[M.r[2]]*(y0-native)[M.r[2]]))]

    if elastic_channel not in ('both','unlike','same') or (without_elastic and elastic_channel!='both'):
        raise ValueError('Use one unambiguous elastic channel control.')
    if without_elastic:
        elastic=None
    elif elastic_channel=='unlike':
        elastic=C
    else:
        if elastic_channel=='same':
            C.strength=0.
        elastic=ElasticChannels(C,saved['grid'],1.)
    add=(ComptonKernel(M),ConversionKernel(M))
    k=M.r[1].start
    charts=M.cells[k:]
    fits=M.weights(y0,True)[1][k:]
    anchor=np.r_[y0[:k].ravel(),fits.ravel(),n0,0.,0.]
    cond=[]
    coordinate_defect=[]

    def decode(flat,derivatives):
        st=flat[:C.j0].reshape(y0.shape)
        out=st.copy()
        jac=[]
        for j,c in enumerate(charts,start=k):
            f=expit(st[j]@c.B)
            out[j]=c.I@f
            if derivatives:
                jac.append((c.I*(f*(1-f)))@c.B.T)
        return out,jac

    np.testing.assert_allclose(decode(anchor,False)[0],y0,rtol=2e-12,atol=2e-14)
    def totals(z):
        Y,_=decode(z,False)
        _,_,_,h,s=M.weights(Y,True)
        _,r,b,*_=C.scalar(z[C.j0:-2])
        return h+r,s+b,float(Y[:,0].sum())

    calls=0
    def rhs(t,flat):
        nonlocal calls
        calls+=1
        Y,jac=decode(flat,True)
        z=np.r_[Y.ravel(),flat[C.j0:]]
        d=Cold.rhs(C,z,True)
        E,F,G,*_=M.weights(Y,True)
        for p in add:
            dr,S,W,*_=p.collision(E,F,G)
            d[:C.j0]+=dr.ravel()
            d[-2:]+=[S,W]
        if elastic is not None:
            dr,S,W,_=elastic.collision(G,F,G)
            d[:C.j0]+=dr.ravel()
            d[-2:]+=[S,W]
        dy=d[:C.j0].reshape(y0.shape)
        for j,J in enumerate(jac,start=k):
            x=dy[j].copy()
            dy[j]=np.linalg.solve(J,x)
            if calls%8==1:
                coordinate_defect.append(float(np.max(abs(J@dy[j]-x))))
                cond.append(float(np.linalg.cond(J)))
        return d
    initial=anchor
    before=0.
    checkstate=checkpoint_input is not None or emit_checkpoint or adaptive is not None
    provenance=None
    if checkstate:
        if any((prepared_input is not None, retune_local_daughter, compensation_only,
                number_only_local_retune, retune_donor_support is not None)):
            raise ValueError('Checkpoint pilot only retains the native and exact local archive preparations.')
        H=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
        modules={}
        for name,mod in sorted(sys.modules.items()):
            src=getattr(mod,'__file__',None)
            if name.startswith('early_pole_') and src and src.endswith('.py'):
                modules[name]=H(src)
        inputs={}
        for grid in {saved['grid'],saved['aperture_grid']}:
            inp=f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
            inputs[grid]=H(inp)
        provenance=dict(initializer=H(bundle),program=H(__file__),state_schema='joined-logit-f64-v1',module_sources=modules,
                        spectra=inputs,python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__,
                        longdouble_bits=int(np.finfo(np.longdouble).nmant),projection=projection,depth=depth,
                        charged_swap=swap,elastic_channel=elastic_channel,without_elastic=without_elastic,
                        local_daughter_fraction=local_daughter_fraction,rtol=rtol,base_atol=base_atol,
                        step_bound=duration/panels)
        if checkpoint_input is not None:
            prior=json.loads(Path(checkpoint_input).read_text())['checked_end_state']
            if prior['provenance']!=provenance:
                raise ValueError('Checkpoint came from different source, solver, or matched scope.')
            initial=np.asarray(prior['state'])
            before=float(prior['time'])
            if initial.shape!=anchor.shape or not np.all(np.isfinite(initial)) or not before>0.:
                raise ValueError('Checkpoint has invalid state or exposure.')
            ng=initial[C.j0:-2]
            if not np.all(ng>0.):
                raise ValueError('Invalid restarted cold occupations.')
            E,S,N=totals(initial)
            A,B,N0=totals(anchor)
            if (abs(E-A)>3e-12 or abs(E-A-initial[-1])>3e-12
                    or abs(S-B-initial[-2])>3e-11 or abs(N-N0)>3e-12):
                raise ValueError('Restart input does not reproduce the original global ledgers.')
    e0,s0,nH=totals(initial)
    atol=np.r_[np.full(2*k,base_atol),np.full(C.j0-2*k,2e-3*rtol),
               np.full(len(n0),base_atol/10),base_atol,base_atol/10]
    sample=before+np.linspace(0.,duration,panels+1)
    if adaptive is None:
        result=solve_ivp(rhs,(before,before+duration),initial,method='DOP853',rtol=rtol,atol=atol,
                         max_step=duration/panels,t_eval=sample)
    else:
        from types import SimpleNamespace
        from accepted_dop853_checkpoint_reference import evolve,dense_matrix
        mode,path,cadence,stop=adaptive
        if checkpoint_input is not None:
            raise ValueError('Accepted-step restart cannot be combined with the legacy end-state input')
        path=Path(path)
        def save():
            return dict(calls=calls,coordinate_defect=coordinate_defect,condition=cond)
        def restore(payload):
            nonlocal calls
            count=payload['calls']
            err,conditions=payload['coordinate_defect'],payload['condition']
            if not isinstance(count,int) or count<1 or not np.all(np.isfinite(err+conditions)):
                raise ValueError('Invalid stateful diagnostics in accepted-step file')
            calls=count
            coordinate_defect[:]=err
            cond[:]=conditions
        mode_opts=dict(save_model_state=save,save_every=cadence,
                       stop=stop if mode=='suspend' else None,
                       start=path if mode=='resume' else None,
                       load_model_state=restore if mode=='resume' else None)
        tag=dict(provenance=provenance,record_mode_history=record_mode_history)
        control=evolve(rhs,(before,before+duration),initial,sample,duration/panels,rtol,atol,
                       tag,path,**mode_opts)
        if mode=='suspend':
            if control is not None:
                raise ArithmeticError('Requested adaptive save did not happen')
            return dict(suspended=True)
        result=SimpleNamespace(y=dense_matrix(control),success=True,message='validated resume')
    if not result.success:
        raise RuntimeError(result.message)
    z=result.y[:,-1]
    e1,s1,n1=totals(z)
    en,sn=e1-e0,s1-s0-(z[-2]-initial[-2])
    assert abs(en)<3e-12 and abs(en-(z[-1]-initial[-1]))<3e-12
    assert abs(sn)<3e-11 and abs(n1-nH)<3e-12
    cold_modes=result.y[C.j0:-2]
    count=cold_modes.sum(axis=0)
    modes=(dict(cold_mode_count_history=cold_modes.T.tolist(),
                cold_momentum_nodes_over_T=C.q.tolist())
           if record_mode_history else {})
    return dict(projection=projection,duration=duration,report_panels=panels,rtol=rtol,base_atol=base_atol,
                calls=calls,mean_cold_rate=float((count[-1]-count[0])/duration),
                cold_count_history=count.tolist(),history_exposure=sample.tolist(),**modes,
                largest_sample_coordinate_residual=max(coordinate_defect),
                largest_sample_chart_condition_number=max(cond),
                energy_defect=en,entropy_defect=sn,hard_number_defect=n1-nH,
                **(dict(prepared_charged_swap=swap) if swap is not None else {}),
                **(dict(complement_retune=retune_check) if retune_check is not None else {}),
                **(dict(elastic_channels=False) if without_elastic else {}),
                **(dict(selected_elastic_channel=elastic_channel) if elastic_channel!='both' else {}),
                **(dict(aperture_subdivision=depth) if depth!=1 else {}),
                **(dict(local_archive_support=local_archive_meta) if local_archive_meta is not None else {}),
                **(dict(checked_end_state=dict(provenance=provenance,time=float(sample[-1]),state=z.tolist())) if emit_checkpoint else {}),
                scope='Unsplit frozen unlike and same-flavor response; finite explicit time reference in charged logits.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    p.add_argument('--duration',type=float,default=3.3333333333333335e-7)
    p.add_argument('--panels',type=int,default=16)
    p.add_argument('--rtol',type=float,default=3e-9)
    p.add_argument('--base-atol',type=float,default=2e-14)
    group=p.add_mutually_exclusive_group()
    group.add_argument('--prepared-input')
    group.add_argument('--retune-local-daughter',action='store_true')
    group.add_argument('--compensation-only',action='store_true')
    group.add_argument('--number-only-local-retune',action='store_true')
    group.add_argument('--local-from-paired-archive',action='store_true')
    group.add_argument('--retune-donor-support',choices=('below','no_neighbors'))
    elast=p.add_mutually_exclusive_group()
    elast.add_argument('--without-elastic',action='store_true')
    elast.add_argument('--elastic-channel',choices=('both','unlike','same'),default='both')
    p.add_argument('--record-mode-history',action='store_true')
    p.add_argument('--local-daughter-fraction',type=float,default=1.)
    p.add_argument('--checkpoint-input')
    p.add_argument('--emit-checkpoint',action='store_true')
    p.add_argument('--adaptive-mode',choices=('run','suspend','resume'))
    p.add_argument('--adaptive-path',type=Path)
    p.add_argument('--adaptive-save-every',type=int,default=1)
    p.add_argument('--adaptive-stop-after',type=int,default=1)
    a=p.parse_args()
    opts=vars(a)
    mode,path=opts.pop('adaptive_mode'),opts.pop('adaptive_path')
    gap,stop=opts.pop('adaptive_save_every'),opts.pop('adaptive_stop_after')
    if min(a.duration,a.panels,a.rtol,a.base_atol,gap,stop)<=0.:
        p.error('Use positive scope, adaptive cadence and accuracy controls.')
    if (mode is None)!=(path is None):
        p.error('Adaptive operation needs both --adaptive-mode and --adaptive-path')
    adaptive=None if mode is None else (mode,path,gap,stop)
    print(json.dumps(run(**opts,adaptive=adaptive),indent=2,allow_nan=False))
