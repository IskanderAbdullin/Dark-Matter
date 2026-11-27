"""Unsplit finite early joined reference in admissible charged-logit coordinates."""
import argparse
import gzip
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


def run(bundle,projection,duration,panels,rtol,base_atol,prepared_input=None,retune_local_daughter=False,compensation_only=False,number_only_local_retune=False,without_elastic=False,retune_donor_support=None,elastic_channel='both'):
    with gzip.open(bundle,'rt') as f:
        branches=json.load(f)['branches']
    saved=branches[projection]
    shot=saved['snapshots'][0]
    y0=np.asarray(shot['hard_moments'])
    n0=np.asarray(shot['cold_bins'])
    C=Screened(saved['grid'],1,len(n0),duration,(1,1,2),.08,1.,saved['aperture_grid'],projection)
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
    e0,s0,nH=totals(anchor)
    atol=np.r_[np.full(2*k,base_atol),np.full(C.j0-2*k,2e-3*rtol),
               np.full(len(n0),base_atol/10),base_atol,base_atol/10]
    sample=np.linspace(0.,duration,panels+1)
    result=solve_ivp(rhs,(0.,duration),anchor,method='DOP853',rtol=rtol,atol=atol,
                     max_step=duration/panels,t_eval=sample)
    if not result.success:
        raise RuntimeError(result.message)
    z=result.y[:,-1]
    e1,s1,n1=totals(z)
    en,sn=e1-e0,s1-s0-z[-2]
    assert abs(en)<3e-12 and abs(en-z[-1])<3e-12
    assert abs(sn)<3e-11 and abs(n1-nH)<3e-12
    count=result.y[C.j0:-2].sum(axis=0)
    return dict(projection=projection,duration=duration,report_panels=panels,rtol=rtol,base_atol=base_atol,
                calls=calls,mean_cold_rate=float((count[-1]-count[0])/duration),
                cold_count_history=count.tolist(),history_exposure=sample.tolist(),
                largest_sample_coordinate_residual=max(coordinate_defect),
                largest_sample_chart_condition_number=max(cond),
                energy_defect=en,entropy_defect=sn,hard_number_defect=n1-nH,
                **(dict(prepared_charged_swap=swap) if swap is not None else {}),
                **(dict(complement_retune=retune_check) if retune_check is not None else {}),
                **(dict(elastic_channels=False) if without_elastic else {}),
                **(dict(selected_elastic_channel=elastic_channel) if elastic_channel!='both' else {}),
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
    group.add_argument('--retune-donor-support',choices=('below','no_neighbors'))
    elast=p.add_mutually_exclusive_group()
    elast.add_argument('--without-elastic',action='store_true')
    elast.add_argument('--elastic-channel',choices=('both','unlike','same'),default='both')
    a=p.parse_args()
    if min(a.duration,a.panels,a.rtol,a.base_atol)<=0.:
        p.error('Use positive scope and accuracy controls.')
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
