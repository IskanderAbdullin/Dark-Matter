"""Exact two-packet realization of stored cohort birth moments."""
import numpy as np

from frw_exterior_cohort_reference import Exterior
from frw_exterior_restart_reference import export_restart,import_restart
from joint_prompt_two_moment_evolution_reference import Cells
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import m,mp
from warm_prompt_high_parent_ledger_reference import high
from warm_prompt_reentry_shape_audit_reference import get_mom,pad


def mom(s):
    return get_mom(np.asarray(s["lab_momentum_edges_GeV"]),
                   [s["entropy_variable_scalar"],s["entropy_variable_fermion"]])


def seed(s,hard):
    p=np.asarray(s["lab_momentum_edges_GeV"])
    j0=len(p)-1
    wide=np.r_[p,p[-1]*np.geomspace(1.,2.,6)[1:]]
    assert p[-1]<.5<hard<wide[-1]
    th=pad(s,wide)
    ext=import_restart(s["positive_exterior_restart"])
    a,dm,_,side,birth_floor,n,mean=ext.array()
    assert np.all(side==1)
    # A lower packet sits above both the original birth edge and the
    # present q32 truncation. Stay below the mean, without placing all
    # low moment weight at a shell endpoint.
    floor=np.maximum(birth_floor,np.hypot(dm,p[-1]/a))
    low=floor+.2*(mean-floor)
    upper=np.hypot(dm,hard/a)
    assert np.all((floor<mean)&(mean<upper))
    frac=(mean-low)/(upper-low)
    wn=np.array([n*(1-frac),n*frac])
    E0=np.array([low,upper])
    np.testing.assert_allclose((wn*E0).sum(axis=0),n*mean,rtol=4e-15,atol=0.)
    pp=a*np.sqrt((E0-dm)*(E0+dm))
    E=np.hypot(pp,dm)
    shell=np.searchsorted(wide,pp,side="right")-1
    assert np.all((shell>=j0)&(shell<len(wide)-1))
    ledger=[]
    for j,d in enumerate((mp,m)):
        ok=dm==d
        keys=shell[:,ok].ravel()
        NN=np.bincount(keys,weights=wn[:,ok].ravel(),minlength=len(wide)-1)
        KK=np.bincount(keys,weights=(wn[:,ok]*pp[:,ok]**2/(E[:,ok]+d)).ravel(),minlength=len(wide)-1)
        C=Cells(wide,d,(1,4)[j])
        t,err,_=C.reconstruct([NN,KK])
        assert err<1e-9
        th[j][:,j0:]=np.asarray(t)[:,j0:]
        ledger.append(dict(added_count=NN.sum(),hard_packet_count=wn[1,ok].sum(),reconstruction_error=err))
    N,K,Eout=get_mom(wide,th)
    ref,K2,E2=mom(s)
    extN=np.array([n[dm==d].sum() for d in (mp,m)])
    np.testing.assert_allclose(N.sum(axis=1),ref.sum(axis=1)+extN,rtol=4e-14,atol=0.)
    initial_energy=E2.sum()+ext.energy(1.)
    realized=(wn*E).sum()
    assert ext.energy(1.)<=realized*(1+2e-14) and realized<=ext.upper_bound(1.)*(1+2e-14)
    assert abs(Eout.sum()/(E2.sum()+realized)-1)<2e-12
    S=dict(s,lab_momentum_edges_GeV=wide.tolist(),entropy_variable_scalar=th[0].tolist(),
           entropy_variable_fermion=th[1].tolist(),
           represented_species_count_a3=(N.sum(axis=1)*s["scale_factor"]**3).tolist(),
           positive_exterior_restart=export_restart(Exterior()))
    return S,dict(cohorts=len(n),anchor_hard_p_GeV=hard,relative_total_energy_above_jensen=float(Eout.sum()/initial_energy-1.),
                  cohort_energy_between_Jensen_realization_and_upper=[ext.energy(1.),realized,ext.upper_bound(1.)],
                  per_species=ledger)


def highrate(s,b,eos):
    _,_,ee=mom(s)
    R=ee.sum()
    h=np.sqrt((TabulatedRadiation(eos,s["T_visible_GeV"]).state(1.)[1]+R)/3)/MPL
    vals=np.array([high(s,1./b,n)[1:3] for n in (8,14)])/(R*h)
    assert np.all(vals>0.)
    assert np.max(abs(vals[-1]/vals[0]-1))<.01
    return vals.tolist()


