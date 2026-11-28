"""Separate scalar abundance and shape response under fixed-peak clock controls."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from archived_collisionless_streaming_reference import Background
from empty_lepton_pair_spectrum_reference import Source
from finite_cooling_relativistic_los_reference import trajectory
from kloe_pion_band_portal_reference import table
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M
from warm_direct_vacuum_scalar_psd_reference import pion_parents,scalar_cumulative
from warm_direct_vacuum_fixed_metric_response_reference import kernel,layout
from warm_light_occupation_stage_audit_reference import read


def build(eos,kloe,edges,n):
    bg=Background(eos)
    rad=TabulatedRadiation(eos,.05)
    data=table(kloe)
    src=Source((12,40,28),tail=75.,dark_mass=M)
    x,w=roots_legendre(n)
    all_C=np.zeros((4,2,len(edges)))
    for lo,hi in zip((.012,.025,.035,.04,.045),(.025,.035,.04,.045,.05),strict=True):
        a,b=np.log([lo,hi]);h=(b-a)/2
        for T,dln in zip(np.exp((a+b)/2+h*x),w*h,strict=True):
            A=bg.scale(T)
            _,rho=rad.state(A/bg.scale(.05))
            ent=2*np.pi**2/45*T**3*np.exp(bg.gs(np.log(T)))
            d=(3+bg.gs.derivative()(np.log(T)))/3
            omega=dln*d/(ent*np.sqrt(rho/3)/MPL)
            alt=np.array([1.,1/np.sqrt(1+(T/.03)**2),1/np.sqrt(1+(T/.02)**2),.2 if T>.045 else 1.])
            parents=(src.particles(T,42),pion_parents(T,data,7,32,22))
            C=np.asarray([scalar_cumulative(E,W,edges*bg.T0/A) for E,W in parents])
            all_C+=alt[:,None,None]*C*omega
    return all_C


def run(eos,kloe,source,metric,prior,reference):
    bg=Background(eos)
    N=read(source)
    P=read(prior)
    input_edge=N["conserved_momentum_edges_q_over_T0"]
    edges=np.array([x if x is not None else np.inf for x in input_edge])
    j=int(np.flatnonzero(edges==.1)[0])
    fine=build(eos,kloe,edges,12)
    coarse=build(eos,kloe,edges,8)
    assert np.max(abs(fine[:, :, -1]/coarse[:, :, -1]-1))<2e-4
    assert np.max(abs(fine[:, :, j]/coarse[:, :, j]-1))<.003
    old=np.asarray(N["controls"][-1]["conditional_total_scalar_yield"])
    assert np.max(abs(fine[0,:,-1]/old-1))<.001
    target=np.asarray([r["controls"][-1]["combined_clock_control_ratios"] for r in P["cases"]])[0]
    true_cold=fine[:,:,j].sum(axis=1)/fine[0,:,j].sum()
    assert max(abs(target-true_cold))<.002
    Y=fine.sum(axis=1)
    h=np.maximum(np.diff(Y,axis=1),0.)
    H=h/h.sum(axis=1)[:,None]
    rate=Y[:,-1]/Y[0,-1]
    f=read(reference)["held_and_actual_direct_scalar_rest_fractions_with_same_archived_stable_count"][1]
    fclock=rate*f/(1-f+rate*f)
    local_cdf=np.max(abs(np.cumsum(H-H[0],axis=1)),axis=1)
    ep=np.where(np.isfinite(edges),edges,36.)
    q,w=layout(ep,12,"q")
    out=[]
    for raw in metric:
        m=trajectory(raw)
        E,_,_,K=kernel(bg,m,q,13000)
        weight=(E*K).reshape(len(ep)-1,12)@w
        norm=E.reshape(len(ep)-1,12)@w
        R=(H@weight)/(H@norm)
        shape=f*(R-R[0])
        actual=fclock*(R-1)-f*(R[0]-1)
        np.testing.assert_allclose(actual,(fclock-f)*(R[0]-1)+fclock*(R-R[0]),atol=2e-20)
        out.append(dict(k_per_Mpc=m[0],scalar_shape_only_response_over_CDM=R.tolist(),
                        scalar_response_change_vs_standard=(R-R[0]).tolist(),
                        normalized_matter_proxy_difference_shape_only=shape.tolist(),
                        normalized_matter_proxy_difference_shape_and_count=actual.tolist()))
    return dict(scope="Same thermal source kernels, peak temperature and entropy-to-q map. Alternate clock "
                      "changes only production duration; final scalar uses same fixed metric. No "
                      "self-consistent alternative expansion history or physical transfer.",
                scalar_rest_fractions_in_same_stable_background=fclock.tolist(),clock_labels=P["clock_columns"],combined_total_yield=Y[:,-1].tolist(),
                yield_ratios_vs_standard=(Y[:,-1]/Y[0,-1]).tolist(),
                equality_cut_cold_yield_ratios_vs_standard=true_cold.tolist(),
                largest_normalized_CDF_shift_vs_standard=local_cdf.tolist(),metrics=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","kloe","source","prior","reference"):
        p.add_argument(k,type=Path)
    p.add_argument("metric",nargs="+",type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.kloe,a.source,a.metric,a.prior,a.reference),indent=2,allow_nan=False))
