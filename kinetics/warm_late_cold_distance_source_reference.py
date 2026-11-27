"""Ask how much of the direct vacuum lepton source can have a small travelled distance."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp
from archived_collisionless_streaming_reference import Background
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import P,mp
from warm_late_cold_scalar_prompt_partition_reference import counts
from warm_late_cold_visible_offpole_dose_reference import source
from warm_late_cold_visible_recoil_box_upper_reference import log_slope_minimum
from warm_light_occupation_stage_audit_reference import digest,read


def q_limit(bg, ai, length):
    def r(x,n):
        return float(bg.streams(np.array([x*bg.T0]),mp,ai,bg.aeq,n)[0])
    root=brentq(lambda z:r(z,5)/length-1,1e-8,10.,xtol=1e-14)
    lower,upper=root*np.array([1-2e-7,1+2e-7])
    for n in (4,7,10):
        assert r(lower,n)<length<r(upper,n)
    return root*bg.T0,upper*bg.T0


def analytic_late_q_enclosure(bg, ending_birth, length, n):
    # Drop motion before the last EOS knot. From there the reference
    # radiation factor is constant and the rest of a^4 H^2 increases.
    ak=np.exp(bg.v[0])
    assert ending_birth<=ak<bg.aeq
    a=np.geomspace(ak,bg.aeq,n+1)
    denom=np.sqrt(bg.or0+bg.om*a[1:]+bg.ol*a[1:]**4)
    def lower(x):
        q=x*bg.T0
        d=q/mp*np.diff(np.arcsinh(mp*a/q))
        return bg.D0*np.sum(d/denom)
    x=brentq(lambda z:lower(z)/length-1,1e-8,10.,xtol=1e-14)*(1+2e-7)
    assert lower(x)>length
    return float(x*bg.T0)


def run(eos,archive):
    bg=Background(eos)
    monotonic=[log_slope_minimum(bg.gs,3),log_slope_minimum(bg.gr,4)]
    z=read(archive);S=z["snapshot"]
    assert z["snapshot_sha256"]==digest(S)
    Ta=S["T_visible_GeV"]
    aa=bg.scale(Ta)
    ends=[.05,.0001]
    alate=bg.scale(ends[1])
    s_anchor=(2*np.pi**2/45)*Ta**3*np.exp(bg.gs(np.log(Ta)))
    edges=np.unique(np.r_[[0,.5,1,2,4,7,9,10,np.inf],np.linspace(0,10,41)])
    rows=[]
    for D in (.0556,.1,.3,1.):
        a,am=q_limit(bg,aa,D)
        b,check=q_limit(bg,alate,D)
        caps=[analytic_late_q_enclosure(bg,alate,D,n) for n in (64,256,1024)]
        bm=caps[-1]
        assert caps[0]>=caps[1]>=bm>check>b>a
        held=float(counts(S,np.array([a/aa]),24)[0]/s_anchor)
        vals=[]
        for n in (96,384):
            T=np.geomspace(*ends,n+1)
            A=np.array([bg.scale(t) for t in T])
            p=bm/A[:-1]
            assert np.all(p<P)
            ang=p*p/(2*P*(P+np.sqrt((P-p)*(P+p))))
            t=T[1:];u=np.log(t)
            denom=(2*np.pi**2/45)*t**5*np.exp(bg.gs(u)+.5*bg.gr(u))*np.sqrt(np.pi**2/90)/MPL
            sub=np.array([source(t,q,edges) for t,q in zip(T[:-1],p,strict=True)])
            log=float(logsumexp(np.log(2*np.diff(np.log(A))*ang/denom)+sub))
            vals.append(dict(time_bins=n,log10_direct_lepton_yield_upper=log/np.log(10.),
                             log10_upper_relative_to_held_distance_selected_yield=(log-np.log(held))/np.log(10.)))
        assert vals[-1]["log10_direct_lepton_yield_upper"]<=vals[0]["log10_direct_lepton_yield_upper"]
        rows.append(dict(max_path_Mpc=D,archive_momentum_cut_over_T=a/(aa*Ta),
                         conservative_late_placement_cut_q_over_today_T=bm/bg.T0,
                         upper_q_over_today_T_using_post_EOS_only=[d/bg.T0 for d in caps],
                         numeric_late_placement_crossing_q_over_today_T=b/bg.T0,
                         archive_held_distance_selected_yield=held,controls=vals))
    return dict(scope="Reference ballistic travel to equality, not a transfer scale. Later birth allows a shorter "
                      "path at fixed comoving momentum, so use the latest tested birth as a conservative "
                      "momentum cut for the whole 50 to .1 MeV direct thermal charged-lepton vacuum interval. "
                      "No initial parent, hadronic, later scattering or growing perturbations covered.",
                archive=z["snapshot_sha256"],monotone_EOS_lower=monotonic,rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path);p.add_argument("archive",type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.archive),indent=2,allow_nan=False))
