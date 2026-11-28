"""Localize the positive cold source in temperature and change only its exposure clock."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from archived_collisionless_streaming_reference import Background
from kloe_pion_band_portal_reference import table
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL
from warm_late_cold_distance_source_reference import q_limit
from warm_late_cold_visible_exact_recoil_fold_reference import local
from warm_late_cold_visible_kloe_recoil_reference import rates
from warm_light_occupation_stage_audit_reference import read


def run(eos,tablefile,eqdata,paths):
    bg=Background(eos)
    rad=TabulatedRadiation(eos,.05)
    g,cs,e=table(tablefile)
    eq=read(eqdata)
    d=read(paths)
    compare={"eq_p_over_T_0.1":[eq["rows"][-1]["controls"][-1]["FD_yield"],None],
             "travel_at_most_0.3_Mpc":[d["rows"][0]["controls"][-1]["lepton_FD"],d["rows"][0]["controls"][-1]["pion_BE"]],
             "travel_at_most_1_Mpc":[d["rows"][2]["controls"][-1]["lepton_FD"],d["rows"][2]["controls"][-1]["pion_BE"]]}
    Tedge=np.array([.012,.025,.035,.040,.045,.0475,.049,.05])
    rows=[]
    for name,D in (("eq_p_over_T_0.1",None),("travel_at_most_0.3_Mpc",.3),("travel_at_most_1_Mpc",1.)):
        c=[]
        for nt,no in ((5,8),(9,12)):
            x,w=roots_legendre(nt)
            moments=[]
            for t0,t1 in zip(Tedge[:-1],Tedge[1:],strict=True):
                a,b=np.log([t0,t1])
                jac=(b-a)/2
                out=[]
                for T in np.exp((a+b)/2+jac*x):
                    A=bg.scale(T)
                    _,rho=rad.state(A/bg.scale(.05))
                    S=2*np.pi**2/45*T**3*np.exp(bg.gs(np.log(T)))
                    H=np.sqrt(rho/3)/MPL
                    dt=(3+bg.gs.derivative()(np.log(T)))/3
                    q=.1*bg.T0 if D is None else q_limit(bg,A,D)[0]
                    p=q/A
                    base=np.array([local(T,p,no)[0],rates(T,p,g,cs,e,no)[1]])*dt/(H*S)
                    # Fixed entropy and momentum/temperature mapping, only
                    # the interval's proper time is modified in these controls.
                    clock=np.array([1.,1/np.sqrt(1+(T/.03)**2),1/np.sqrt(1+(T/.02)**2),
                                    .2 if T>.045 else 1.])
                    out.append(clock[:,None]*base)
                moments.append(jac*np.einsum('i,ijk->jk',w,np.asarray(out)))
            result=np.asarray(moments)
            net=result.sum(axis=0)
            for i,v in enumerate(compare[name]):
                if v is not None:
                    assert abs(net[0,i]/v-1)<.004
            assert np.all(net>0.) and np.all(net[1:]<=net[0])
            c.append(dict(time_rule=nt,rate_rule=no,
                standard_FD_lepton_and_BE_measured_pion_yield=net[0].tolist(),
                rate_fraction_per_interval_standard=(result[:,0,:]/net[0]).tolist(),
                combined_fraction_from_top5MeV=float(result[4:,0,:].sum()/net[0].sum()),
                combined_fraction_from_top2p5MeV=float(result[5:,0,:].sum()/net[0].sum()),
                combined_fraction_from_top1MeV=float(result[-1,0,:].sum()/net[0].sum()),
                combined_clock_control_ratios=(net.sum(axis=1)/net[0].sum()).tolist()))
        assert abs(c[-1]["combined_fraction_from_top5MeV"]-c[0]["combined_fraction_from_top5MeV"])<.001
        rows.append(dict(selection=name,controls=c))
    return dict(scope="Positive vacuum direct FD lepton plus measured-window BE pion source, 12–50 MeV. "
                      "All alternate histories only change the amount of proper time at fixed temperatures; "
                      "no new reheating or dark collision model solved.",
                KLOE_sha256=hashlib.sha256(tablefile.read_bytes()).hexdigest(),temperature_edges_GeV=Tedge.tolist(),
                clock_columns=["radiation","stiff_extra_fluid_equal_at_30_MeV","stiff_extra_fluid_equal_at_20_MeV",
                               "imposed_20_percent_residence_above_45_MeV"],cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","kloe","eqdata","paths"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.kloe,a.eqdata,a.paths),indent=2,allow_nan=False))
