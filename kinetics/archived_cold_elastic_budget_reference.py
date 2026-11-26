"""Instantaneous archived cold-shell elastic budget on paired supports."""
import argparse,json
from pathlib import Path
import numpy as np

from compact_elastic_panel_reference import CompactCombined
from continuum_dark_normal_capacity_reference import from_row
from joint_below_pole_elastic_linear_modes_reference import setup
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from late_archived_shape_continuation_reference import Late
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from three_species_decay_onsager_reference import m,mp


def test(eos,table,row,shape):
    base=from_row(row,eos,20)
    n=np.array(base["species_counts_GeV3"])
    T=base["common_T_GeV"]
    rho=TabulatedRadiation(table,base["actual_visible_T_GeV"]).state(1.)[1]
    H=np.sqrt((rho+np.array([mp,m])@n+base["retained_dark_kinetic_energy_GeV4"])/3)/MPL
    edges=np.array(shape["lab_momentum_edges_GeV"])
    th=[np.array(shape["entropy_variable_"+s]) for s in ("scalar","fermion")]
    boundary=16.*row["initial_visible_T_GeV"]/row["final_a"]
    nr=int(np.argmin(abs(edges-boundary)))
    assert nr<len(edges)-1 and abs(edges[nr]-boundary)<5e-13
    outcomes=[]
    for label,q,state in (("wide",edges,th),("paired inner only",edges[:nr+1],[t[:,:nr] for t in th])):
        model=Late(Collision(q),H,[],False)
        U=model.moment(pack(state))[0].reshape(2,model.n,2)
        ref,_R,Td=model.reference(pack(state),T)
        for order in ((3,2,2,8),(4,3,2,10),(6,4,3,12)):
            panels=[]
            setup(Td,model.k,Td,model.eta,*order,32.,.82,cache=panels,cache_only=True)
            b=CompactCombined(model.k,H,panels)
            b.t=Td;b.eta=model.eta
            F,_=b.elastic(pack(state)-ref)
            d=F.reshape(2,model.n,2)
            sums=d[:,:,0].sum(axis=1)
            energy=np.array([d[j,:,0]@(c.mass+c.low)+d[j,:,1]@c.dk for j,c in enumerate(model.k.c)])
            norms=np.array([np.sum(abs(d[:,:,0])),sum(np.sum(abs(d[j,:,0])*(c.mass+c.low)+abs(d[j,:,1])*c.dk)
                                     for j,c in enumerate(model.k.c))])
            defect=[float(np.max(abs(sums))/norms[0]),float(abs(energy.sum())/norms[1])]
            assert max(defect)<2e-11
            cp=model.k.c[0]
            k0=U[0,0,0]*cp.low[0]+U[0,0,1]*cp.dk[0]
            kdot=d[0,0,0]*cp.low[0]+d[0,0,1]*cp.dk[0]
            outcomes.append(dict(label=label,order=order,panels=len(panels),
                    cold_number_fractional_slope_per_Hubble=float(d[0,0,0]/U[0,0,0]),
                    cold_kinetic_fractional_slope_per_Hubble=float(kdot/k0),
                    scalar_energy_slope_over_total_dark_kinetic_per_Hubble=float(energy[0]/base["retained_dark_kinetic_energy_GeV4"]),
                    invariant_defects=defect))
    return dict(rows=outcomes,Hubble_GeV=float(H),
      scope="common H for both domains; physical archived shapes within each momentum cutoff; instantaneous below-pole slopes, not evolved rates")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos");p.add_argument("table");p.add_argument("run");p.add_argument("shape")
    a=p.parse_args()
    r=json.loads(Path(a.run).read_text())["rows"][-1]
    s=json.loads(Path(a.shape).read_text())["snapshots"][-1]
    print(json.dumps(test(StandardModel(a.eos),a.table,r,s),indent=2,allow_nan=False))
