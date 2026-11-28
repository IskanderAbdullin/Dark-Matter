"""Open the final-state support of the frozen below-pole elastic channel."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import expit

from joint_below_pole_elastic_linear_modes_reference import setup
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from late_radiation_background_reference import MPL, TabulatedRadiation
from warm_coupled_serial_restart_reference import atomic, sha
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom, pad


class Ledger:
    def __init__(self, theta, cut, boundary):
        self.theta, self.cut, self.boundary = theta, cut, boundary
        self.closed = np.zeros(3)
        self.opened = np.zeros(3)
        self.exterior = np.zeros((2,3))
        self.transition_energy = 0.
        self.incoming_outside = 0.
        self.exchange = np.zeros(3)
        self.maximum_support_output = np.zeros(2)
        self.panels = self.cross_panels = 0

    def append(self,data):
        E1,E2,L,F,raw = data
        self.panels += 1
        n = len(self.theta)//2
        a,b = F[:, :n]@self.theta[:n], F[:, n:]@self.theta[n:]
        assert np.min(a)>0.
        E = np.stack((E1,E2))
        high = E>self.boundary[:,None]
        old = ~np.any(high,axis=0)
        block = expit(b)/(-np.expm1(-a))
        lr = -a-b
        block[high[0]] = expit(b[high[0]])
        block[high[1]] = 1./(-np.expm1(-a[high[1]]))
        block[np.all(high,axis=0)] = 1.
        rate = raw*block[:,None]*block[None,:]*np.where(old,np.exp(lr),0.)[:,None]
        # Symmetric quadrature in the incoming and outgoing angle.
        assert np.all(rate[~old] == 0.)
        W = np.maximum(E-self.cut,0.)
        total_W = W.sum(axis=0)
        change = total_W[None,:]-total_W[:,None]
        value = np.stack((np.maximum(change,0.),np.maximum(-change,0.),change))
        self.closed += np.sum(value*rate*old[None,:],axis=(1,2))
        self.opened += np.sum(value*rate,axis=(1,2))
        if np.all(old):
            return
        self.cross_panels += 1
        out = rate.sum(axis=0)
        tail = out[None,:]*high
        self.exterior[:,0] += np.sum(tail,axis=1)
        self.exterior[:,1] += np.sum(E*tail,axis=1)
        self.exterior[:,2] += np.sum(W*tail,axis=1)
        self.transition_energy += L*np.sum(out[~old])
        self.incoming_outside += L*np.sum(rate[~old])
        for j,x in enumerate(tail):
            if np.any(x>0.):
                self.maximum_support_output[j] = max(self.maximum_support_output[j],np.max(E[j,x>0.]))
        event_to_outside=rate[:,~old]
        self.exchange += [np.sum(event_to_outside), np.sum(event_to_outside*change[:,~old]),
                          np.sum(event_to_outside*total_W[:,None])]


def run(eos,tapes,out,order_mode):
    rows=[]
    for path in tapes:
        r=read(path)
        S=r["snapshot"]
        assert r["hash"]==digest(S) and r["plan"]["eos_sha256"]==sha(eos)
        p=np.asarray(S["lab_momentum_edges_GeV"])
        T=S["T_visible_GeV"]
        z=pad(S,p)
        rho=get_mom(p,z)[2].sum()
        H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+rho)/3.)/MPL
        base=Collision(p)
        edge=np.array([np.hypot(c.mass,p[-1]) for c in base.c])
        L=edge.sum()
        # Even an outgoing massless particle cannot exceed the total
        # energy of two resolved incoming particles. Extra shells are
        # exact vacuum inputs and allow the entire outgoing range.
        q=np.r_[p,np.geomspace(p[-1],np.nextafter(L,np.inf),7)[1:]]
        original=pack(z)
        extended=pack(pad(S,q))
        k=Collision(q)
        cut=r["plan"]["normalized_cut_GeV"]
        assert cut<edge.min()
        assert np.all(get_mom(q,pad(S,q))[0][:,len(p)-1:]==0.)
        order=((3,2,2,8),(4,3,3,12))[order_mode]
        led=Ledger(extended,cut,edge)
        setup(T,k,T,(0.,0.),*order,L/T,.82,cache=led,cache_only=True,energy_observable_cuts=(cut,))
        assert led.incoming_outside==0.
        gross=max(led.opened[:2])
        assert abs(led.opened[0]-led.opened[1]-led.opened[2])<2e-12*max(gross,1e-290)
        assert abs(led.opened[2]-led.closed[2]-led.exchange[1]) < 2e-12*max(gross,1e-290)
        assert np.all(led.maximum_support_output<=L)
        # A separate uncropped-old-grid root verifies comparison
        # sensitivity to quadrature changes from the novel p boundaries.
        old=Ledger(original,cut,edge)
        setup(T,base,T,(0.,0.),*order,L/T,.82,cache=old,cache_only=True,energy_observable_cuts=(cut,))
        assert np.array_equal(old.closed,old.opened) and not np.any(old.exterior)
        rows.append(dict(hard_GeV=r["plan"]["hard_GeV"],time=r["time"],source_hash=r["hash"],
                         rho_GeV4=float(rho),H_GeV=float(H),cut_GeV=cut,
                         old_and_new_shells=[len(p)-1,len(q)-1],
                         physical_incoming_lab_ceiling_GeV=float(L),
                         old_chart_gross_and_signed=(old.closed/(rho*H)).tolist(),
                         matched_chart_closed_gross_signed=(led.closed/(rho*H)).tolist(),
                         inclusive_open_gross_signed=(led.opened/(rho*H)).tolist(),
                         newly_outside_by_species_N_E_W_over_rhoH=(led.exterior/(rho*H)).tolist(),
                         maximum_sampled_outgoing_energy_GeV=led.maximum_support_output.tolist(),
                         total_lab_energy_of_events_with_outgoing_exterior_over_rhoH=float(led.transition_energy/(rho*H)),
                         crossing_events_N_signed_W_incoming_W_over_rhoH=(led.exchange/(rho*H)).tolist(),
                         panels=led.panels,eligible_panels=led.cross_panels,
                         comparison_prompt_signed_over_rhoH=r["signed_native_flux"]))
        atomic(out,dict(order=order,scope="Frozen vacuum-exterior outgoing test, not re-entry or time evolution.",rows=rows))
        print(r["plan"]["hard_GeV"],r["time"],order,file=sys.stderr,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("tapes",type=Path,nargs="+")
    p.add_argument("--order-mode",type=int,choices=(0,1),default=0)
    a=p.parse_args()
    run(a.eos,a.tapes,a.output,a.order_mode)
