"""Frozen prompt, elastic and pair current on the exact late reactivation phases."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import expit

from joint_below_pole_elastic_linear_modes_reference import setup
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from late_radiation_background_reference import MPL, TabulatedRadiation
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import observable_weights
from warm_coupled_serial_restart_reference import atomic, sha
from warm_light_occupation_stage_audit_reference import digest, read
from warm_pair_common_fermion_chart_reference import rate
from warm_prompt_reentry_shape_audit_reference import get_mom


class ElasticLedger:
    def __init__(self, theta, cut):
        self.theta = theta
        self.cut = cut
        self.sum = np.zeros(3)
        self.projection = np.zeros_like(theta)
        self.panels = self.high_lab_panels = 0
        self.high_lab = np.zeros(3)
        self.max_lab = 0.

    def append(self, item):
        x, y, lab, P, raw = item
        n = len(self.theta)//2
        a, b = P[:, :n] @ self.theta[:n], P[:, n:] @ self.theta[n:]
        assert np.min(a)>0.
        block = expit(b)/(-np.expm1(-a))
        log = -a-b
        ends = np.maximum(x-self.cut, 0.)+np.maximum(y-self.cut, 0.)
        change = ends[None, :]-ends[:, None]
        delta = log[:, None]-log[None, :]
        mat = raw*block[:, None]*block[None, :]
        ordered = mat*np.exp(log[:, None])
        net = mat*np.sign(delta)*np.exp(np.maximum(log[:, None], log[None, :]))*(-np.expm1(-abs(delta)))
        v = np.array([np.sum(ordered*np.maximum(change, 0.)),
                      np.sum(ordered*np.maximum(-change, 0.)), .5*np.sum(net*change)])
        self.sum += v
        self.projection -= P.T@net.sum(axis=1)
        self.panels += 1
        self.max_lab = max(self.max_lab, lab)
        if lab>self.cut:
            self.high_lab += v
            self.high_lab_panels += 1


def run(eos, tapes, out, channel, order_mode):
    rows = []
    A = read(tapes[0])
    phase0 = A["snapshot"]
    assert A["hash"]==digest(phase0)
    assert A["plan"]["eos_sha256"]==sha(eos)
    T = phase0["T_visible_GeV"]
    p = np.array(phase0["lab_momentum_edges_GeV"])
    keys = ("entropy_variable_scalar", "entropy_variable_fermion")
    rho = get_mom(p, [np.array(phase0[k]) for k in keys])[2].sum()
    H = np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+rho)/3.)/MPL
    k = Collision(p)
    op = LocalAbsolutePrompt(k,H,[],True)
    cut = A["plan"]["normalized_cut_GeV"]
    weights = observable_weights(op,cut)
    limit = sum(np.hypot(cell.mass,p[-1]) for cell in k.c)
    common_order = [(3,2,2,8),(4,3,3,12)][order_mode]
    pair_order = [(3,7,6,2),(4,10,10,4)][order_mode]
    for path in tapes:
        r = read(path)
        S = r["snapshot"]
        assert r["hash"]==digest(S) and r["plan"]==A["plan"]
        assert np.array_equal(S["lab_momentum_edges_GeV"],p)
        th = pack([S[key] for key in keys])
        f = (weights@op.collision(th)[0])/rho
        if r["signed_native_flux"]:
            assert abs(f/r["signed_native_flux"]-1.) < 2e-11
        row = dict(time=r["time"], input_sha256=r["hash"], prompt_signed_over_rhoH=float(f))
        if channel == "elastic":
            controls = []
            for extent in (32., limit/T):
                ledger = ElasticLedger(th,cut)
                setup(T,k,T,(0.,0.),*common_order,extent,.82,cache=ledger,cache_only=True,
                      energy_observable_cuts=(cut,))
                norm = max(np.max(ledger.sum[:2]),1e-290)
                assert abs(ledger.sum[0]-ledger.sum[1]-ledger.sum[2])<2e-12*norm
                assert abs(weights@ledger.projection-ledger.sum[2])<2e-10*norm
                # Closed support: this is the exact energy-conservation
                # ceiling for either two-particle side, not a thermal tail.
                assert 0. < ledger.max_lab <= limit*(1.+2e-14)
                controls.append(dict(boost_extent_in_T=float(extent),max_lab_GeV=float(ledger.max_lab),
                                     panels=ledger.panels,panels_above_cut=ledger.high_lab_panels,
                                     rates_over_rhoH=(ledger.sum/(rho*H)).tolist(),
                                     lab_above_cut_rates_over_rhoH=(ledger.high_lab/(rho*H)).tolist()))
            assert controls[0]["rates_over_rhoH"] == [0.,0.,0.]
            row.update(elastic=controls)
        else:
            B,_,F,_ = snapshot_functions(S)
            excess = {}
            _,_,gross,err,high = rate(p,F,B,pair_order,excess_cut=(cut,),excess_result=excess)
            row.update(pair_rates_over_rhoH=(np.array(excess["rate_GeV5"])/(rho*H))[0].tolist(),
                       pair_columns=excess["columns"],pair_ledger_defect=float(err),
                       gross_pair_energy_rates_over_rhoH=(gross[:,1]/(rho*H)).tolist(),
                       empty_outgoing_fermion_limit_GeV=float(high))
        rows.append(row)
        atomic(out,dict(channel=channel,order=common_order if channel=="elastic" else pair_order,
                        normalization_rho_GeV4=float(rho),H_GeV=float(H),cut_GeV=cut,
                        scope="Frozen reactivated phase: closed resolved elastic, open-outgoing pair, closed prompt.",
                        rows=rows))
        print(channel,order_mode,r["time"],file=sys.stderr,flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("channel",choices=("elastic","pair"))
    p.add_argument("tapes",type=Path,nargs="+")
    p.add_argument("--order-mode",type=int,choices=(0,1),default=0)
    a=p.parse_args()
    run(a.eos,a.tapes,a.output,a.channel,a.order_mode)
