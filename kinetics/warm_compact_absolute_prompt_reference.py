"""Assemble the direct-occupation Jacobian within each parent phase panel."""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from compact_prompt_jacobian_reference import compile_prompt,iter_prompt
from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_frw_split_reference import pack
from late_radiation_background_reference import MPL,TabulatedRadiation
from warm_direct_prompt_stress_reference import DirectPrompt
from warm_light_occupation_stage_audit_reference import read
from warm_prompt_reentry_shape_audit_reference import get_mom


class LocalAbsolutePrompt(DirectPrompt):
    def __init__(self,k,H,panels,prompt,stream_local=False):
        super().__init__(k,H,panels,prompt)
        assert prompt and not panels
        self.parents=None if stream_local else compile_prompt(self)

    def local_panels(self):
        return iter_prompt(self) if self.parents is None else iter(self.parents)

    def collision(self,q):
        k=self.k
        B,F=self.arrays(q)
        f=k.c[0].f(B,k.ph,k.I)
        g=k.c[1].f(F,k.ch,k.J)
        a=k.pre*f*g
        b=k.pre*(1+f)*(1-g)
        def sum_parent(z):
            return np.bincount(k.y,weights=z,minlength=k.n)[k.y]
        A=sum_parent(a); B=sum_parent(b); G=A+B
        r=np.divide(k.fac,G*self.H,out=np.zeros_like(G),where=G>0.)
        val=b*A-a*B
        out=self.proj(r*val)
        jac=None
        def derive(z):
            nonlocal jac
            if jac is None:
                jac=np.zeros((self.d,self.d))
                for s,ix,P,sf in self.local_panels():
                    p=a[s,None]
                    q=b[s,None]
                    f0=f[s,None]
                    g0=g[s,None]
                    da=-p*P*np.where(sf,1+f0,1-g0)
                    db=q*P*np.where(sf,-f0,g0)
                    dA=da.sum(axis=0)
                    dB=db.sum(axis=0)
                    x=s.start
                    rat=(dA+dB)/G[x] if G[x]>0 else np.zeros_like(dA)
                    moving=db*A[x]+q*dA-da*B[x]-p*dB-val[s,None]*rat
                    jac[np.ix_(ix,ix)]+=P.T@(r[s,None]*moving)
            return jac@z
        return out,derive


def run(eos,paths):
    out=[]
    for p in paths:
        row=read(p)
        S=row["snapshot"]
        edge=np.asarray(S["lab_momentum_edges_GeV"])
        th=[np.asarray(S["entropy_variable_"+k]) for k in ("scalar","fermion")]
        U=pack(th)
        rho=get_mom(edge,th)[2].sum()
        H=np.sqrt((TabulatedRadiation(eos,S["T_visible_GeV"]).state(1.)[1]+rho)/3)/MPL
        ker=Collision(edge)
        full=DirectPrompt(ker,H,[],True)
        local=LocalAbsolutePrompt(ker,H,[],True)
        values=[]
        ticks=[]
        for model in (full,local):
            t=time.perf_counter()
            x,j=model.collision(U)
            J=np.stack([j(v) for v in np.eye(model.d)],axis=1)
            ticks.append(time.perf_counter()-t)
            values.append((x,J))
        f,j=values[0]
        g,k=values[1]
        np.testing.assert_array_equal(f,g)
        whole=np.linalg.norm(k-j)/np.linalg.norm(j)
        cols=np.linalg.norm(k-j,axis=0)/np.maximum(np.linalg.norm(j,axis=0),1e-250)
        # Include rate contributions for initially exactly empty cells.
        v=np.sin(.55*np.arange(len(U))+.2)
        h=2e-5
        exact=(full.collision(U+h*v)[0]-full.collision(U-h*v)[0])/(2*h)
        fd=np.linalg.norm(exact-k@v)/np.linalg.norm(exact)
        assert whole<1e-12 and np.max(cols)<1e-10 and fd<4e-7
        out.append(dict(source_snapshot_sha256=row.get("hash"),shells=full.n,
                        rates_bit_for_bit_equal=bool(np.array_equal(f,g)),jacobian_whole_relative_error=whole,
                        max_column_relative_error=float(np.max(cols)),relative_finite_difference_check=fd,
                        per_call_timings_seconds_direct_and_local=ticks))
    return dict(scope="Same physical flux and direct entropy-variable Jacobian, not a new physics approximation. "
                      "Wall timings in one cloud process, not a resource guarantee.",cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("paths",type=Path,nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.paths),indent=2,allow_nan=False))
