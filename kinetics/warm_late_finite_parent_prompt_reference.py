"""Evolve the parent on the native prompt chart, eliminating it in each implicit step."""
import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.special import xlogy

from joint_prompt_two_moment_evolution_reference import Collision
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import M,P,width
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import observable_weights
from warm_coupled_serial_restart_reference import atomic, sha
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


class FiniteParent(LocalAbsolutePrompt):
    def __init__(self,k,H):
        super().__init__(k,H,[],True)
        p=np.sqrt((k.E-M)*(k.E+M))
        assert np.all(p>0.)
        self.relax_factor=width*M*M/(2.*P*k.E*p*H)
        self.measure=k.w/(H*self.relax_factor)
        self.old=np.zeros(k.n)
        self.dt=0.

    def reactions(self,u):
        k=self.k
        n,g=self.arrays(u)
        f=k.c[0].f(n,k.ph,k.I)
        g=k.c[1].f(g,k.ch,k.J)
        a=k.pre*f*g
        b=k.pre*(1.+f)*(1.-g)
        A=np.bincount(k.y,weights=a,minlength=k.n)
        B=np.bincount(k.y,weights=b,minlength=k.n)
        return f,g,a,b,A,B

    def equilibrium(self,u):
        _,_,_,_,A,B=self.reactions(u)
        return np.divide(A,A+B,out=np.zeros_like(A),where=A+B>0.)

    def occupation(self,u):
        _,_,_,_,A,B=self.reactions(u)
        s=self.dt*self.relax_factor
        return (self.old+s*A)/(1.+s*(A+B))

    def collision(self,u):
        k=self.k
        f,g,a,b,A,B=self.reactions(u)
        s=self.dt*self.relax_factor[k.y]
        old=self.old[k.y]
        A,B=A[k.y],B[k.y]
        den=1.+s*(A+B)
        V=s*(b*A-a*B)+(a+b)*old-a
        coef=k.fac/(self.H*den)
        ans=self.proj(coef*V)
        jac=None
        def derive(z):
            nonlocal jac
            if jac is None:
                jac=np.zeros((self.d,self.d))
                for panel,ix,G,sf in self.parents:
                    ap,bp=a[panel,None],b[panel,None]
                    da=-ap*G*np.where(sf,1.+f[panel,None],1.-g[panel,None])
                    db=bp*G*np.where(sf,-f[panel,None],g[panel,None])
                    dA,dB=da.sum(axis=0),db.sum(axis=0)
                    q=panel.start
                    dv=s[q]*(db*A[q]+bp*dA-da*B[q]-ap*dB)+(da+db)*old[q]-da
                    Jac=coef[panel,None]*(dv-V[panel,None]*(s[q]/den[q])*(dA+dB))
                    jac[np.ix_(ix,ix)]+=G.T@Jac
            return jac@z
        return ans,derive

    def heavy_entropy(self,h):
        return self.measure@(-xlogy(h,h)-xlogy(1.-h,1.-h))


def run(eos,paths,out,initial):
    raw=[read(p) for p in paths]
    assert all(x["hash"]==digest(x["snapshot"]) and x["plan"]==raw[0]["plan"] for x in raw)
    assert raw[0]["plan"]["eos_sha256"]==sha(eos)
    def theta(S):
        return [S["entropy_variable_"+k] for k in ("scalar","fermion")]
    S=raw[0]["snapshot"]
    T=S["T_visible_GeV"]
    p=np.array(S["lab_momentum_edges_GeV"])
    C=raw[0]["plan"]["normalized_cut_GeV"]
    rho=get_mom(p,theta(S))[2].sum()
    H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+rho)/3.)/MPL
    m=FiniteParent(Collision(p),H)
    E=m.k.E
    w=m.measure
    W=observable_weights(m,C)
    u=pack(theta(S))
    m.old=m.equilibrium(u) if initial=="equilibrium" else np.zeros_like(w)
    h0=m.old.copy()
    rows=[]
    def totals(u,hv):
        th=m.arrays(u)
        n,_,e=get_mom(p,th)
        return n.sum(axis=1)+w@hv,e.sum()+(w*E)@hv,m.moment(u)[1]+m.heavy_entropy(hv)
    N0,E0,S0=totals(u,h0)
    for a,b in zip(raw,raw[1:]):
        assert np.array_equal(p,b["snapshot"]["lab_momentum_edges_GeV"]) and b["snapshot"]["T_visible_GeV"]==T
        m.dt=b["time"]-a["time"]
        assert m.dt>0.
        Tbefore=totals(u,m.old)[2]
        start=W@m.moment(u)[0]/rho
        guess=pack(theta(b["snapshot"]))
        if not rows:
            # Test the actual condensed residual derivative at the first
            # accepted root's guess before applying the time step.
            z=.12*np.sin(.57*np.arange(m.d)+.55)
            dx=1e-5
            f,j=m.collision(guess)
            df=(m.collision(guess+dx*z)[0]-m.collision(guess-dx*z)[0])/(2*dx)
            jerr=float(np.linalg.norm(df-j(z))/np.linalg.norm(df))
            assert jerr<2e-7
        def entropy_guard(v,after,before):
            return after+m.heavy_entropy(m.occupation(v))>=before+m.heavy_entropy(m.old)
        v,log=advance(m,u,m.dt,T,initial_theta=guess,vacuum_scaled=True,balanced_merit=True,
                      entropy_guard=entropy_guard)
        next_h=m.occupation(v)
        assert np.all((next_h>=0.)&(next_h<=1.))
        N,En,Sn=totals(v,next_h)
        assert Sn>=Tbefore-3e-15*S0
        led=np.r_[N/N0-1.,En/E0-1.]
        assert max(abs(led))<2e-12
        last=W@m.moment(v)[0]/rho
        f=W@m.collision(v)[0]/rho
        defect=last-start-m.dt*f
        assert abs(defect)<1e-10*max(last,start,1e-250)
        # Parent equation weighted directly, independently of the daughter projection.
        _,_,_,_,AA,BB=m.reactions(v)
        dh=m.dt*m.relax_factor*((1.-next_h)*AA-next_h*BB)
        closure=(next_h-m.old)-dh
        scale=np.sum(w*abs(next_h-m.old))
        # Relative-to-change errors here amplify cancellation by the
        # stiff response. Judge the direct equation against rounding of
        # its two large RHS terms, independent of the analytical root.
        residue=w@abs(closure)
        cancellation_budget=np.finfo(float).eps*(w@(abs(next_h)+abs(m.old)+
             m.dt*m.relax_factor*((1.-next_h)*AA+next_h*BB)))
        parent_precision=residue/max(scale,1e-260)
        rounding_ratio=residue/max(cancellation_budget,1e-260)
        assert rounding_ratio<12.
        prompt=float(W@m.moment(guess)[0]/rho)
        rows.append(dict(time=b["time"],actual_W=float(last),stationary_W=prompt,
                         actual_signed_flux=float(f),stationary_signed_flux=b.get("signed_native_flux"),
                         total_invariant_relative_defects=led.tolist(),parent_weighted_residual_relative=float(parent_precision),
                         parent_weighted_residual_over_binary64_rounding_scale=float(rounding_ratio),
                         tail_projection_residual=float(defect),solver=log,seeding_phase_hash=b["hash"],
                         total_entropy_change_relative=float(Sn/S0-1.),
                         signed_heavy_energy_change_over_rho=float((w*E)@(next_h-h0)/rho)))
        u=v
        m.old=next_h
        out.mkdir(exist_ok=True,parents=True)
        state=dict(S,**{"entropy_variable_"+k:r.tolist() for k,r in zip(("scalar","fermion"),m.arrays(u),strict=True)})
        atomic(out/f"step-{len(rows)}.json",dict(snapshot=state,hash=digest(state),row=rows[-1],
                                              parent_occupations=next_h.tolist()))
        print(initial,b["time"],log["newton_iterations"],file=sys.stderr,flush=True)
    atomic(out/"summary.json",dict(scope="Fixed-background, closed, finite parent chart, not a general thermodynamic trajectory.",
                                  initial_parent=initial,initial_hash=raw[0]["hash"],derivative_relative_error=jerr,
                                  parent_energy_nodes=len(E),eos_sha256=sha(eos),rows=rows))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("paths",type=Path,nargs="+")
    p.add_argument("--initial",choices=("equilibrium","vacuum"),default="equilibrium")
    a=p.parse_args()
    run(a.eos,a.paths,a.output,a.initial)
