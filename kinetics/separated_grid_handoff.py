"""Late equal-number FD closure on the SAME resolved collision grid.

No remapping of scalar or chi2 and no switching to a different collision
kernel. The parent is replaced by the equal-number discrete FD distribution;
bath temperature is adjusted to preserve the instantaneous total energy.
The reduced equations preserve parent+scalar and chi2-scalar source balances
and FRW continuity. Handoff criteria are local, not a global accuracy bound.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import expit
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from parent_angle_kernel import ParentAngle,SeparatedRectangle
from laine_portal_rates import PortalTable,make_models
from radiation_fixed_bath import omega_from_number
from hybrid_parent_closure import handoff_passes


class GridClosure:
    def __init__(self,model):
        self.m=model;self.n=model.n;self.nd=2*model.n

    def thermo(self,x,N,zt):
        m=self.m;m.geometry(x);theta=np.exp(zt);E=m.e1;v=m.v1
        if N>0:
            xi=brentq(lambda u:v@expit(u-E/theta)-N,x/theta-750,x/theta+80,xtol=1e-12)
            g=expit(xi-E/theta);s=g*(1-g);A=v@s
        else:
            xi=-np.inf;g=np.zeros_like(E);s=g;A=0.
        if A>1e-250:
            rn=float((v*s)@E/A)
            gn=s/A
            gt=s*(E-rn)/theta
            ht=float((v*E)@gt)
            force=float((v*m.t1)@gt)
        else:
            b=np.exp((x-E)/theta);gn=b/(v@b);rn=float((v*E)@gn)
            gt=np.zeros_like(g);ht=force=0.
        rb=m.rad0*theta**4;rho=float((v*E)@g)
        return dict(g=g,gn=gn,gt=gt,rhoN=rn,heat=4*rb+ht,force=force,theta=theta,rb=rb,rho=rho)

    def evaluate(self,z,y,jac=False):
        m=self.m;n=self.nd;N,zt,W,S1,S2=y[n:];th=self.thermo(np.exp(z),N,zt)
        fy=np.r_[y[:n],th['g'],zt,W,S1,S2]
        f=m.rhs(z,fy)
        D=f[:n];dn=float(m.v1@f[n:m.dim])
        dt=(th['force']-th['rhoN']*dn-m.energy_weights[:n]@D)/th['heat']
        out=np.r_[D,dn,dt,f[-3:]]
        if not jac:return out
        J=m.jac(z,fy)
        # Chain rule for the reduced state. Conservation rows use the exact
        # full-grid source derivatives; two temperature derivatives are FD.
        def reduced(A):
            return np.column_stack((A[:,:n],A[:,n:m.dim]@th['gn'],
                A[:,n:m.dim]@th['gt']+A[:,-4],A[:,-3:]))
        JD=reduced(J[:n])
        JN=reduced((m.v1@J[n:m.dim])[None,:])[0]
        JT=(-th['rhoN']*JN-m.energy_weights[:n]@JD)/th['heat']
        for k in [n,n+1]:
            eps=max(1e-10,abs(y[k])*1e-5) if k==n else 1e-6
            a=y.copy();b=y.copy();a[k]+=eps;b[k]-=eps
            if k==n and b[k]<0:JT[k]=(self.evaluate(z,a)[n+1]-out[n+1])/eps
            else:JT[k]=(self.evaluate(z,a)[n+1]-self.evaluate(z,b)[n+1])/(2*eps)
        return np.vstack((JD,JN,JT,reduced(J[-3:])))

    def rhs(self,z,y):return self.evaluate(z,y)
    def jac(self,z,y):return self.evaluate(z,y,True)

    def validate(self):
        m=self.m;x=.5;th=self.thermo(x,0.,0.)
        N=float(m.v1@expit(-m.e1))*.8
        rng=np.random.default_rng(73)
        y=np.r_[rng.uniform(.01,.3,self.nd),N,.01,0.,0.,0.]
        dr=rng.normal(size=len(y));dr/=np.linalg.norm(dr);eps=1e-6
        a=(self.rhs(np.log(x),y+eps*dr)-self.rhs(np.log(x),y-eps*dr))/(2*eps)
        b=self.jac(np.log(x),y)@dr
        err=float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300));assert err<1e-5,err
        return dict(jacobian_relative_error=err)


def model_from(row,table):
    cls=ParentAngle if row['kernel']=='parent-angle' else SeparatedRectangle
    kw=dict(n=row['n'],nparent=row['nparent'],scalar_ratio=row['scalar_ratio'],fermion_ratio=row['fermion_ratio'],
        kappa=row['kappa'],kappa_sm1=row['kappa_sm1'],kappa2=row['kappa2'],g_bath=row['g_bath'],
        g_fermion=row.get('g_fermion',4),qmin=row['frozen_cell_spectrum']['edges'][0],
        qmax=row['frozen_cell_spectrum']['edges'][-1],parent_grid=row['parent_grid'],
        fermion_grid=row['fermion_grid'],parent_floor=row['parent_floor'],gamma=row['gamma'],scatter=row['scatter'])
    if table:
        _,R=make_models()
        class T(cls,R):pass
        return T(table=table,parent_GeV=row['parent_GeV'],**kw)
    return cls(**kw)


def run(row,table=None,xend=80.,rtol=2e-6,require_guard=True):
    f0=row['final'];x0=f0['x'];spec=row['frozen_cell_spectrum']
    if not xend>x0:raise ValueError('handoff requires unfinished integration')
    if row.get('rate_source_sha256') and table is None:raise ValueError('Original published rate table required')
    m=model_from(row,table);C=GridClosure(m);checks=C.validate();m.geometry(x0)
    assert np.allclose(m.edges,spec['edges'],rtol=1e-12,atol=0)
    assert np.allclose(m.fermion_edges,spec['chi2_edges'],rtol=1e-12,atol=0)
    assert np.allclose(m.parent_edges,spec['parent_edges'],rtol=1e-12,atol=0)
    p=np.asarray(spec['f_phi']);f=np.asarray(spec['f_chi2']);g=np.asarray(spec['f_parent'])
    yd=np.r_[p,f];N=float(m.v1@g);Np=float(m.vol@p);N2=float(m.v2@f)
    theta=f0['theta'];target=f0['total_rho_Tref4']
    t0=C.thermo(x0,N,np.log(theta))
    def mismatch(zt):
        th=C.thermo(x0,N,zt)
        return th['rb']+th['rho']+m.energy_weights[:2*m.n]@yd-target
    zt=brentq(mismatch,np.log(theta)-.1,np.log(theta)+.1,xtol=1e-14)
    th=C.thermo(x0,N,zt)
    old=m.dark(np.r_[yd,g])[:m.n];new=m.dark(np.r_[yd,th['g']])[:m.n]
    a=float(m.vol@old);b=float(m.vol@new)
    diag=dict(x=x0,temperature_relative_change=float(np.exp(zt)/theta-1),
        shape_TV=float(.5*m.v1@np.abs(g-t0['g'])/N),
        parent_energy_relative_difference=float((m.v1*m.e1)@g/t0['rho']-1),
        relative_source_mismatch=float((b-a)/max(abs(a),abs(b),1e-30)),
        l1_source_mismatch_over_sum=float(m.vol@np.abs(new-old)/max(m.vol@(np.abs(new)+np.abs(old)),1e-30)),
        energy_matching_residual=float(mismatch(zt)/target))
    diag['accuracy_guard_passed']=handoff_passes(diag)
    print('same_grid_handoff',json.dumps(diag),flush=True)
    if require_guard and not diag['accuracy_guard_passed']:raise ValueError('same-grid handoff guard failed')
    y0=np.r_[yd,N,zt,0.,0.,0.]
    ts=sorted(set([x0]+[x for x in [3,4,5,6,8,10,15,20,30,40,60,80] if x0<x<xend]+[xend]))
    atol=np.r_[np.full(m.n,1e-12),np.full(m.n,1e-13),1e-14,1e-13,1e-12,1e-13,1e-13]
    sol=solve_ivp(C.rhs,[np.log(x0),np.log(xend)],y0,method='BDF',jac=C.jac,rtol=rtol,
        atol=atol,max_step=.035,t_eval=np.log(ts))
    if not sol.success:raise RuntimeError(sol.message)
    history=[];n=m.n;l,h=m.edges[:-1],m.edges[1:]
    for x,y in zip(ts[1:],sol.y.T[1:]):
        pp=y[:n];ff=y[n:2*n];N1,zt,W,S1,S2=y[2*n:];th=C.thermo(x,N1,zt);theta=th['theta']
        s=m.snapshot(y[:2*n]);s.pop('number_mismatch',None)
        np_=s['number'];n2=float(m.v2@ff);rho1=th['rho']
        rho2=float((m.v2*m.energy_chi)@ff);rhop=float((m.vol*m.energy_phi)@pp)
        total=th['rb']+rho1+rho2+rhop
        e1=(N1+np_-N-Np-S1)/max(np_,1e-12);e2=(n2-np_-(N2-Np)-S2)/max(np_,1e-12)
        cont=(total-target-W)/target
        cold=np.maximum(np.minimum(h,.1*theta)**3-l**3,0)/(h**3-l**3)
        assert max(abs(e1),abs(e2))<1e-6 and abs(cont)<5e-5
        assert s['min_occupation']>-1e-7 and s['max_fermion']<1.00001 and N1>-1e-9
        s.update(x=x,theta=float(theta),mother_number=float(N1),mother_rho_Tref4=rho1,
            fermion_number=n2,rho_chi2_Tref4=rho2,rho_phi_Tref4=rhop,total_rho_Tref4=float(total),
            continuity_residual=float(cont),number_source_parent_residual=float(e1),number_source_residual=float(e2),
            integrated_parent_bath_source=float(S1+f0['integrated_parent_bath_source']),
            integrated_fermion_bath_source=float(S2+f0['integrated_fermion_bath_source']),
            matter_fraction_chi2=float(rho2/total),mean_p_over_actual_T=float(s['mean_q']/theta),
            hot_fraction_p_over_actual_T_gt_0p1=float((pp*m.vol)@(1-cold)/np_),
            omega_phi_h2_m86=float(omega_from_number(np_/theta**3,gs=m.g_bath)),phase='same_grid_closure')
        history.append(s)
    out={k:v for k,v in row.items() if k not in ('history','final','frozen_cell_spectrum','success','message','xend')}
    out.update(hybrid=True,fully_resolved_until=x0,xend=xend,success=True,closure_checks=checks,
        closure_nfev=sol.nfev,handoff=diag,history=[dict(s,phase='resolved') for s in row['history']]+history,
        final=history[-1],frozen_cell_spectrum=dict(edges=m.edges.tolist(),chi2_edges=m.fermion_edges.tolist(),
            f_phi=sol.y[:n,-1].tolist(),f_chi2=sol.y[n:2*n,-1].tolist(),g_phi=1,g_chi2=m.g_fermion,
            theta_freezeout=history[-1]['theta']))
    return out


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--table',default='.runtime/mesh_imPiR.dat');p.add_argument('--xend',type=float,default=80.)
    p.add_argument('--exploratory',action='store_true');a=p.parse_args()
    r=json.loads(Path(a.input).read_text())[0]
    t=PortalTable(a.table) if r.get('rate_source_sha256') else None
    out=run(r,t,a.xend,require_guard=not a.exploratory)
    Path(a.output).write_text(json.dumps([out],indent=2,allow_nan=False)+'\n')
    print('final',json.dumps(out['final']),flush=True)

if __name__=='__main__':main()
