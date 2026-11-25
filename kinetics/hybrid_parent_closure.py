"""Controlled but approximate handoff from momentum-resolved to number closure.

The switch preserves scalar and chi2 PSDs, parent number and total energy.
It does not claim a fully momentum-resolved late universe. Quantify shape,
energy/temperature and instantaneous collision mismatches at the handoff.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit
from yukawa_rate_bath import YukawaRateBath
from resolved_parent_bath import ResolvedParent
from radiation_fixed_bath import omega_from_number
from laine_portal_rates import PortalTable,make_models

def handoff_passes(d):
    return (d['x']>=1.5 and d['shape_TV']<1e-3 and
            abs(d['parent_energy_relative_difference'])<1e-3 and
            abs(d['relative_source_mismatch'])<.01 and
            d['l1_source_mismatch_over_sum']<.005 and
            abs(d['energy_matching_residual'])<1e-10)

def run(row,table=None,xend=80.,rtol=2e-6,require_guard=False):
    spec=row['frozen_cell_spectrum'];f0=row['final'];x0=f0['x'];n=row['n']
    if not xend>x0: raise ValueError('Handoff requires unfinished production')
    kw=dict(n=n,qmin=spec['edges'][0],qmax=spec['edges'][-1],scalar_ratio=row['scalar_ratio'],
            fermion_ratio=row['fermion_ratio'],kappa=row['kappa'],kappa_sm1=row['kappa_sm1'],
            kappa2=row['kappa2'],g_bath=row['g_bath'],g_fermion=row.get('g_fermion',4),
            gamma=row['gamma'],scatter=row['scatter'])
    if row.get('rate_source_sha256'):
        if table is None: raise ValueError('Original rate table required')
        C,R=make_models();kw.update(table=table,parent_GeV=row['parent_GeV'])
    else:C,R=YukawaRateBath,ResolvedParent
    m=C(**kw);r=R(nparent=row['nparent'],parent_floor=row['parent_floor'],parent_grid=row.get('parent_grid','legacy'),**kw)
    m.geometry(x0);r.geometry(x0)
    p=np.asarray(spec['f_phi']);f=np.asarray(spec['f_chi2']);g=np.asarray(spec['f_parent'])
    assert np.allclose(m.edges,spec['edges'],rtol=1e-12,atol=0)
    assert np.allclose(r.parent_edges,spec['parent_edges'],rtol=1e-12,atol=0)
    N=float(r.v1@g);theta=f0['theta'];target=f0['total_rho_Tref4']
    # Preserve total energy, not just temperature; kinetic re-equilibration
    # transfers any parent energy mismatch into the common radiation bath.
    yd=np.r_[p,f]
    ztheta=brentq(lambda z:m.thermo(x0,N,z)['rb']+m.thermo(x0,N,z)['rho']+m.ew@yd-target,
                   np.log(theta)-.1,np.log(theta)+.1,xtol=1e-14)
    th=m.thermo(x0,N,ztheta);fit=expit(th['xi']-r.e1/th['theta'])
    fit_same_temp=expit(brentq(lambda u:r.v1@expit(u-r.e1/theta)-N,-750,100)-r.e1/theta)
    Np=float(m.vol@p);N2=float(m.g_fermion*m.vol@f)
    total=r.rad0*theta**4+r.energy_weights@np.r_[p,f,g]
    Hf=np.sqrt(r.rad0/total)
    old=r.dark(np.r_[p,f,g],Hf)
    old_eq=r.dark(np.r_[p,f,fit_same_temp],Hf)
    y0=np.r_[yd,N,ztheta,0.,0.]
    m.set_state(x0,y0[:-1]);new=m.collision_rhs(yd)
    a=float(m.vol@old[:n]);b=float(m.vol@old_eq[:n]);c=float(m.vol@new[:n])
    diagnostics=dict(x=x0,temperature_before=theta,temperature_after=float(np.exp(ztheta)),
        temperature_relative_change=float(np.exp(ztheta)/theta-1),
        shape_TV=float(.5*r.v1@np.abs(g-fit_same_temp)/N),
        parent_energy_relative_difference=float((r.v1*r.e1)@g/((r.v1*r.e1)@fit_same_temp)-1),
        dark_number_source_resolved=a,dark_number_source_same_grid_FD=b,dark_number_source_closure=c,
        relative_source_mismatch=float((c-a)/max(abs(c),abs(a),1e-30)),
        l1_source_mismatch_over_sum=float(m.vol@np.abs(new[:n]-old[:n])/max(m.vol@(np.abs(new[:n])+np.abs(old[:n])),1e-30)),
        energy_matching_residual=float((th['rb']+th['rho']+m.ew@yd-target)/target))
    diagnostics['accuracy_guard_passed']=handoff_passes(diagnostics)
    print('handoff',json.dumps(diagnostics),flush=True)
    if require_guard and not diagnostics['accuracy_guard_passed']:
        raise ValueError('Cross-kernel handoff accuracy guard failed')
    xs=sorted(set([x0]+[x for x in [.5,1,2,3,5,10,15,20,30,40,60,80] if x0<x<xend]+[xend]))
    atol=np.r_[np.full(n,1e-12),np.full(n,1e-15),1e-14,1e-13,1e-12,1e-13]
    sol=solve_ivp(m.rhs_full,[np.log(x0),np.log(xend)],y0,method='BDF',jac=m.jac_full,
                  rtol=rtol,atol=atol,max_step=.035,t_eval=np.log(xs))
    if not sol.success:raise RuntimeError(sol.message)
    history=[]
    l,h=m.edges[:-1],m.edges[1:]
    for x,y in zip(xs[1:],sol.y.T[1:]):
        m.geometry(x);yd=y[:2*n];p,f=np.split(yd,2);N,zt,W,S=y[-4:]
        th=m.thermo(x,N,zt);theta=th['theta'];s=m.snapshot(yd);s.pop('number_mismatch',None)
        Nphi=s['number'];Nf=float(m.g_fermion*m.vol@f);rho1=th['rho']
        rho2=float((m.g_fermion*f*m.vol)@m.energy_chi);rhop=float((p*m.vol)@m.energy_phi)
        total=th['rb']+rho1+rho2+rhop
        cold=np.maximum(np.minimum(h,.1*theta)**3-l**3,0)/(h**3-l**3)
        cont=(total-target-W)/target
        err=(Nf-Nphi-(N2-Np)-S)/max(Nphi,1e-12)
        s.update(x=x,theta=float(theta),mother_number=float(N),mother_rho_Tref4=rho1,
            fermion_number=Nf,rho_chi2_Tref4=rho2,rho_phi_Tref4=rhop,total_rho_Tref4=float(total),
            continuity_residual=float(cont),number_source_residual=float(err),
            integrated_fermion_bath_source=float(S+f0['integrated_fermion_bath_source']),
            matter_fraction_chi2=float(rho2/total),mean_p_over_actual_T=float(s['mean_q']/theta),
            hot_fraction_p_over_actual_T_gt_0p1=float((p*m.vol)@(1-cold)/Nphi),
            omega_phi_h2_m86=float(omega_from_number(Nphi/theta**3,gs=m.g_bath)),phase='number_closure')
        assert abs(cont)<5e-5 and abs(err)<1e-6 and s['min_occupation']>-1e-7 and s['max_fermion']<1.00001
        history.append(s)
    out={k:v for k,v in row.items() if k not in ('history','final','frozen_cell_spectrum','success','message','xend')}
    out.update(hybrid=True,fully_resolved_until=x0,xend=xend,success=True,closure_nfev=sol.nfev,
               handoff=diagnostics,history=[dict(s,phase='resolved') for s in row['history']]+history,
               final=history[-1],frozen_cell_spectrum=dict(edges=m.edges.tolist(),f_phi=sol.y[:n,-1].tolist(),
                 f_chi2=sol.y[n:2*n,-1].tolist(),g_phi=1,g_chi2=m.g_fermion,theta_freezeout=history[-1]['theta']))
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',nargs='+',required=True)
    p.add_argument('--table',default='.runtime/mesh_imPiR.dat');p.add_argument('--output',required=True)
    p.add_argument('--xend',type=float,default=80.);p.add_argument('--require-guard',action='store_true');a=p.parse_args()
    out=[]
    for path in a.input:
        row=json.loads(Path(path).read_text())[0]
        table=PortalTable(a.table) if row.get('rate_source_sha256') else None
        r=run(row,table,a.xend,require_guard=a.require_guard);out.append(r)
        print('final',path,json.dumps(r['final']),flush=True)
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
