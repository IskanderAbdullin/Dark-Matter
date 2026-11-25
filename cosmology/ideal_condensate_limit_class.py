"""Conditional sudden ideal Bose equilibration + early elastic decoupling.

Keep the total scalar number and kinetic energy fixed, assume saturation
with conserved number, put the zero mode into CDM, and free-stream only
the thermal excited component. This is a hypothetical alternative to the
original nonthermal PSD, not an evolved interacting scalar prediction.
"""
import argparse,json,tempfile
from pathlib import Path
from importlib.metadata import version,PackageNotFoundError
import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import CubicSpline
from kinetic_psd_class import transformed_quadrature_check

def run(row,eq,folder,bins=120):
    from classy import Class
    theta=row['final']['theta']
    omega=row['final']['omega_phi_h2_m86']
    warmfrac=eq['excited_number_fraction']
    T=eq['temperature_over_Tref']
    # Tref redshifts exactly from the production endpoint.
    tr=(3.909/106.75)**(1/3)
    Trel=T/theta*tr
    zeta=row['scalar_ratio']*row['xend']/T
    q=np.unique(np.r_[np.geomspace(1e-8,.1,2500),np.geomspace(.1,60,4500)])
    K=q*q/(np.hypot(q,zeta)+zeta)
    f=np.maximum(1/np.expm1(K),1e-100*np.exp(-q))/(2*np.pi)**3
    file=Path(folder)/'saturated.dat'
    np.savetxt(file,np.column_stack((q,f)),fmt='%.17e')
    # Physical q units relative to SM endpoint bath T (as in other exports).
    scale=T/theta
    moments=[float(4*np.pi*scale**(k+3)*trapezoid(q**(k+2)*f,q)) for k in (0,1,2)]
    c=dict(file=str(file),scale=scale,number=moments[0],energy=moments[1],q2=moments[2])
    expnum=row['final']['number']/theta**3*warmfrac
    assert abs(moments[0]/expnum-1)<1e-4
    quad=transformed_quadrature_check(dict(components=[c]),bins)
    base=dict(h=.674,omega_b=.02237,T_cmb=2.7255,N_ur=3.044,A_s=2.1e-9,n_s=.965,tau_reio=.0544,
              gauge='synchronous',output='mPk',z_max_pk=5.1)
    base['P_k_max_h/Mpc']=34.5
    p=dict(base,N_ncdm=1,m_ncdm='86',T_ncdm=str(Trel),deg_ncdm='1',
           use_ncdm_psd_files='1',ncdm_psd_filenames=str(file),omega_cdm=omega*(1-warmfrac))
    try:comm=version('classy-community')
    except PackageNotFoundError:comm=None
    if comm and int(comm.split('.')[0])>=26:
        p.update(quadrature_strategy_ncdm_standard='2',N_momentum_bins_ncdm_standard=str(bins),
                 N_momentum_bins_bg_ncdm_standard=str(bins))
    else:p.update(ncdm_quadrature_strategy='2',ncdm_N_momentum_bins=str(bins))
    ks=np.geomspace(.003,20.,125);zs=(0.,5.)
    def solve(par):
        cl=Class()
        try:
            cl.set(par);cl.compute()
            def s(k):
                v=getattr(cl,k);return float(v() if callable(v) else v)
            return dict(omega_ncdm=s('Omega_nu')*s('h')**2,Neff=s('Neff'),
                        power={str(z):np.array([cl.pk_lin(float(k),z) for k in ks]) for z in zs})
        finally:
            try:cl.struct_cleanup();cl.empty()
            except Exception:pass
    a=solve(p)
    assert abs(a['omega_ncdm']/(omega*warmfrac)-1)<.005,(a['omega_ncdm'],omega*warmfrac)
    ref=solve(dict(base,omega_cdm=p['omega_cdm']+a['omega_ncdm']))
    rat={z:(a['power'][z]/ref['power'][z]) for z in a['power']}
    for r in rat.values():
        assert np.all(np.isfinite(r)) and np.all(r>0) and np.max(abs(r[ks<.03]-1))<.02
    return dict(bins=bins,version=comm,omega_expected=omega,omega_condensate=p['omega_cdm'],
                omega_ncdm_actual=a['omega_ncdm'],thermal_fraction_expected=warmfrac,temperature_ratio=Trel,
                Neff=a['Neff'],quadrature=quad,thermal_moment_number_relative_error=moments[0]/expnum-1,
                k_per_Mpc=ks.tolist(),ratios={z:r.tolist() for z,r in rat.items()},
                summary={z:dict(P10=float(np.interp(np.log(10),np.log(ks),r)),minimum=float(r.min())) for z,r in rat.items()})
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--screen',required=True)
    ap.add_argument('--output',required=True);a=ap.parse_args()
    row=json.loads(Path(a.input).read_text())[0]
    screen=json.loads(Path(a.screen).read_text())[0]
    eq=screen['ideal_equilibrium_limit'];assert eq['saturated']
    result=dict(input=a.input,scope='Instant ideal saturation at production end followed by immediate decoupling; pressureless condensate, collisionless excited bosons, equal primordial adiabatic modes. Not a solution of finite-rate condensation, not valid at cross sections requiring later kinetic decoupling.',runs=[])
    with tempfile.TemporaryDirectory() as f:
        for n in (120,240):
            r=run(row,eq,f,n);result['runs'].append(r)
            Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
            print(n,r['summary'],r['omega_ncdm_actual'],r['Neff'],flush=True)
if __name__=='__main__':main()
