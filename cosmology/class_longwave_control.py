"""Diagnose unacceptable large-scale normalization of the two-tracer CLASS run."""
import argparse,json,tempfile
from pathlib import Path
import numpy as np
from energy_closed_class import canonical
from kinetic_psd_class import export_components


def main():
    from classy import Class
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    raw=json.loads(Path(a.input).read_text())[0];row=canonical(raw)
    k=np.array([.001,.003,.01,.03,.1,.3,.5,1.,1.5]);z=0.
    base=dict(h=.674,omega_b=.02237,T_cmb=2.7255,N_ur=3.044+row['final']['delta_Neff_chi2_if_stable_free_streaming'],
              A_s=2.1e-9,n_s=.965,tau_reio=.0544,output='mPk,mTk',z_max_pk=.1,P_k_max_h_Mpc=3.)
    base['P_k_max_h/Mpc']=base.pop('P_k_max_h_Mpc')
    def solve(p):
        c=Class()
        try:
            c.set(p);c.compute();h=float(c.h);bg=c.get_background();b={v:float(bg[v][-1]) for v in bg if v.startswith('(.)rho_')}
            tr=c.get_transfer(z);kt=np.asarray(tr['k (h/Mpc)'])*h
            T={v:np.interp(np.log(k),np.log(kt),tr[v]).tolist() for v in tr if v.startswith('d_') or v in ('phi','psi')}
            return dict(power=[float(c.pk_lin(float(x),z)) for x in k],omega_m=float(c.Omega_m)*h*h,
                        omega_ncdm=float(c.Omega_nu)*h*h,Neff=float(c.Neff),sigma8=float(c.sigma8()),background=b,transfer=T)
        finally:
            c.struct_cleanup();c.empty()
    with tempfile.TemporaryDirectory() as folder:
        e=export_components(row,folder);cold,hot=e['components'];T=(3.909/106.75)**(1/3)
        ncdm=dict(N_ncdm=2,m_ncdm='86,86',T_ncdm=','.join(str(v['scale']*T) for v in [cold,hot]),deg_ncdm='1,1',
                   use_ncdm_psd_files='1,1',ncdm_psd_filenames=','.join(v['file'] for v in [cold,hot]),ncdm_quadrature_strategy='2,2',ncdm_N_momentum_bins='120,120',omega_cdm=1e-8)
        out={}
        for gauge in ['newtonian','synchronous']:
            b=dict(base,gauge=gauge);warm=solve(dict(b,**ncdm));om=warm['omega_ncdm'];out[gauge+'_two_tracer']=warm
            out[gauge+'_CDM']=solve(dict(b,omega_cdm=om+1e-8))
            frac=hot['number']/(cold['number']+hot['number']);oh=om*frac
            hp=dict(N_ncdm=1,m_ncdm='86',T_ncdm=str(T),deg_ncdm='1',use_ncdm_psd_files='1',ncdm_psd_filenames=hot['file'],
                    ncdm_quadrature_strategy='2',ncdm_N_momentum_bins='120',omega_cdm=om-oh+1e-8)
            out[gauge+'_hybrid']=solve(dict(b,**hp))
            # Same physical background but use finite positive uniform-q quadrature on the cold component.
            qp=dict(ncdm,ncdm_quadrature_strategy='3,2',ncdm_N_momentum_bins='400,120',ncdm_maximum_q='12.,45.')
            out[gauge+'_uniform_cold']=solve(dict(b,**qp))
        ratios={}
        for key,v in out.items():
            b=out[key.split('_')[0]+'_CDM'];ratios[key]=(np.array(v['power'])/b['power']).tolist()
            print(key,json.dumps(dict(omega_m=v['omega_m'],omega_ncdm=v['omega_ncdm'],Neff=v['Neff'],sigma8=v['sigma8'],ratios=ratios[key])),flush=True)
    Path(a.output).write_text(json.dumps(dict(k_per_Mpc=k.tolist(),runs=out,ratios=ratios),indent=2)+'\n')

if __name__=='__main__':main()
