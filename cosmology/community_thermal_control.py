"""Independent reference check of CLASS-community 26 against archived CLASS 3.4.3."""
import argparse,json,tempfile
from pathlib import Path
from importlib.metadata import version
import numpy as np
from thermal_bose_class import thermal_background,write_psd,params,first_cross


def main():
    from classy import Class
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    b=thermal_background();ks=np.geomspace(.002,1.,80)
    def solve(d):
        c=Class()
        try:
            c.set(d);c.compute()
            def val(k):
                v=getattr(c,k);return float(v() if callable(v) else v)
            return val('Omega_nu')*val('h')**2,val('Neff'),[c.pk_lin(float(k),5.) for k in ks]
        finally:
            c.struct_cleanup();c.empty()
    with tempfile.TemporaryDirectory() as d:
        file=Path(d)/'thermal.dat';write_psd(file)
        warm=params(b,file);warm['P_k_max_h/Mpc']=2.
        om,ne,wp=solve(warm)
        cold=params(b,cold_omega=om+1e-8);cold['P_k_max_h/Mpc']=2.
        _,_,cp=solve(cold)
    ratio=np.array(wp)/np.array(cp)
    valid=np.isfinite(ratio)&(ratio>0)
    T=np.full_like(ratio,np.nan);T[valid]=np.sqrt(ratio[valid]);cross=first_cross(ks,T)
    ref=json.loads((Path(__file__).parent/'results/CLASS-thermal-summary.json').read_text())
    ref_cross=ref['half_amplitude']['z5_k']
    result=dict(classy_community=version('classy-community'),omega=om,Neff=ne,analytic=b,
        z5_half_amplitude_k=cross,archived_CLASS_3_4_3_cross=ref_cross,
        relative_cross_shift=abs(cross/ref_cross-1),k_per_Mpc=ks.tolist(),amplitude=[float(v) if np.isfinite(v) else None for v in T],
        invalid_power_modes=int(np.count_nonzero(~valid)),invalid_k_per_Mpc=ks[~valid].tolist())
    Path(a.output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('k_per_Mpc','amplitude')}),flush=True)
    assert abs(om/b['omega_h2']-1)<.015
    assert abs(ne-3.044-b['delta_Neff_early'])<.003
    assert result['relative_cross_shift']<.03

if __name__=='__main__':main()
