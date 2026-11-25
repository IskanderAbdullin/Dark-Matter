"""Conditional linear transfer of actual frozen external-bath kinetic PSD.
No rescaling to target abundance; two numerical collisionless tracers, one scalar.
"""
import argparse,json,math,tempfile
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import PchipInterpolator, CubicSpline


def export_components(row,folder,tag='fine',split_multiple=10.,points=6000):
    raw=row['frozen_cell_spectrum'];edges=np.array(raw['edges']);f=np.array(raw['f_phi'])
    l,h=edges[:-1],edges[1:];cent=np.sqrt(l*h)
    assert np.min(f)>=-1e-7
    f=np.maximum(f,0.)
    N=np.sum(f*(h**3-l**3)/(6*np.pi**2))
    E=np.sum(f*(h**4-l**4)/(8*np.pi**2))
    Q2=np.sum(f*(h**5-l**5)/(10*np.pi**2))
    assert abs(N/row['final']['number']-1)<1e-8
    scale=float(cent[np.argmax(f*cent**3)])
    split=split_multiple*scale
    # Reconstruct log occupation at log-momentum cell midpoints. The raw
    # number and energy are independently checked below; never normalize.
    interp=PchipInterpolator(np.log(cent),np.log(np.maximum(f,1e-100)),extrapolate=False)
    def evaluate(q):
        out=np.full_like(q,1e-100)
        within=(q>=cent[0])&(q<=cent[-1])
        out[within]=np.exp(interp(np.log(q[within])))
        above=q>cent[-1]
        slope=min(-.01,(np.log(max(f[-1],1e-100))-np.log(max(f[-2],1e-100)))/(cent[-1]-cent[-2]))
        out[above]=np.maximum(f[-1],1e-100)*np.exp(slope*(q[above]-cent[-1]))
        return out
    components=[]
    for is_cold,s,limit in [(True,scale,60.),(False,1.,45.)]:
        u=np.unique(np.r_[np.geomspace(1e-7,.1,points//3),np.geomspace(.1,limit,2*points//3)])
        q=s*u;v=evaluate(q);exponent=(q/split)**8
        w=np.exp(-exponent) if is_cold else -np.expm1(-exponent)
        occupation=np.maximum(v*w,1e-100*np.exp(-u))
        psd=occupation/(2*np.pi)**3
        file=Path(folder)/(tag+('_cold.dat' if is_cold else '_hot.dat'))
        np.savetxt(file,np.column_stack((u,psd)),fmt='%.17e')
        # Local spline positivity stress check. A tiny tail overshoot is a
        # diagnostic, not a proof of CLASS's own internal interpolation.
        probe=(u[:-1]+u[1:])/2
        cv=CubicSpline(u,psd)(probe)
        neg_mass=4*np.pi*s**3*trapezoid(probe**2*np.maximum(-cv,0.),probe)
        moments=[float(4*np.pi*s**(k+3)*trapezoid(u**(k+2)*psd,u)) for k in (0,1,2)]
        components.append(dict(file=str(file),scale=s,points=len(u),number=moments[0],energy=moments[1],q2=moments[2],spline_negative_number=neg_mass))
    sums=[sum(c[k] for c in components) for k in ('number','energy','q2')]
    errors=[float(a/b-1) for a,b in zip(sums,[N,E,Q2])]
    assert max(map(abs,errors))<.015,errors
    assert sum(c['spline_negative_number'] for c in components)<N*1e-9
    return dict(raw_number=float(N),raw_energy=float(E),raw_q2=float(Q2),
                split=split,split_multiple=split_multiple,components=components,moment_relative_errors=errors,
                reconstructed_hot_number_fraction=components[1]['number']/sums[0],
                method='positive log-PCHIP with smooth partition; no abundance rescaling')


def transformed_quadrature_check(export,bins):
    # CLASS strategy 2: q=1/t-1 with equally spaced t in (0,1).
    errors=[]
    for c in export['components']:
        u,f=np.loadtxt(c['file']).T
        t=np.arange(1,bins+1)/(bins+1);x=1/t-1
        values=np.zeros_like(x);inside=(x>=u[0])&(x<=u[-1])
        values[inside]=CubicSpline(u,f)(x[inside]);values[x<u[0]]=f[0]
        assert f[-1]>0 and f[-2]>=f[-1],('invalid CLASS PSD extrapolation',c.get('label'),f[-2:])
        values[x>u[-1]]=f[-1]*np.exp((x[x>u[-1]]-u[-1])*(f[-1]-f[-2])/f[-1]/(u[-1]-u[-2]))
        # CLASS computes q/f0 * df0/dq; an invisible subnormal tail can
        # overflow the reciprocal even when all physical moments converge.
        assert np.all(np.isfinite(values)),('nonfinite PSD',c.get('label'))
        nz=values!=0
        with np.errstate(over='ignore',divide='ignore',invalid='ignore'):
            assert np.all(np.isfinite(x[nz]/values[nz])),('CLASS dlnf reciprocal overflow',c.get('label'),bins)
        w=values/(bins+1)/t**2
        errors.append([float(4*np.pi*c['scale']**(k+3)*(w@x**(k+2))/c[key]-1)
                       for k,key in enumerate(['number','energy','q2'])])
    assert np.max(np.abs(errors))<.001,errors
    return errors


def transfer(row,export,ks,zs,bins=120,precision=None,kmax_h=45.):
    from classy import Class
    base={'h':.674,'omega_b':.02237,'T_cmb':2.7255,'N_ur':3.044+row['final']['delta_Neff_chi2_if_stable_free_streaming'],
          'A_s':2.1e-9,'n_s':.965,'tau_reio':.0544,'gauge':'synchronous','output':'mPk',
          'P_k_max_h/Mpc':kmax_h,'z_max_pk':5.1}
    def solve(p):
        c=Class()
        try:
            c.set(p);c.compute()
            def scalar(name):
                v=getattr(c,name)
                return float(v() if callable(v) else v)
            h=scalar('h')
            bg=c.get_background()
            keys=[k for k in bg if k.startswith('(.)rho_ncdm[') or k.startswith('(.)rho_ncdm__')]
            keys.sort(key=lambda k:int(k.split('[')[-1].rstrip(']')) if '[' in k else int(k.split('__')[-1]))
            if not keys and 'N_ncdm' in p:
                raise RuntimeError(('Unknown NCDM background keys',list(bg)))
            species=[float(bg[k][-1]) for k in keys]
            return {'omega_ncdm_h2':scalar('Omega_nu')*h*h,'Neff':scalar('Neff'),'rho_species_z0':species,
                    'power':{str(z):[float(c.pk_lin(float(k),float(z))) for k in ks] for z in zs}}
        finally:
            try:c.struct_cleanup();c.empty()
            except Exception:pass
    comps=export['components'];tr=(3.909/106.75)**(1/3)
    nc=len(comps)
    assert nc>=1
    p=dict(base,N_ncdm=nc,m_ncdm=','.join(['86']*nc),T_ncdm=','.join(str(c['scale']*tr) for c in comps),
           deg_ncdm=','.join(['1']*nc),use_ncdm_psd_files=','.join(['1']*nc),
           ncdm_psd_filenames=','.join(c['file'] for c in comps),omega_cdm=1e-8)
    # CLASSpp 26 uses distinct legacy-to-instance aliases for standard NCDM.
    # Its old CLASS input names otherwise silently trigger automatic sampling.
    from importlib.metadata import version,PackageNotFoundError
    try:
        community=version('classy-community')
    except PackageNotFoundError:
        community=None
    if community and int(community.split('.')[0])>=26:
        p.update(quadrature_strategy_ncdm_standard=','.join(['2']*nc),
                 N_momentum_bins_ncdm_standard=','.join([str(bins)]*nc),
                 N_momentum_bins_bg_ncdm_standard=','.join([str(bins)]*nc))
    else:
        p.update(ncdm_quadrature_strategy=','.join(['2']*nc),ncdm_N_momentum_bins=','.join([str(bins)]*nc))
    if precision:
        precision=dict(precision)
        if community and int(community.split(".")[0])>=26 and "tol_perturbations_integration" in precision:
            precision["tol_perturb_integration"]=precision.pop("tol_perturbations_integration")
        p.update(precision)
    quadrature_errors=transformed_quadrature_check(export,bins)
    warm=solve(p)
    target=row['final']['omega_phi_h2_m86_g106p75']
    assert abs(warm['omega_ncdm_h2']/target-1)<.025,(warm['omega_ncdm_h2'],target)
    assert len(warm['rho_species_z0'])==nc,warm['rho_species_z0']
    fractions=np.asarray(warm['rho_species_z0'])/sum(warm['rho_species_z0'])
    expected=np.asarray([c['number'] for c in comps])/sum(c['number'] for c in comps)
    assert np.max(np.abs(fractions-expected))<.015,(fractions,expected)
    if nc>=3:
        assert np.max(np.abs(fractions/expected-1))<.02,(fractions,expected)
    ratio=float(fractions[1:].sum())
    assert abs(ratio-export['reconstructed_hot_number_fraction'])<.015,(ratio,export['reconstructed_hot_number_fraction'])
    delta_phi=export['raw_energy']/(7*np.pi**2/120)*tr**4/(4/11)**(4/3)
    assert abs(warm['Neff']-base['N_ur']-delta_phi)<.003,(warm['Neff'],base['N_ur'],delta_phi)
    cold=solve(dict(base,omega_cdm=warm['omega_ncdm_h2']+1e-8))
    result={}
    for z in zs:
        key=str(z);r=np.asarray(warm['power'][key])/np.asarray(cold['power'][key]);assert np.all(np.isfinite(r)) and np.all(r>0)
        # Equal-density adiabatic models must approach the same large-scale power.
        # The previous Newtonian ultra-cold-tracer run failed this condition.
        large=ks<=.03
        assert not np.any(large) or np.max(np.abs(r[large]-1))<.05,(z,r[large].tolist())
        crossings=np.flatnonzero((r[:-1]>=.5)&(r[1:]<.5))
        half=None
        if len(crossings):
            i=crossings[0];half=float(np.exp(np.interp(.5,[r[i+1],r[i]],[np.log(ks[i+1]),np.log(ks[i])])))
        result[key]={'power_ratio':r.tolist(),'half_power_k_per_Mpc':half,
                     'P_ratio_k1':float(np.interp(0.,np.log(ks),r)) if ks[0]<=1.<=ks[-1] else None,
                     'P_ratio_k10':float(np.interp(np.log(10),np.log(ks),r)) if ks[0]<=10.<=ks[-1] else None}
    return {'warm_background':{k:v for k,v in warm.items() if k!='power'},'expected_omega':target,
            'delta_phi_early':float(delta_phi),'chi2_Nur':row['final']['delta_Neff_chi2_if_stable_free_streaming'],
            'hot_tracer_mass_fraction':ratio,'tracer_mass_fractions':fractions.tolist(),'ratios':result,'quadrature_strategy':2,
            'momentum_bins':bins,'classy_community_version':community,'gauge':'synchronous','quadrature_moment_relative_errors':quadrature_errors}


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',nargs='+',required=True);p.add_argument('--output',required=True)
    p.add_argument('--analytic-only',action='store_true');a=p.parse_args()
    ks=np.geomspace(.003,20.,125);zs=[0.,5.];rows=[]
    with tempfile.TemporaryDirectory() as folder:
        for i,path in enumerate(a.input):
            row=json.loads(Path(path).read_text())[0]
            export=export_components(row,folder,tag=str(i),points=6000)
            alt=export_components(row,folder,tag=str(i)+'_alternate',split_multiple=8.,points=8000)
            result={'input':path,'n':row['n'],'x0':row['x0'],'kappa':row['kappa'],'export':export,'alternate_export':alt,
                    'quadrature_120':transformed_quadrature_check(export,120),'quadrature_240':transformed_quadrature_check(alt,240)}
            if not a.analytic_only:
                result['CLASS']=transfer(row,export,ks,zs)
                result['alternate_CLASS']=transfer(row,alt,ks,zs,240)
                for z in zs:
                    k=str(z);r=np.asarray(result['CLASS']['ratios'][k]['power_ratio']);r2=np.asarray(result['alternate_CLASS']['ratios'][k]['power_ratio'])
                    pre=ks<=min(10.,result['CLASS']['ratios'][k]['half_power_k_per_Mpc'] or 10.)
                    result.setdefault('partition_and_sampling_comparison',{})[k]={'max_relative_difference_before_half':float(np.max(np.abs(r2[pre]/r[pre]-1)))}
            # File paths are temporary; retain only numerical metadata.
            for ex in (export,alt):
                for c in ex['components']:c.pop('file')
            rows.append(result)
            print(json.dumps({k:v for k,v in result.items() if k not in ('CLASS','alternate_CLASS')}),flush=True)
    Path(a.output).write_text(json.dumps({'k_per_Mpc':ks.tolist(),'results':rows,'scope':'Conditional external-bath freeze-out + collisionless linear perturbations; no observation likelihood'},indent=2)+'\n')

if __name__=='__main__':main()
