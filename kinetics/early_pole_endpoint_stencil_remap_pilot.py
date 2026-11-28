"""Cross the archived frozen states and native decay meshes without re-evolving."""
import json
from math import fsum

import numpy as np

from early_pole_decay_mode_clock_reference import X
from early_pole_massive_daughter_energy_stencil_reference import MassiveFermion
from finite_mass_qed_timelike_pole_batched_absorption_reference import load


SPECS={'fine':(65,90,72),'third':(95,125,108)}
LABELS=('scalar','fermion','parent')


def parts(d):
    sc=d.d.ed
    return [sc,d.d.fed,d.pe]


def archives(a):
    return (a['scalar_momentum_edges'],a['endpoint_fermion_momentum_edges'],a['endpoint_parent_momentum_edges'])


def remap(old,occ,new):
    # Work in the exact number-measure coordinate p^3 so a nonnegative
    # histogram stays nonnegative and any Pauli cap remains preserved.
    a=np.asarray(old)**3
    b=np.asarray(new)**3
    if a[0]!=b[0] or a[-1]!=b[-1]:
        raise ArithmeticError('Cannot number-remap different supports')
    W=np.maximum(0.,np.minimum(a[1:,None],b[None,1:])-np.maximum(a[:-1,None],b[None,:-1]))
    before=fsum(np.diff(a)*occ)
    q=np.asarray(occ)@W/np.diff(b)
    after=fsum(q*np.diff(b))
    if abs(after-before)>3e-15*max(abs(before),1e-20):
        raise ArithmeticError('Remap failed the number integral')
    return q,(after-before)/max(abs(before),1e-30)


def run():
    refs={g:load(g)[0] for g in SPECS}
    if refs['fine']['x_end']!=refs['third']['x_end']:
        raise ArithmeticError('Cannot compare the absorption clock at different scale factors')
    x=refs['fine']['x_end']
    t=np.log(x/X)
    out=[]
    for dest,shape in SPECS.items():
        r=refs[dest]
        d=MassiveFermion(*shape,r['mass_ratio'],r['elastic_rate'],True)
        if not all(np.array_equal(a,b) for a,b in zip(parts(d),archives(r),strict=True)):
            raise ArithmeticError(('Native archived axes do not match driver',dest))
        for source,v in refs.items():
            vals=[np.asarray(v['endpoint_'+label+'_occupations']) for label in LABELS]
            if source==dest:
                offers=vals;errs=[0.]*3
            else:
                offers,errs=zip(*(remap(o,f,a) for o,f,a in zip(archives(v),vals,parts(d),strict=True)),strict=True)
            if any(np.min(z)<0 for z in offers) or max(np.max(offers[1]),np.max(offers[2]))>1+1e-14:
                raise ArithmeticError('The occupation bounds were lost')
            y=np.concatenate(offers)
            n=fsum(d.d.w*y[:d.ns])
            change=d.w*d.rhs(t,y)
            s,fb,p=np.split(change,(d.ns,d.iP))
            vnet=fsum(s)
            pair=max(abs(vnet-fsum(fb)),abs(vnet+fsum(p)))
            if pair>3e-12*max(abs(vnet),1e-15):
                raise ArithmeticError('Remapped frozen decay no longer closes number')
            M=d.d.gap(x,y[:d.ns])
            if source==dest:
                raw=v['checkpoints'][-1]
                if abs(vnet/n/raw['net_log_growth_rate']-1)>3e-12 or abs(M/raw['scalar_pole_over_T']-1)>3e-13:
                    raise ArithmeticError('Could not reconstruct the saved endpoint')
            lo,hi=np.asarray(parts(d)[0][:-1]),np.asarray(parts(d)[0][1:])
            fraction_below=fsum(s[hi<=1e-5])/vnet
            if (not 1e-5 in parts(d)[0]) or not .9<fraction_below<1.01:
                raise ArithmeticError('Cold source localization not resolved')
            ell=.75*(hi**4-lo**4)/(hi**3-lo**3)
            out.append(dict(mesh=dest,input_histogram=source,q_faces=d.d.ed.tolist(),
                            net_scalar_by_q=s.tolist(),scalar_N=n,scalar_mass=M,
                            number_conservation=pair,remap_number_defects=errs,
                            normalized_net_log_clock=vnet/n,
                            flux_fraction_below_q_1e_5=fraction_below,
                            volume_weighted_mean_flux_q=float(s@ell/vnet),
                            norm_positive_negative_flux=[fsum(s[s>0])/n,-fsum(s[s<0])/n]))
    return out


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
