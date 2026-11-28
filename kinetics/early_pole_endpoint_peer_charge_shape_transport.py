"""Compare scalar-distribution shape between different native grids at matched charge levels."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_endpoint_scalar_shape_transport import CUT,cdf,pair
from early_pole_endpoint_scalar_spectral_moment_budget import LABELS,integrate
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


def unpack(c,g,n):
    if (c['grid'],c['subcells_per_archived_cell'],c['x_end'],c['reference_x'],c['steps'])!=(g,n,.05,X,2048):
        raise ArithmeticError('Chart and clock provenance changed')
    d=CellChart(load(g)[0],SPECS[g],n,reference_x=X)
    a=np.asarray(c['endpoint_occupation'])
    if a.shape!=(d.n,) or c['counts']!=[d.ns,d.nf,d.np] or not np.all(np.isfinite(a)) or min(a)<0 or max(a[d.ns:])>1:
        raise ArithmeticError('Invalid endpoint occupations')
    return d,a,d.state(.05,a)[0]


def quantile_support_bounds(left,right):
    # The quantile in each probability span lies within its containing cell
    # for every positive distribution with the same archived cell counts.
    a,fa,*_=left
    b,fb,*_=right
    if a[0]!=b[0] or a[-1]!=b[-1]:
        raise ArithmeticError('Count-only cold comparison requires common exact terminal edges')
    z=np.unique(np.clip(np.r_[fa,fb,0.,1.],0.,1.))
    u=np.searchsorted(fa,z[:-1],side='right')-1
    v=np.searchsorted(fb,z[:-1],side='right')-1
    u=np.minimum(u,len(a)-2)
    v=np.minimum(v,len(b)-2)
    if np.any(u<0) or np.any(v<0):
        raise ArithmeticError('Invalid native probability tails')
    A,B,C,D=a[u],a[u+1],b[v],b[v+1]
    gap=np.maximum(0.,np.maximum(C-B,A-D))
    diam=np.maximum(abs(A-D),abs(B-C))
    dm=np.diff(z)
    witness=None
    if np.all(np.maximum(A,C)<np.minimum(B,D)):
        # Positive atoms inside the common support at each probability interval:
        # both sets of archived counts hold exactly, to floating summation error.
        q=(np.maximum(A,C)+np.minimum(B,D))/2
        if np.any(np.diff(q)<0):
            raise ArithmeticError('Overlapping quantile support did not give an ordered witness')
        residual=0.
        for ed,f in ((a,fa),(b,fb)):
            k=np.searchsorted(ed,q,side='right')-1
            projected=np.bincount(k,weights=dm,minlength=len(ed)-1)
            residual=max(residual,max(abs(projected-np.diff(f))))
        if residual>3e-14:
            raise ArithmeticError('Common positive count-only witness violates the archived shell counts')
        witness=dict(atoms=len(q),max_absolute_cell_mass_defect=residual)
    return fsum(dm*gap),fsum(dm*diam),witness


def run(fine,third):
    if len(fine)!=len(third) or len(fine) not in (3,4):
        raise ArithmeticError('Need both entire independent ladders of length 3 or 4')
    d=(1,2,4,8)[:len(fine)]
    mass=unpack(fine[0],'fine',1)[-1]
    out=[]
    for f,t,n in zip(fine,third,d,strict=True):
        (a,u,ma),(b,v,mb)=unpack(f,'fine',n),unpack(t,'third',n)
        if np.array_equal(a.ds.ed,b.ds.ed):
            raise ArithmeticError('Need independent scalar supports on the two chart families')
        moment=[]
        defect=[]
        supports=[]
        for h,z in ((a,u),(b,v)):
            lo,hi=[integrate(h,z,mass,j) for j in (24,48)]
            err=np.max(abs(lo-hi)/np.maximum(abs(hi),1e-105))
            if err>2e-12:
                raise ArithmeticError('Unresolved common-mass moment measure')
            defect.append(err)
            moment.append(np.array([fsum(x) for x in hi.T]))
            supports.append([h.ds.ed[0],h.ds.ed[-1]])
        P,Q=moment
        if min(P)<=0 or min(Q)<=0:
            raise ArithmeticError('Nonpositive reconstructed moment')
        z=pair(cdf(a,u,a.ds.ed[-1]),cdf(b,v,b.ds.ed[-1]))
        if CUT not in a.ds.ed or CUT not in b.ds.ed:
            raise ArithmeticError('Cold cutoff is not a shared archived edge')
        cold=cdf(a,u,CUT),cdf(b,v,CUT)
        c=pair(*cold)
        lo,hi,witness=quantile_support_bounds(*cold)
        (qf,pf,*_),(qt,pt,*_)=cold
        common=np.intersect1d(qf,qt)[1:-1]
        shared=[dict(q=float(q),third_minus_fine_strict_lower_mass=float(pt[np.searchsorted(qt,q)]-pf[np.searchsorted(qf,q)]))
                for q in common]
        tvmin=max((abs(v['third_minus_fine_strict_lower_mass']) for v in shared),default=0.)
        if not lo-2e-13*hi<=c['wasserstein1_q']<=hi*(1+2e-13) or c['total_variation']+2e-14<tvmin:
            raise ArithmeticError('Piecewise-constant numerical fit escaped its count-only transport interval')
        if np.max(abs(np.array(z['means'])-np.array([P[1]/P[0],Q[1]/Q[0]])))>3e-14*max(z['means']):
            raise ArithmeticError('Transport mean did not agree with independent quadrature')
        out.append(dict(charge_subcells=n,scalar_bins=[a.ns,b.ns],supports=supports,
                        conditional_cold=c,cold_count_only_wasserstein1_interval=[lo,hi],cold_common_count_only_atom_witness=witness,cold_exact_shared_edge_mass_differences=shared,
                        cold_count_only_tv_lower_bound=tvmin,full_distribution=z,third_relative_moment=(Q/P-1).tolist(),
                        third_relative_per_particle_moment=((Q[1:]/Q[0])/(P[1:]/P[0])-1).tolist(),
                        scalar_mass_third_relative=mb/ma-1,maximum_moment_check=max(defect)))
    return dict(x=.05,fixed_measurement_mass=mass,cold_support_upper=CUT,columns=LABELS,
                chart_order=('fine','third'),comparisons=out,
                scope='Only cross-grid numerical histogram transport and common-mass moments. Narrower number spread does not imply more closely matching spectra.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fine',type=Path,nargs='+',required=True)
    p.add_argument('--third',type=Path,nargs='+',required=True)
    a=p.parse_args()
    print(json.dumps(run([read(x) for x in a.fine],[read(x) for x in a.third]),indent=2,allow_nan=False))
