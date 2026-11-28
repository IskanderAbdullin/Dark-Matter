"""Resolve both explicit-population replacement orders after direct edge reintegration."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_normalized_profile_envelope import fit_state,geometry
from early_pole_endpoint_cold_shared_reference_scan import fit_refined_left
from early_pole_endpoint_cold_spatial_population_factor_bounds import FIELDS,field


def decompose(raw):
    if raw.ndim!=3 or raw.shape[:2]!=(2,2) or np.any(raw<0.) or not np.all(np.isfinite(raw)):
        raise ArithmeticError('Expected a positive two-axis rectangle')
    total=np.array([[fsum(row) for row in x] for x in raw])
    if np.min(total)<=0.:
        raise ArithmeticError('Nonpositive integrated rectangle entry')
    p=raw/total[:,:,None]
    scalar=(p[0,1]-p[0,0],p[1,1]-p[1,0])
    kernel=(p[1,0]-p[0,0],p[1,1]-p[0,1])
    native=p[1,1]-p[0,0]
    u=(scalar[0]+scalar[1])/2
    v=(kernel[0]+kernel[1])/2
    log=np.log(total)
    lpop=((log[0,1]-log[0,0])+(log[1,1]-log[1,0]))/2
    lk=((log[1,0]-log[0,0])+(log[1,1]-log[0,1]))/2
    L=log[1,1]-log[0,0]
    def tv(x):
        return fsum(abs(x))/2
    err=max(abs(lpop+lk-L),float(np.max(abs(u+v-native))))
    if err>2e-14 or abs(tv(native)-tv(u))>tv(v)+2e-14:
        raise ArithmeticError('Independent two-order decomposition failed')
    return dict(response_rectangle=total.tolist(),normalized_profile_rectangle=p.tolist(),
                native_log_ratio=float(L),explicit_scalar_log=float(lpop),bare_kernel_log=float(lk),
                explicit_scalar_share_of_log=float(lpop/L) if L else None,
                four_corner_log_interaction=float(log[1,1]-log[1,0]-log[0,1]+log[0,0]),
                native_profile_TV=tv(native),both_same_kernel_population_TVs=[tv(z) for z in scalar],
                both_same_population_bare_kernel_TVs=[tv(z) for z in kernel],
                symmetric_scalar_profile_contrast=u.tolist(),symmetric_kernel_profile_contrast=v.tolist(),
                symmetric_scalar_profile_TV=tv(u),symmetric_kernel_profile_TV=tv(v),
                normalized_profile_interaction_TV=tv(scalar[1]-scalar[0]),
                maximum_additive_defect=err)


def neutrality(out,masses):
    def matrix(k):
        return np.asarray(out[k]['response_rectangle'])
    an,bn,ax,bx=(matrix(k) for k in ('forward_number','inverse_number','forward_excess','inverse_excess'))
    m=np.asarray(masses)[:,None]
    dn=an-bn
    dx=ax-bx
    if min(an.min(),bn.min(),dn.min(),ax.min(),bx.min(),dx.min(),m.min())<=0.:
        raise ArithmeticError('Cannot form positive conditional neutrality roots')
    num=bn/dn
    energy=(m*bn+bx)/(m*dn+dx)
    # Stable separation of the two tiny conditional roots, with its direct ratio check.
    split=(bx/bn-ax/an)/(1-bn/an)/(m+dx/dn)
    if np.max(abs(energy/num-1.-split))>8e-16:
        raise ArithmeticError('Centered and direct neutrality roots disagree')
    return dict(number_root=num.tolist(),energy_root=energy.tolist(),
                energy_root_over_number_root_minus_one=split.tolist())


def run(old,new,right,states,archive,levels,plan,fixture):
    if tuple(levels) not in ((2,4),(4,8),(8,16)) or plan['levels']!=levels:
        raise ArithmeticError('Need single validated adjacent archived spatial pair')
    n,m=levels
    a=fit_refined_left(new,old,states[0],archive,plan['targeted_parent_panel_splits'],level=n)
    b=fit_state(states[1],archive,m,right)
    ix=[fixture['levels'].index(i) for i in levels]
    if [a[-1],b[-1]]!=[fixture['states_sha256_uncompressed'][i] for i in ix] or [a[-1],b[-1]]!=plan['states_sha256_uncompressed']:
        raise ArithmeticError('The refinement, reference partition, and endpoints are not bound')
    edges=np.asarray(fixture['common_q_edges'])
    if [edges[0],edges[-1]]!=[1e-6,1.5e-5]:
        raise ArithmeticError('The original fixed display band changed')
    out=[]
    for t in (0,1):
        O,U,V=[geometry(y[t]) for y in (old,new,right)]
        if not all(np.all(np.isin(edges,e[0])) for e in (O,U,V)):
            raise ArithmeticError('Original common display edge is absent')
        row={}
        for name,i,j in (*FIELDS,('inverse_excess',2,1)):
            data=field((U,V),(a,b),edges,i,j)
            C=np.asarray(data['native_and_crossed_bin_response_bounds'])
            if data['crossed_panels_with_nonconstant_factor']!=[[0,0],[0,0]] or not np.array_equal(C[:,:,:,0],C[:,:,:,1]):
                raise ArithmeticError('At least one population substitution is unresolved')
            out_diag=field((O,V),(a,b),edges,i,j)
            old_diag=np.asarray(out_diag['native_and_crossed_bin_response_bounds'])[0,0,:,0]
            ratio=float(np.max(abs(C[0,0,:,0]/old_diag-1.)))
            if not ratio<6e-5:
                raise ArithmeticError('Direct split changed original coarse bin response')
            row[name]=dict(maximum_native_coarse_bin_relative_recomposition=ratio,**decompose(C[:,:,:,0]))
        out.append(row)
    defect=max(float(np.max(abs(np.asarray(out[0][k]['response_rectangle'])/v['response_rectangle']-1.))) for k,v in out[1].items())
    profile=max(float(np.max(abs(np.asarray(out[0][k]['normalized_profile_rectangle'])-v['normalized_profile_rectangle']))) for k,v in out[1].items())
    if defect>6e-5 or profile>2e-9:
        raise ArithmeticError('Full crossed rectangle did not close across both q meshes')
    split_a,split_b=[neutrality(u,[a[2],b[2]]) for u in out]
    split_defect=float(np.max(abs(np.asarray(split_a['energy_root_over_number_root_minus_one'])-
                                   split_b['energy_root_over_number_root_minus_one'])))
    if split_defect>1e-10:
        raise ArithmeticError('The centered neutrality separation did not resolve')
    return dict(levels=levels,states_sha256_uncompressed=[a[-1],b[-1]],common_q_edges=edges.tolist(),
                conditional_neutrality_rectangle=split_b,maximum_paired_neutral_split_defect=split_defect,
                observables=out[1],maximum_relative_paired_gross_defect=defect,maximum_paired_profile_defect=profile,
                scope='Fixed-state panel factor substitutions in both orders; remaining combined kernel is not a dynamical attribution.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ('left-baseline','left-refined','right','states'):
        parser.add_argument('--'+k,nargs=2,type=Path,required=True)
    for k in ('archive','plan','display-fixture'):
        parser.add_argument('--'+k,type=Path,required=True)
    parser.add_argument('--levels',nargs=2,type=int,required=True)
    a=parser.parse_args()
    print(json.dumps(run(a.left_baseline,a.left_refined,a.right,a.states,a.archive,a.levels,
                         json.loads(a.plan.read_text()),json.loads(a.display_fixture.read_text())),indent=2,allow_nan=False))
