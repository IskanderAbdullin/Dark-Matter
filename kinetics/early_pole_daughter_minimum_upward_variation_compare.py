"""Cross-check finite upward-variation witnesses and conservative numerical lower bounds."""
import argparse
import json
from pathlib import Path


def compare(raw,coarse_duals,refined_duals):
    if len(raw)!=6 or len(coarse_duals)!=2 or len(refined_duals)!=2:
        raise ValueError('Missing archive fixture or moment controls.')
    out=[]
    for a,b,c,d,e,f,u,w,v,z in zip(*raw,*coarse_duals,*refined_duals,strict=True):
        allv=(a,b,c,d,e,f,u,w,v,z)
        key=(a['hard_depth'],a['cold_modes'])
        if any((g['hard_depth'],g['cold_modes'])!=key for g in allv):
            raise ValueError('Different physical fixtures.')
        for x in (a,b,c,d,e,f):
            p=x['positive_witness']
            if p['verified_paired_rate']<=0. or p['total_upward_variation']<p['numerical_dual_lower']-1e-8:
                raise ValueError('Invalid positive witness.')
        terms=[g['numerical_lower_total_upward_variation_for_nonnegative_pair'] for g in (u,w,v,z)]
        if any(p is None or p<=0. for p in terms):
            raise ValueError('The monotone certificate did not give a positive floor.')
        if terms[2]>f['positive_witness']['total_upward_variation'] or terms[3]>b['positive_witness']['total_upward_variation']:
            raise ValueError('Unexpected conflict with a finite verified upper witness.')
        out.append(dict(hard_depth=key[0],cold_modes=key[1],three_moment_controls=[a,b,c,d],
                        two_moment_controls=[e,f],two_moment_numerical_continuum_lower=terms[2],
                        three_moment_numerical_continuum_lower=terms[3]))
    return dict(scope='Total upward variation is a sum of increases, not relative occupation. Numerical lower bounds are conditional on the retained geometry.',
                fixtures=out,dual_cutoff_certificates=coarse_duals+refined_duals)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--witnesses',nargs=6,required=True)
    p.add_argument('--coarse-duals',nargs=2,required=True)
    p.add_argument('--refined-duals',nargs=2,required=True)
    a=p.parse_args()
    print(json.dumps(compare(*([json.loads(Path(f).read_text()) for f in getattr(a,k)]
                                for k in ('witnesses','coarse_duals','refined_duals'))),
                     indent=2,allow_nan=False))
