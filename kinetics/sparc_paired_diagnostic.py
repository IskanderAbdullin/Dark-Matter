"""Paired descriptive galaxy bootstrap; no formal likelihood claim."""
import argparse,json
from pathlib import Path
import numpy as np

def paired(a,b,rng):
    keys=sorted(set(a['rows'])&set(b['rows']))
    assert len(keys)==len(a['rows'])==len(b['rows'])
    dd=np.array([a['rows'][k]['test']-b['rows'][k]['test'] for k in keys])
    n=np.array([a['rows'][k]['num_test'] for k in keys])
    assert np.all(n==np.array([b['rows'][k]['num_test'] for k in keys]))
    ix=rng.integers(0,len(keys),size=(10000,len(keys)))
    boot=dd[ix].sum(axis=1)/n[ix].sum(axis=1)
    return dict(A=[a['kind'],a['mass_eV']],B=[b['kind'],b['mass_eV']],galaxies=len(keys),
        total_test_delta=float(dd.sum()),per_point_delta=float(dd.sum()/n.sum()),median_galaxy_per_point_delta=float(np.median(dd/n)),
        A_wins=int(sum(dd<0)),B_wins=int(sum(dd>0)),tie=int(sum(dd==0)),
        paired_galaxy_bootstrap_95pct=np.quantile(boot,[.025,.975]).tolist(),
        largest_abs_contributions=[dict(galaxy=keys[i],delta_chi2=float(dd[i])) for i in np.argsort(abs(dd))[-5:][::-1]])

def run(data):
    cs=data['cases'];rng=np.random.default_rng(20261125)
    best={k:min([x for x in cs if x['kind']==k],key=lambda z:z['objective']) for k in set(x['kind'] for x in cs)}
    a86=next(x for x in cs if x['kind']=='one' and x['mass_eV']==86)
    b86=next(x for x in cs if x['kind']=='two' and x['mass_eV']==86)
    b74=next(x for x in cs if x['kind']=='two' and x['mass_eV']==74)
    prs=[(b86,a86),(b74,a86),(best['two'],best['one']),(best['one'],best['pis']),(best['two'],best['pis']),(best['burkert'],best['pis'])]
    if 'classical' in best:prs.extend([(best['one'],best['classical']),(best['two'],best['classical']),(best['classical'],best['pis'])])
    return dict(sample=data['sample'],prior=data['prior'],best_training=[{k:v for k,v in x.items() if k!='rows'} for x in best.values()],
        paired=[paired(a,b,rng) for a,b in prs],
        scope='Descriptive interleaved radial holdout; galaxy resampling keeps within-galaxy points paired. Not independent data, calibrated confidence intervals, Bayes factors, or population mass inference.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--quality',required=True);p.add_argument('--pilot',required=True);p.add_argument('--prior',required=True);p.add_argument('--classical-quality');p.add_argument('--classical-prior');p.add_argument('--output',required=True);a=p.parse_args()
    q=json.loads(Path(a.quality).read_text());pr=json.loads(Path(a.prior).read_text());pi=json.loads(Path(a.pilot).read_text())
    for dat,extra in [(q,a.classical_quality),(pr,a.classical_prior)]:
        if extra:dat['cases']+=json.loads(Path(extra).read_text())['cases']
    v=[run(x) for x in [q,pi,pr]]
    Path(a.output).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    for x in v:
        print(x['sample'],x['prior'],flush=True)
        for z in x['paired']:print(z['A'],z['B'],z['total_test_delta'],z['A_wins'],z['B_wins'],z['paired_galaxy_bootstrap_95pct'],flush=True)
