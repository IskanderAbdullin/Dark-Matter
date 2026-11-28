"""Prepare an occupation ladder using checked near-packed CLASS configurations."""
import argparse
import json
from pathlib import Path


def read(path):
    return dict((k.strip(),v.strip()) for k,v in (
        line.split('=',1) for line in path.read_text().splitlines()))


def build(source, out, caps, tier):
    out=out.resolve()
    out.mkdir(parents=True,exist_ok=True)
    f=source/f'{tier}-packed.dat'
    (out/f'{tier}-packed.dat').write_bytes(f.read_bytes())
    conf=read(source/f'{tier}-F1.ini')
    st=conf.copy()
    st['ncdm_psd_filenames']=str(out/f'{tier}-packed.dat')
    rows=[]
    for cap in caps:
        tag=f'{tier}-F{cap:g}'
        st.update(T_ncdm=f'{float(conf["T_ncdm"])*cap**(-1/3):.17g}',
                  root=str(out/(tag+'_')))
        (out/(tag+'.ini')).write_text(''.join(f'{k} = {v}\n' for k,v in st.items()))
        rows.append(dict(label=tag,cap=cap,tier=tier,
                         T_ncdm_ratio=float(st['T_ncdm'])))
    return dict(cases=rows,scope="Fixed selected half-abundance and logistic-edge distribution family; "
                "no microphysical model, extremal linear-power proof, or observational likelihood.")


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('tier',choices=('compact','dense'))
    p.add_argument('caps',type=float,nargs='+')
    a=p.parse_args()
    assert all(f>0 for f in a.caps) and len(set(a.caps))==len(a.caps)
    print(json.dumps(build(a.source,a.output,a.caps,a.tier),indent=2))
