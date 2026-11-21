"""SPARC Table 1 parsing and quality audit. Requires Python standard library only.
python sparc_quality.py SPARC_Lelli2016c.mrt
"""
import json,sys
from pathlib import Path
def parse(path):
    data={}
    for line in Path(path).read_text().splitlines():
        a=line.split()
        if len(a)!=19 or not a[0][0].isalpha():continue
        try:
            q=int(a[17]);inc=float(a[5]);D=float(a[2])
            if not (1<=q<=3 and inc>0 and D>0):continue
            data[a[0]]=dict(quality=q,inclination_deg=inc,
                inclination_error_deg=float(a[6]),distance_Mpc=D,
                distance_error_Mpc=float(a[3]),ref=a[18])
        except(ValueError,IndexError):continue
    counts={q:sum(x['quality']==q for x in data.values()) for q in (1,2,3)}
    s=dict(galaxies=len(data),quality_counts=counts,
           low_inclination=sum(x['inclination_deg']<30 for x in data.values()),
           quality_le2_inc_ge30=sum(x['quality']<=2 and x['inclination_deg']>=30
                   for x in data.values()))
    assert s['galaxies']==175 and counts=={1:99,2:64,3:12}
    return dict(summary=s,per_galaxy=data)
if __name__=='__main__':print(json.dumps(parse(sys.argv[1]),indent=2))
