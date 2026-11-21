"""Direct SPARC table QC. Run: python sparc_audit.py MassModels_Lelli2016c.mrt.
Requires Python stdlib only; produces descriptive counts, not a model fit.
"""
import sys,collections,json
from pathlib import Path
def audit(filename):
    rows=[]
    for line in Path(filename).read_text(encoding="utf-8").splitlines():
        a=line.strip().split()
        if len(a)!=10 or not a[0][:1].isalpha():continue
        try: n=list(map(float,a[1:]))
        except ValueError:continue
        rows.append(dict(id=a[0],D=n[0],r=n[1],v=n[2],dv=n[3],
                         vg=n[4],vd=n[5],vb=n[6]))
    by=collections.defaultdict(list)
    for x in rows: by[x["id"]].append(x)
    def dm2(x):return x["v"]**2-x["vg"]*abs(x["vg"])-.5*x["vd"]**2-.7*x["vb"]**2
    counts=sorted(len(g) for g in by.values())
    by_galaxy=[]
    for gid,rs in sorted(by.items()):
        f=sorted(dm2(x)/x["v"]**2 for x in rs)
        by_galaxy.append(dict(id=gid,n=len(rs),min_r=min(x["r"] for x in rs),
           max_r=max(x["r"] for x in rs),negative_gas=sum(x["vg"]<0 for x in rs),
           negative_dm2=sum(dm2(x)<0 for x in rs),
           median_dm2_fraction=f[len(f)//2]))
    summary=dict(rows=len(rows),galaxies=len(by),points_min=counts[0],
       points_median=counts[len(counts)//2],points_max=counts[-1],
       invalid=sum(x["dv"]<=0 or x["r"]<=0 or x["v"]<0 for x in rows),
       negative_gas=sum(x["vg"]<0 for x in rows),
       negative_dm2=sum(dm2(x)<0 for x in rows),
       galaxies_negative_dm2=sum(p["negative_dm2"]>0 for p in by_galaxy),
       galaxies_median_dm2_fraction_gt_07=sum(p["median_dm2_fraction"]>.7 for p in by_galaxy))
    assert summary["rows"]==3391 and summary["galaxies"]==175
    assert summary["negative_gas"]==361
    return dict(summary=summary,by_galaxy=by_galaxy)
if __name__=="__main__":print(json.dumps(audit(sys.argv[1]),indent=2))
