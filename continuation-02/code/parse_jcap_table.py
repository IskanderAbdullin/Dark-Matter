"""Parse the 22-row JCAP table directly from uploaded source TeX.
Only published reduced chi^2; no original velocity-point likelihood reconstructed.
"""
import re,sys,json,statistics
from pathlib import Path
raw=Path(sys.argv[1]).read_bytes().decode('latin-1')
start=re.search(r'1\s*&\s*\\text\{F563-1\}',raw)
if start is None: raise ValueError('Missing original first row')
block=raw[start.start():].split(r'\end{tabular}')[0]
rows=[]
for r in block.split(r'\\'):
    if not re.match(r'^\s*\d+\s*&',r): continue
    r=re.sub(r'\\(?:text|textbf)\{([^{}]*)\}',r'\1',r)
    a=[v.replace('\\','').strip() for v in r.split('&')]
    assert len(a)==12,a
    rows.append(dict(galaxy=a[1],nu=float(a[4]),m=None if a[5]=='---' else float(a[5]),
       ch_min=float(a[6]),ch86=float(a[7]),chIS=float(a[8]),chPIS=float(a[9]),chBur=float(a[10]),chNFW=float(a[11])))
assert len(rows)==22,len(rows)
summary=dict(n=len(rows),nu_near_crit=sum(r['nu']>=2.6 for r in rows),
 nu_zero=sum(r['nu']==0 for r in rows),
 ch86_better_IS=sum(r['ch86']<r['chIS'] for r in rows),
 ch86_equal_IS=sum(r['ch86']==r['chIS'] for r in rows),
 ch86_better_Bur=sum(r['ch86']<r['chBur'] for r in rows),
 ch86_better_PIS=sum(r['ch86']<r['chPIS'] for r in rows),
 median_delta_IS=statistics.median(r['ch86']-r['chIS'] for r in rows))
print(json.dumps({'summary':summary,'rows':rows},ensure_ascii=False,indent=2))