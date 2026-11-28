"""Draw the count CDF and resolved contrast for exact S2 and S4 state files."""
import argparse
import json
from math import fsum
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,load


def histogram(p,s):
    card=read(p)
    g=card['grid']
    if g!='fine' or (card['charge_div'],card['scalar_div'],card['x_end'])!=(2,s,.05):
        raise ArithmeticError('Expected the exact fixed-charge scalar-refined endpoint')
    d=ScalarSubshellChart(load(g)[0],SPECS[g],2,s,card['cut'])
    y=np.array(card['endpoint_occupation'])
    if y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1:
        raise ArithmeticError('Invalid saved phase space')
    counts=d.ds.w*y[:d.ns]
    counts/=np.sum(counts)
    return d.ds.ed,np.r_[0.,np.cumsum(counts)]


def render(paths,output,receipt=None):
    both=[histogram(p,s) for p,s in zip(paths,(2,4),strict=True)]
    for e,c in both:
        if not np.all(np.diff(c)>=0) or abs(c[-1]-1)>1e-14:
            raise ArithmeticError('Invalid CDF')
    union=np.union1d(both[0][0],both[1][0])
    patches=[np.interp(union**3,e**3,c) for e,c in both]
    delta=patches[1]-patches[0]
    j=np.argmax(abs(delta))
    if receipt is not None:
        record=dict(x=.05,chart='fine charge d2, scalar S2 and S4',
                    union_shells=len(union)-1,union_cdf_total_variation=fsum(abs(np.diff(delta)))/2,
                    maximum_absolute_cdf_difference=float(abs(delta[j])),
                    q_at_maximum=float(union[j]),
                    scope='Exact given the piecewise-constant occupation reconstruction in each chart.')
        receipt.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    q=np.union1d(np.linspace(2e-6,5.8e-6,480),union)
    q=q[(q>=2e-6)&(q<=5.8e-6)]
    cdf=[np.interp(q**3,e**3,c) for e,c in both]
    plt.rcParams.update({'svg.hashsalt':'early-pole-cdf-s2-s4','svg.fonttype':'none','font.size':11,
                         'axes.spines.right':False,'axes.spines.top':False})
    fig,(ax,bx)=plt.subplots(2,1,sharex=True,figsize=(7.8,5.6),
                             gridspec_kw={'height_ratios':(3,1.5)},layout='constrained')
    colors=['#176475','#b75131']
    for y,col,label in zip(cdf,colors,('S2: 93 scalar shells','S4: 145 scalar shells'),strict=True):
        ax.plot(q*1e6,y*100,color=col,label=label,linewidth=2.2)
    ax.set_ylim(0,101)
    ax.set_ylabel('Cumulative scalar count (%)')
    ax.legend(loc='upper left',frameon=False)
    ax.set_title('Cold scalar distribution at x = 0.05 (same charge grid)', loc='left',weight='semibold',pad=9)
    diff=(cdf[1]-cdf[0])*100
    bx.plot(q*1e6,diff,color='#4c4565',linewidth=2.1)
    bx.fill_between(q*1e6,diff,0,color='#4c4565',alpha=.11)
    bx.axhline(0,color='#777579',linewidth=.8)
    bx.set_ylabel('S4 − S2\n(percentage points)')
    bx.set_xlabel('Dimensionless scalar momentum q (×10⁻⁶)')
    bx.set_xlim(2,5.8)
    for a in (ax,bx):
        a.grid(alpha=.15)
    fig.savefig(output,metadata={'Date':None})
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('states',type=Path,nargs=2)
    p.add_argument('output',type=Path)
    p.add_argument('--receipt',type=Path)
    a=p.parse_args()
    render(a.states,a.output,a.receipt)
