"""Positive resolved freeze-in distributions on a prescribed cooling history."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from finite_cooling_reheat_portal_reference import TF, SourceFold
from three_species_decay_onsager_reference import Es, M, P, m, mp


GL3=roots_legendre(3)


def bins(edge, mass, center, width, isotropic, quadrupole, red):
    """Exact partial-angle quadrature of the count and p² in each bin."""
    center, width = np.ravel(center), np.ravel(width)
    n, a = np.ravel(isotropic), np.ravel(quadrupole)
    total = np.zeros((2, len(edge)-1))
    exc = width == 0.
    if exc.any():
        p2 = np.maximum(0., (center[exc]-mass)*(center[exc]+mass))*red**2
        k = np.searchsorted(edge, np.sqrt(p2), side="right")-1
        total[0] += np.bincount(k, 2*n[exc], minlength=len(edge)-1)
        total[1] += np.bincount(k, 2*n[exc]*p2, minlength=len(edge)-1)
    keep = ~exc
    c, w, n, a = (x[keep, None] for x in (center, width, n, a))
    # Energy thresholds at the production epoch for the fixed output p/T grid.
    cuts = np.hypot(mass, edge/red)
    x = np.clip((cuts[None,:]-c)/w, -1., 1.)
    d = (x[:,1:]-x[:,:-1])/2
    z = (x[:,1:]+x[:,:-1])/2
    ii,jj = np.nonzero(d)
    z, d, c, w, n, a = z[ii,jj], d[ii,jj], c[ii,0], w[ii,0], n[ii,0], a[ii,0]
    for t, b in zip(*GL3):
        y = z+d*t
        g = (n+a*(3*y*y-1)/2)*d*b
        assert np.all(g>=-1e-200)
        E = c+w*y
        p2 = np.maximum(0.,(E-mass)*(E+mass))*red**2
        total[0] += np.bincount(jj, g, minlength=len(edge)-1)
        total[1] += np.bincount(jj, g*p2, minlength=len(edge)-1)
    assert abs(total[0].sum()/(2*np.sum(isotropic))-1.)<2e-12
    return total


def positive_range_sum(lo,hi,weight,nb):
    """Add nonnegative weights on half-open integer ranges without cancellation."""
    size=1<<(nb-1).bit_length()
    if not np.all((0<=lo)&(lo<=hi)&(hi<=nb)&(weight>=0.)):
        raise ValueError("Invalid positive half-open range")
    tree=np.zeros(2*size)
    a,b=lo+size,hi+size
    while True:
        keep=a<b
        if not keep.any():
            break
        a,b,weight=a[keep],b[keep],weight[keep]
        left=(a%2)==1
        right=(b%2)==1
        np.add.at(tree,a[left],weight[left])
        np.add.at(tree,b[right]-1,weight[right])
        a=(a+1)//2
        b=b//2
    k=1
    while k<size:
        tree[2*k:4*k]+=np.repeat(tree[k:2*k],2)
        k*=2
    return tree[size:size+nb]


def decay_bins(edge,mass,center,width,W,red):
    """Exact positive angle bins in near O(n log nb), not O(n nb)."""
    nb=len(edge)-1
    total=np.zeros((2,nb))
    cuts=np.hypot(mass,edge/red)
    low=np.maximum(mass,center-width)
    high=np.maximum(mass,center+width)
    ia=np.searchsorted(cuts,low,side="right")-1
    ib=np.searchsorted(cuts,high,side="right")-1

    def add(i,c,n,w,x,y):
        d=w*(y-x)
        E=c+w*x
        v=np.maximum(0.,(E-mass)*(E+mass)+E*d+d*d/3)*red**2
        g=n*(y-x)/2
        total[0]+=np.bincount(i,g,minlength=nb)
        total[1]+=np.bincount(i,g*v,minlength=nb)

    same=ia==ib
    add(ia[same],center[same],W[same],width[same],-1.,1.)
    a,b,c,w,n=(x[~same] for x in (ia,ib,center,width,W))
    if np.any(w<=0.):
        raise ValueError("Zero-width interval split between momentum bins.")
    L=np.clip((cuts[a+1]-c)/w,-1.,1.)
    R=np.clip((cuts[b]-c)/w,-1.,1.)
    add(a,c,n,w,-1.,L)
    add(b,c,n,w,R,1.)
    if nb>1 and a.size:
        z=positive_range_sum(a+1,b,n/(2*w),nb)
        assert z[-1]==0.
        l=cuts[:-2]
        d=np.diff(cuts[:-1])
        g=z[:-1]*d
        v=np.maximum(0.,(l-mass)*(l+mass)+l*d+d*d/3)*red**2
        total[0,:-1]+=g
        total[1,:-1]+=g*v
    assert abs(total[0].sum()/W.sum()-1.)<2e-12
    return total


def deposit(fold,T,edge,scale,angle_order):
    q,p,v,n,a=fold.C[0].event_moments(T)
    red=scale/TF
    direct=bins(edge,m,q/2,
                p*np.sqrt(v)[:,None]/2,n,a,red)
    E, W=fold.C[1].particles(T,angle_order)
    width=np.sqrt(np.maximum(0.,(E-M)*(E+M)))*P/M
    scalar=decay_bins(edge,mp,E*Es/M,width,W,red)
    paired=decay_bins(edge,m,E*(M-Es)/M,width,W,red)
    assert abs(scalar[0].sum()/paired[0].sum()-1.)<2e-12
    return np.array([scalar,direct+paired])


def resolve(fold,edge,tr,tmax,nt,steps,ang):
    r,l,hi=np.log([tr,TF,tmax])
    Y=np.zeros((2,2,2,len(edge)-1))
    x,w=roots_legendre(nt)
    for i,(lo,end) in enumerate(((l,r),(r,hi))):
        if end<=lo:
            continue
        grid=np.unique(np.clip(np.r_[np.linspace(lo,end,steps),fold.bg.sg.x],lo,end))
        for low,up in zip(grid[:-1],grid[1:]):
            for u,b in zip((low+up)/2+(up-low)/2*x,(up-low)/2*w):
                if i==0:
                    scale=(fold.s(l)/fold.s(u))**(1/3)
                    factor=b*(3+fold.grad(u))/(3*fold.H(u)*fold.s(u))
                else:
                    scale=(fold.s(l)/fold.s(r))**(1/3)*np.exp(-8*(u-r)/3)
                    factor=b*(8/3)*np.exp(-12*(u-r))/(fold.H(r)*fold.s(r))
                Y[i]+=factor*deposit(fold,np.exp(u),edge,scale,ang)
    return Y


def stream_bounds(bg, edge, count, p2, mass, order):
    ai=float(bg.scale(TF))
    finite=edge[:-1]
    distances=np.r_[0.,bg.streams(finite[1:]*TF*ai,mass,ai,bg.aeq,order)]
    n=count.sum()
    ok=count>0.
    local=np.divide(p2,count,out=np.zeros_like(p2),where=ok)
    # D is monotone and concave in x=(p/T)^2, by the pointwise
    # velocity sqrt(x/(x+C)); thus a chord is a lower bound in each
    # finite bin and Jensen's inequality gives an upper bound.
    weight=np.clip((local[:-1]-finite[:-1]**2)/np.diff(finite**2),0.,1.)
    lower=distances.copy()
    lower[:-1]+=weight*np.diff(distances)
    upper=bg.streams(np.maximum(np.sqrt(local[ok]),1e-30)*TF*ai,mass,ai,bg.aeq,order)
    whole=bg.streams(np.array([np.sqrt(p2.sum()/n)*TF*ai]),mass,ai,bg.aeq,order)[0]
    cf=np.cumsum(count)/n
    idx=np.searchsorted(cf,[.1,.5,.9],side="left")
    assert np.all(idx+1<len(distances))
    return dict(mean_Mpc=[float(lower@count/n),float(upper@count[ok]/n)],
                global_Jensen_upper_Mpc=float(whole),
                quantile_brackets_Mpc=np.column_stack([distances[idx],distances[idx+1]]).tolist())


def check_partial_angles():
    # Uniform-in-E angle and a quadratic angular density, checked by
    # exact primitives on the three disjoint subintervals.
    edges=np.array([0.,1.,2.,np.inf])
    m0=.5
    cut=np.clip((np.hypot(edges,m0)-2)/1.5,-1.,1.)
    v=bins(edges,m0,np.array([2.]),np.array([1.5]),np.array([.7]),np.array([.3]),1.)
    density=np.polynomial.Polynomial([.55,0.,.45])
    p2=np.polynomial.Polynomial([3.75,6.,2.25])
    ref=np.array([np.diff(f.integ()(cut)) for f in (density,density*p2)])
    v2=decay_bins(edges,m0,np.array([2.]),np.array([1.5]),np.array([1.4]),1.)
    expected=np.array([.7*np.diff(cut),.7*np.diff(p2.integ()(cut))])
    err=max(np.max(abs(v-ref)),np.max(abs(v2-expected)))
    assert err<1e-14
    return float(err)


def run(eos, reference):
    ref=json.loads(Path(reference).read_text())
    eos_sha=hashlib.sha256(Path(eos).read_bytes()).hexdigest()
    assert ref["eos_sha256"]==eos_sha
    rows=ref["controls"][-1]["rows"]
    F=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    bg=Background(eos)
    probes=np.array([.005,.01,.025,.05,.1,.2,.5,1.,2.,4.,8.,16.])
    edge=np.unique(np.r_[0.,probes,np.geomspace(1e-5,110,65),np.inf])
    descriptors=[((7,18,14),4,19,18),((8,23,17),6,28,28)]
    variants=[]
    for order,nt,steps,ang in descriptors:
        fold=SourceFold(eos,order,121,4)
        variants.append([resolve(fold,edge,r["T_R_GeV"],r["Tmax_GeV"],nt,steps,ang) for r in (rows[0],rows[-1])])
    out=[]
    ci=np.searchsorted(edge,probes)-1
    for row,low,high in zip((rows[0],rows[-1]),*variants):
        total=high.sum(axis=0)
        counts=total[:,0].sum(axis=-1)
        rr=np.sqrt(total[:,1].sum(axis=-1)/counts)
        assert np.max(abs(rr/row["combined_rms_momentum_over_1MeV"]-1.))<6e-5
        calc=counts*np.array([mp,m])*F
        ref_omega=[row["scalar_prompt_fixed_count_Omega_h2"],ref["target_stable_Omega_h2"]]
        assert np.max(abs(calc/ref_omega-1.))<6e-5
        pfs=high[1,:,0].sum(axis=-1)/counts
        assert np.max(abs(pfs-row["pre_fraction_of_scalar_and_stable_number"]))<2e-5
        def cdf(x):
            v=x.sum(axis=0)[:,0]
            return np.cumsum(v,axis=-1)[:,ci]/v.sum(axis=-1)[:,None]
        cdiff=np.max(abs(cdf(low)-cdf(high)),axis=-1)
        assert np.max(cdiff)<.003
        species=[]
        for j,mass in enumerate((mp,m)):
            v,z=total[j]
            s=stream_bounds(bg,edge,v,z,mass,8)
            s0=stream_bounds(bg,edge,v,z,mass,5)
            assert abs(np.array(s["mean_Mpc"])/s0["mean_Mpc"]-1.).max()<2e-7
            assert s["mean_Mpc"][0]<=s["mean_Mpc"][1]<=s["global_Jensen_upper_Mpc"]*(1+1e-12)
            assert v[-1]/counts[j]<1e-9
            species.append(dict(Omega_h2=float(calc[j]),pre_fraction=float(pfs[j]),rms_p_over_T=float(rr[j]),
                max_absolute_probe_CDF_change=float(cdiff[j]),CDF_at_probes=cdf(high)[j].tolist(),
                pre_share_of_all_number_below_probes=(np.cumsum(high[1,j,0])[ci]/counts[j]).tolist(),
                post_share_of_all_number_below_probes=(np.cumsum(high[0,j,0])[ci]/counts[j]).tolist(),
                streaming=s))
        def compact(a):
            return [float(f"{x:.9g}") for x in a]
        out.append(dict(T_R_GeV=row["T_R_GeV"],Tmax_GeV=row["Tmax_GeV"],species=species,
            count_fraction_per_epoch_species_bin=[[compact(a) for a in x] for x in high[:,:,0]/counts[None,:,None]],
            p_over_T_squared_fractional_number_moment_per_epoch_species_bin=[[compact(a) for a in x] for x in high[:,:,1]/counts[None,:,None]]))
    return dict(anchor_temperature_GeV=TF,epoch_order=["post","pre"],species_order=["86 eV scalar","stable fermion"],
        eos_sha256=eos_sha,reference_sha256=hashlib.sha256(Path(reference).read_bytes()).hexdigest(),
        controls=descriptors,partial_angle_polynomial_absolute_error=check_partial_angles(),boundary_p_over_T=edge[:-1].tolist(),final_bin_extends_to_infinity=True,
        cdf_thresholds_p_over_T=probes.tolist(),cases=out,
        scope="Positive quadrature for the prescribed cooling branch of the vacuum, prompt-only reference. No collisions, inverse, scalar-pair production, heating branch or backreaction. Streaming is a conditional geodesic statistic from 1 MeV to equality, not the linear transfer function.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("reference")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.reference),separators=(",",":"),allow_nan=False))
