"""Two-centroid charged deposition at exact energies in the archived 1<->2 aperture."""
import json
from math import fsum

import numpy as np
from scipy.special import expit,roots_legendre

from early_pole_decay_mode_clock_reference import G,KAPPA,X
from early_pole_endpoint_stencil_remap_pilot import LABELS,SPECS,archives,load,parts,remap,MassiveFermion
from early_pole_live_parent_energy_reference import state_moments


def chart(d,x,m,order,*,histogram_cuts=False):
    Es=np.hypot(m,d.d.p)
    Ep=d.parent_energies(x)
    Ed=d.daughter_energies(x)
    mb=d.daughter_mass(x)
    if x<=m+mb:
        raise ValueError('This pilot requires open vacuum decay')
    D=np.sqrt((x*x-(m+mb)**2)*(x*x-(m-mb)**2))/2
    A=(x*x+m*m-mb*mb)/2
    lo=(x*x*Es*Es+D*D)/(A*Es+d.d.p*D)
    hi=lo+2*d.d.p*D/(m*m)
    # Work in parent energy; both projection grids cut every integration
    # interval in which the bracketing daughter and parent nodes stay fixed.
    left=np.maximum.reduce([lo,np.full_like(lo,Ep[0]),Es+Ed[0]])
    right=np.minimum.reduce([hi,np.full_like(hi,Ep[-1]),Es+Ed[-1]])
    valid=left<right
    clip_right=np.maximum(left,right)
    nt=np.concatenate((np.broadcast_to(Ep,(d.ns,len(Ep))),Es[:,None]+Ed),axis=1)
    if histogram_cuts:
        exact_ep=np.hypot(d.pe,x)
        exact_ed=np.hypot(d.d.fed,mb)
        nt=np.column_stack([nt,np.broadcast_to(exact_ep,(d.ns,len(exact_ep))),Es[:,None]+exact_ed])
    raw=np.sort(np.column_stack([left,right,np.clip(nt,left[:,None],clip_right[:,None])]),axis=1)
    low,high=raw[:,:-1],raw[:,1:]
    si,slot=np.where((high>low)&valid[:,None])
    half=(high[si,slot]-low[si,slot])/2
    cen=(low[si,slot]+high[si,slot])/2
    z,w=roots_legendre(order)
    S=np.repeat(si,order)
    P=(cen[:,None]+half[:,None]*z).ravel()
    daughter=P-Es[S]
    J=np.searchsorted(Ed,daughter,side='right')
    K=np.searchsorted(Ep,P,side='right')
    if np.any(J<1)|np.any(K<1)|np.any(J>=len(Ed))|np.any(K>=len(Ep)):
        raise ArithmeticError('A quadrature node escaped its projection span')
    ad=(daughter-Ed[J-1])/(Ed[J]-Ed[J-1])
    ap=(P-Ep[K-1])/(Ep[K]-Ep[K-1])
    if np.min([ad.min(),ap.min()])<0 or np.max([ad.max(),ap.max()])>1:
        raise ArithmeticError('Negative charge deposition weights')
    N=KAPPA*G*x*((x+mb)**2-m*m)/(Es*d.d.p)*d.d.w
    W=(N[si,None]*half[:,None]*w).ravel()
    weightedD=(1-ad)*Ed[J-1]+ad*Ed[J]
    weightedP=(1-ap)*Ep[K-1]+ap*Ep[K]
    proj=float(np.max(abs((Es[S]+weightedD-weightedP)/weightedP)))
    if proj>4e-14 or np.min(W)<0:
        raise ArithmeticError('Event-level energy conservation lost')
    rate=lambda x,y: np.log1p(np.exp(-x))-np.log1p(np.exp(-y))
    thermal=N@rate(lo,hi)
    retained=N@np.where(valid,rate(left,right),0.)
    return (Es,Ed,Ep,S,J,K,ad,ap,W),dict(
        quadrature_events=len(W),maximum_relative_event_energy_error=proj,
        omitted_fraction_for_thermal_parent_seed=float(1.-retained/thermal))


def event(chart,y,d):
    Es,Ed,Ep,i,j,k,a,b,W=chart
    fs,fd,fp=np.split(y,(d.ns,d.iP))
    p,pc=d.parent_products(fp[k-1],fp[k],1.-b,b,False)
    f,fc=d.parent_products(fd[j-1],fd[j],1.-a,a,False)
    forward=W*p*(1.+fs[i])*fc
    backward=W*pc*fs[i]*f
    net=forward-backward
    vals=np.r_[np.bincount(i,net,minlength=d.ns),
               np.bincount(j-1,(1.-a)*net,minlength=d.nf)+np.bincount(j,a*net,minlength=d.nf),
              -np.bincount(k-1,(1.-b)*net,minlength=d.np)-np.bincount(k,b*net,minlength=d.np)]
    return vals,forward,backward



def exact_step_on_common_span(d,x,m,y):
    # Independent point evaluation of a constant-in-momentum native
    # histogram. This computes the scalar source only, without depositing
    # charge energy into cell centroids. It is a reference, not an ODE.
    Es=np.hypot(m,d.d.p)
    Ep=d.parent_energies(x)
    Ed=d.daughter_energies(x)
    ef=np.hypot(d.d.fed,d.daughter_mass(x))
    ep=np.hypot(d.pe,x)
    mb=d.daughter_mass(x)
    D=np.sqrt((x*x-(m+mb)**2)*(x*x-(m-mb)**2))/2
    A=(x*x+m*m-mb*mb)/2
    lo=(x*x*Es*Es+D*D)/(A*Es+d.d.p*D)
    hi=lo+2*d.d.p*D/(m*m)
    l=np.maximum.reduce([lo,np.full_like(lo,Ep[0]),Es+Ed[0]])
    r=np.minimum.reduce([hi,np.full_like(hi,Ep[-1]),Es+Ed[-1]])
    span=r>l
    cuts=np.concatenate((np.broadcast_to(ep,(d.ns,len(ep))),Es[:,None]+ef),axis=1)
    sorted_cuts=np.sort(np.column_stack([l,r,np.clip(cuts,l[:,None],np.maximum(l,r)[:,None])]),axis=1)
    width=np.diff(sorted_cuts,axis=1)
    si,ell=np.where((width>0)&span[:,None])
    P=(sorted_cuts[si,ell]+sorted_cuts[si,ell+1])/2
    fs,fd,fp=np.split(y,(d.ns,d.iP))
    U=fp[np.searchsorted(ep,P)-1]
    V=fd[np.searchsorted(ef,P-Es[si])-1]
    n=fs[si]
    W=d.d.w[si]*width[si,ell]*KAPPA*G*x*((x+mb)**2-m*m)/(Es[si]*d.d.p[si])
    return np.array([fsum(W*U*(1+n)*(1-V)),fsum(W*(1-U)*n*V)])

def case(d,x,y,orders=(4,8)):
    _,unc,pars,m=state_moments(d,x,y)
    a0=d.w*d.rhs(np.log(x/X),y)
    S=fsum(d.d.w*y[:d.ns])
    H=np.sqrt(d.radiation/(d.radiation+unc+pars))
    out=[]
    exact=exact_step_on_common_span(d,x,m,y)*H/S
    for order in orders:
        ch,sz=chart(d,x,m,order)
        es,ed,ep,*_=ch
        thermal=np.r_[1/np.expm1(es),expit(-ed),expit(-ep)]
        _,fo,inv=event(ch,thermal,d)
        bal=np.max(abs(fo-inv)/np.maximum(fo,1e-250))
        vals,f,b=event(ch,y,d)
        Q=np.r_[es,ed,ep]
        residual=abs(Q@vals)/sum(abs(vals*Q))
        parts=[fsum(v) for v in np.split(vals,(d.ns,d.iP))]
        pair=max(abs(parts[0]-parts[1]),abs(parts[0]+parts[2]))
        if bal>7e-13 or residual>3e-13 or pair>3e-13*sum(map(abs,parts)):
            raise ArithmeticError(('Thermal or conserving gate failed',bal,residual,pair))
        out.append(dict(order=order,thermal_event_relative_defect=float(bal),energy_relative_defect=residual,
                        pair_absolute_defect=pair,net_clock=parts[0]*H/S,gross_forward=fsum(f)*H/S,
                        gross_inverse=fsum(b)*H/S,**sz))
    perturb=max(abs(out[-1][k]/out[-2][k]-1) for k in ('net_clock','gross_forward','gross_inverse'))
    if perturb>2e-5:
        raise ArithmeticError(('Opt-in charge overlap did not refine',perturb,out))
    ch,_=chart(d,x,m,orders[-1])
    _,_,_,ii,jj,kk,*_=ch
    boundary=[]
    # Test an actual populated charged cell with both boundary values.
    for species,bound in ((d.ns,np.bincount(jj,minlength=d.nf).argmax()),
                          (d.iP,np.bincount(kk,minlength=d.np).argmax())):
        index=species+bound
        uv=[]
        for limit in (0.,1.):
            p=y.copy()
            p[index]=limit
            dx=event(ch,p,d)[0][index]
            if dx*(1.-2*limit)<0.:
                raise ArithmeticError(('Wrong Pauli boundary direction',index,limit,dx))
            uv.append(float(dx))
        boundary.append(dict(cell=int(bound),sources_at_empty_and_full=uv))
    return dict(native_clock=fsum(a0[:d.ns])/S,double_centroid=out,mesh_relative_defect=perturb,
                exact_step_on_shared_support_gross=exact.tolist(),
                exact_step_on_shared_support_net=float(exact[0]-exact[1]),
                occupied_charge_boundary_controls=boundary)


def run():
    refs={g:load(g)[0] for g in SPECS}
    out=[]
    for mesh,shape in SPECS.items():
        d=MassiveFermion(*shape,refs[mesh]['mass_ratio'],refs[mesh]['elastic_rate'],True)
        x=refs[mesh]['x_end']
        for g,r in refs.items():
            ys=[r['endpoint_'+j+'_occupations'] for j in LABELS]
            y=np.concatenate([remap(o,v,n)[0] for o,v,n in zip(archives(r),ys,parts(d),strict=True)])
            out.append(dict(mesh=mesh,archive=g,result=case(d,x,y)))
    return out


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
