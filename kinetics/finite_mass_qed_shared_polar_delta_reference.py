"""The two positive scalar-latest absorptive cuts share one closed kernel."""
import gzip,json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import Q,intervals
from finite_mass_qed_charged_complex_guard_reference import chart


def audit(d,E,track):
    M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
    Es,lo,hi,*_=intervals(M,s,m,Q)
    A=(M*M-m*m+s*s)/2
    Z=(M*M+m*m-s*s)/2
    H=np.sqrt((A-M*s)*(A+M*s))
    K=np.pi*Z/(2*H)
    C=(M+m)**2/(Z+M*m)
    const=-np.pi*C/(4*Q)
    rho=3e-4
    sides=[]
    def kernel(δ):
        return const+K/abs(δ)
    for name,label,ma,e0 in (('parent','parent',M,E),('fermion','daughter',m,E-Es)):
        dat=chart(M,m,s,E,label)[2]
        grid=np.hypot(ma,d[f'endpoint_{name}_momentum_edges'])-e0
        if track==name:
            ix=np.argmin(abs(grid))
            if abs(grid[ix])>1e-14:
                raise ArithmeticError('The tuned event is not on the requested cell edge.')
            grid[ix]=0.
        f=np.asarray(d[f'endpoint_{name}_occupations'])
        sides.append((grid,f,dat))
    (gp,fp,P),(gb,fb,B)=sides
    def read(grid,f,q,side='right'):
        i=np.searchsorted(grid,q,side=side)-1
        return f[i] if 0<=i<len(f) else 0
    def shape(δ,kind):
        p,b=read(gp,fp,δ),read(gb,fb,δ)
        return {'parent':p,'daughter':b,'net':p-b,'identical':0.}[kind]
    def kernel_cut(dat,δ):
        _,a,x,root,g,R,r=dat(δ)
        return -np.pi/(4*Q)*(C-2*Z/(R*np.hypot(x,a))),root
    def kernel_dif(δ):
        p,rp=kernel_cut(P,δ)
        d,rd=kernel_cut(B,δ)
        if max(abs(rp),abs(rd))>=1.:
            raise ArithmeticError('An interior polar root was unexpectedly lost.')
        s=kernel(δ)
        return float(max(abs(p-s),abs(d-s))/max(1.,abs(s))),float(abs(p-d))
    guard=[*np.geomspace(1e-11,3e-4,9),-1e-11,-1e-6,-2e-4]
    test=[v for v in guard if lo-E<v<hi-E]
    test.extend(lo+u*(hi-lo)-E for u in (.01,.19,.31,.72,.92,.99))
    defect=max(kernel_dif(q)[0] for q in test)
    # Independently perform the azimuth integral in the scalar rest frame.
    pe=H/s
    # Equivalent and carefully scaled: E=Es*A/s²+(Q pe/s)*ce.
    ce=(E-Es*A/(s*s))/(Q*pe/s)
    az=[]
    for u in (.1,.35,.7,.9):
        δ=lo+u*(hi-lo)-E
        ci=ce+δ/(Q*pe/s)
        if abs(ce)>=1 or abs(ci)>=1 or abs(ci-ce)<.03:
            continue
        g=np.sqrt((1-ce*ce)*(1-ci*ci))
        bare=quad(lambda phi:-1/(2*pe*pe*(1-ce*ci-g*np.cos(phi))),
                  0,2*np.pi,epsabs=1e-11,epsrel=2e-12)[0]/(2*np.pi)
        expect=-Q/(2*H*abs(δ))
        az.append(abs(bare/expect-1))
    if defect>2e-9 or not az or max(az)>3e-10:
        raise ArithmeticError(f'Common kernel failed the two-frame independent angular test: {defect}, {az}.')
    support=(lo-E,hi-E)
    exterior=[]
    for δ in (support[0]-.001,support[1]+.001):
        r0,r1=kernel_cut(P,δ)[1],kernel_cut(B,δ)[1]
        exterior.append(float(min(abs(r0),abs(r1))-1))
    if min(exterior)<=0:
        raise ArithmeticError('The two polar cuts did not share their exact exterior boundaries.')
    knots=sorted({*support,0.,*[q for q in (*gp,*gb) if support[0]<q<support[1]]})
    def exact(kind,eps):
        cuts=[]
        for a,b in zip(knots[:-1],knots[1:],strict=True):
            weight=shape((a+b)/2,kind)
            if eps is not None:
                a,b=(a,min(b,-eps)) if b<=0 else (max(a,eps),b)
            if b<=a: continue
            s=-1 if b<=0 else 1
            cuts.append(weight*(const*(b-a)+K*s*np.log((abs(b) or rho)/(abs(a) or rho))))
        return fsum(cuts)
    jdiff=[read(gp,fp,0.,i)-read(gb,fb,0.,i) for i in ('left','right')]
    J=K*sum(jdiff)
    finite=exact('net',None)
    decomposition=abs(exact('parent',None)-exact('daughter',None)-finite)
    if decomposition>2e-11:
        raise ArithmeticError('Separated finite polar-delta integrations did not close.')
    coeff=[]
    n,w=roots_legendre(41)
    for eps in (1e-6,1e-8,1e-10,1e-12):
        raw=exact('net',eps)
        direct=[]
        equal=[]
        for a,b in zip(knots[:-1],knots[1:],strict=True):
            a,b=(a,min(b,-eps)) if b<=0 else (max(a,eps),b)
            if b<=a: continue
            l,h=np.log(min(abs(a),abs(b))),np.log(max(abs(a),abs(b)))
            q=(-1 if b<0 else 1)*np.exp((l+h)/2+(h-l)/2*n)
            vals=[]
            check_equal=[]
            for δ in q:
                p=read(gp,fp,δ)
                d=read(gb,fb,δ)
                P0,_=kernel_cut(P,δ)
                B0,_=kernel_cut(B,δ)
                vals.append(float(abs(δ)*(p*P0-d*B0)))
                eq=np.exp(-((δ/(hi-lo))**2))
                check_equal.append(float(abs(δ)*eq*(P0-B0)))
            direct.append((h-l)/2*(w@vals))
            equal.append((h-l)/2*(w@check_equal))
        # The regular constant contributes the explicitly omitted width.
        expected=finite+J*np.log(rho/eps)-const*eps*sum(jdiff)
        diff=max(abs(fsum(direct)-raw),abs(raw-expected))
        if diff>2e-8:
            raise ArithmeticError(f'Integral of paired causal kernels did not reproduce the closed expression: {diff}.')
        val_equal=abs(fsum(equal))
        if val_equal>2e-8:
            raise ArithmeticError('Independent comparison with identical smooth occupations did not cancel.')
        coeff.append(dict(hard_radial_guard=eps,integrated_raw=raw,
                          independent_radial_replay_defect=diff,identical_smooth_residual=val_equal))
    return dict(parent_energy=E,tuned_step=track,constant_kernel_term=const,photon_odd_inverse_abs_prefactor=K,
                occupation_difference_across_coincidence=jdiff,
                full_aperture_log_coefficient=J,chosen_finite_part=finite,
                normalized_common_kernel_audit=defect,scalar_rest_azimuth_defect=max(az),
                exterior_root_abs_excess=exterior,ladder=coeff)


def run():
    out=[]
    for g in ('fine','third'):
        root=Path(__file__).parent/'results'
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
        Es,lo,hi,*_=intervals(M,s,m,Q)
        name,ma,shift=('fermion',m,Es) if g=='fine' else ('parent',M,0.)
        e=np.hypot(ma,d[f'endpoint_{name}_momentum_edges'][1:-1])+shift
        ok=np.flatnonzero((e>lo+4e-4)&(e<hi-4e-4))
        idx=ok[np.argmax(abs(np.diff(d[f'endpoint_{name}_occupations'])[ok]))]
        vals=[audit(d,lo+t*(hi-lo),None) for t in (.25,.5,.75)]
        vals.append(audit(d,float(e[idx]),name))
        out.append(dict(grid=g,controls=vals))
    return dict(scope='Paired positive charged scalar-latest polar delta over common vacuum Born support; no common physical thermal factor or damping assigned.',radial_finite_part_reference=3e-4,rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
