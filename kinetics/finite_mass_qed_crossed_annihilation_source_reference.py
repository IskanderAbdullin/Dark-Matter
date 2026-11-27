"""Fixed-q charge-crossed source pilot with an explicit anti-daughter proxy."""
from functools import lru_cache
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre, expit
from scipy.optimize import brentq

from early_pole_decay_mode_clock_reference import Q, intervals
from early_pole_hard_yukawa_photoproduction_reference import invariant_form


@lru_cache(None)
def gl(n):
    return roots_legendre(n)


def coefficient(M,m,s2,S,J):
    d=(M*M+m*m-S)/2
    return 2*((M+m)**2-s2), 4*((J+d)**2+d*(d+2*M*m))/J


def azimuth(M,m,s2,S,J,r,K,sp,x):
    """Exact average of the crossed spin sum over the pair azimuth."""
    c0=np.sqrt((1-sp)*(1+sp))
    Ap=K*(M*M/(np.hypot(r,M)+r)+r*(1-c0*x))
    Ab=K*(m*m/(np.hypot(r,m)+r)+r*(1+c0*x))
    lateral=K*r*sp*np.sqrt((1-x)*(1+x))
    P=np.sqrt((Ap-lateral)*(Ap+lateral))
    B=np.sqrt((Ab-lateral)*(Ab+lateral))
    U,V=coefficient(M,m,s2,S,J)
    return U*(M*M*Ap/P**3+m*m*Ab/B**3)+V*(1/P+1/B)


def angle_primitive(M,m,s2,S,J,r,K,sp,a,b):
    """Exact azimuth average followed by exact lab-boost cosine integration.

    This branch assumes nonzero finite masses and that neither shifted peak
    center is inside [-1,1], checked explicitly for the fixed cold q.
    """
    c0=np.sqrt((1-sp)*(1+sp))
    d0=sp*sp/(1+c0)
    beta=K*r
    U,V=coefficient(M,m,s2,S,J)

    def one(mass,left,right):
        alpha=K*np.hypot(r,mass)
        v=K*mass*mass/(np.hypot(r,mass)+r)-alpha*d0
        va=v+beta*(1-left)
        vb=v+beta*(1-right)
        if min(va.min(),vb.min())<=0.:
            raise ArithmeticError('Azimuth peak entered the fast polar chart')
        Da,Db=np.hypot(va,sp*K*mass),np.hypot(vb,sp*K*mass)
        sumv=va+vb
        dx=right-left
        u=dx*(1+sumv/(Da+Db))/(vb+Db)
        # log1p(beta*u)/beta -> u continuously at the pair threshold.
        I1=np.divide(np.log1p(beta*u),beta,out=u.copy(),where=beta>0.)
        I2=dx*sumv/(Da*Db)*(alpha*sp*sp/(va*Db+vb*Da)+c0/(Da+Db))
        return U*mass*mass*I2+V*I1
    return one(M,a,b)+one(m,-b,-a)


def value(data,T,kcap,nk,nphi,nmu,nouter,equilibrium=False):
    M=data['x_end']; m=M*data['mass_ratio']
    s=data['checkpoints'][-1]['scalar_pole_over_T']
    es=np.hypot(Q,s)
    p_edges,b_edges=[np.hypot(x,data[f'endpoint_{n}_momentum_edges'])
                     for x,n in ((M,'parent'),(m,'fermion'))]
    fp,fb=map(lambda k:np.r_[0.,data[f'endpoint_{k}_occupations'],0.],('parent','fermion'))
    fs=data['endpoint_scalar_occupations'][np.searchsorted(data['scalar_momentum_edges'],Q)-1]
    cm,wmu=gl(nmu); x,w=gl(nouter)
    if nphi>0:
        xp,wp=gl(nphi)
        zp=np.cos(np.pi*(xp+1.))
    rows=[]
    if equilibrium:
        fs=1/np.expm1(es/T)
    for ca,wm in zip(cm,wmu,strict=True):
        k0=((M+m)-s)*((M+m)+s)/(2*(es-Q*ca))
        def params(k):
            S=s*s+2*k*(es-Q*ca)
            t=np.sqrt(S)
            r=np.sqrt(np.maximum((S-(M+m)**2)*(S-(M-m)**2),0.))/(2*t)
            K=k*(es-Q*ca)/t
            lab=k+es
            P=np.sqrt((k+Q*ca)**2+Q*Q*(1-ca*ca))
            EC=lab*(S+M*M-m*m)/(2*S)
            R=P*r/t
            lng=Q*np.sqrt(1-ca*ca)*k/(P*K)
            return S,t,r,K,lab,EC,R,lng
        def batched(k):
            # Integrate all constant charged pieces for up to 128 photon
            # nodes at once; remove zero-width cells by their exact weights.
            S,t,r,K,lab,EC,R,sp=params(k[:,None])
            knots=np.sort(np.clip(np.concatenate(
                [np.broadcast_to([-1.,1.],(len(k),2)),(p_edges-EC)/R,(lab-b_edges-EC)/R],
                  axis=1),-1.,1.),axis=1)
            a,b=knots[:,:-1],knots[:,1:]
            e=EC+R*(a+b)/2
            f=fp[np.searchsorted(p_edges,e)]
            g=fb[np.searchsorted(b_edges,lab-e)]
            wt=angle_primitive(M,m,s*s,S,K*t,r,K,sp,a,b)
            if wt.min()<-1e-8:
                raise ArithmeticError('Negative crossed kernel in positive phase measure')
            n=1/np.expm1(k/T)
            gain=np.sum(wt*f*g,axis=1)*(1+n)
            loss=np.sum(wt*(1-f)*(1-g),axis=1)*n
            loss_without_anti=np.sum(wt*(1-f),axis=1)*n
            phase=k*(r/t)[:,0]/(128*np.pi**3*es)
            return phase[:,None]*np.stack([gain,loss,loss_without_anti],axis=1)
        def inner(k):
            S,t,r,K,lab,EC,R,sp=params(k)
            knots=np.array([-1.,1.]) if equilibrium else np.unique(np.clip(np.r_[-1.,1.,(p_edges-EC)/R,(lab-b_edges-EC)/R],-1.,1.))
            a,b=knots[:-1],knots[1:]
            e_mid=EC+R*(a+b)/2
            f=fp[np.searchsorted(p_edges,e_mid)]
            g=fb[np.searchsorted(b_edges,lab-e_mid)]
            n=1/np.expm1(k/T)
            if nphi<0 and not equilibrium:
                wt=angle_primitive(M,m,s*s,S,K*t,r,K,sp,a,b)
            else:
                q0=(a[:,None]+b[:,None])/2+(b-a)[:,None]*x[None,:]/2
                if nphi>0:
                    c=np.sqrt(1-sp*sp)*q0[:,:,None]+sp*np.sqrt((1-q0)*(1+q0))[:,:,None]*zp[None,None,:]
                    aa=-K*(M*M/(np.hypot(r,M)+r)+r*(1-c))
                    cc=K*(m*m/(np.hypot(r,m)+r)+r*(1+c))
                    amp=-invariant_form(aa,cc,(M*M+m*m-S)/2,M,m)
                    az=amp@(wp/2)
                else:
                    az=azimuth(M,m,s*s,S,K*t,r,K,sp,q0)
                if equilibrium:
                    ep=EC+R*q0
                    fp_node,fb_node=expit(-ep/T),expit(-(lab-ep)/T)
                    m0=(b-a)/2
                    gain=m0@((az*fp_node*fb_node)@w)*(1+n)
                    loss=m0@((az*(1-fp_node)*(1-fb_node))@w)*n
                    loss_if_all_anti_empty=m0@((az*(1-fp_node))@w)*n
                else:
                    wt=(b-a)/2*(az@w)
            if not equilibrium:
                if wt.min()<-1e-8:
                    raise ArithmeticError('Negative crossed kernel')
                gain=wt@(f*g*(1+n))
                loss=wt@((1-f)*(1-g)*n)
                loss_if_all_anti_empty=wt@((1-f)*n)
            # Fixed scalar collision clock, one unaveraged charge assignment.
            phase=k*r/t/(128*np.pi**3*es)
            return phase*np.array([gain,loss,loss_if_all_anti_empty])
        # A log-mapped radial ladder resolves the threshold, bath tail, and
        # the broad archive occupation structure. Mesh convergence is explicit.
        z,wk=gl(nk)
        bounds=np.r_[0.,np.geomspace(1e-10,kcap-k0,95)]
        # Resolve topology changes where an allowed charged energy endpoint
        # crosses a histogram edge. Internal parent/daughter edge-pair
        # order swaps are additionally considered in the support.
        edges=set(bounds)
        for species, bar in (() if equilibrium else ((0,p_edges),(1,b_edges))):
            for sg in (-1.,1.):
                def chart(delta):
                    _,_,_,_,lab,EC,R,_=params(k0+delta)
                    return (lab-EC if species else EC)+sg*R
                ch=[chart(v) for v in bounds]
                for v,u,a0,b0 in zip(ch[:-1],ch[1:],bounds[:-1],bounds[1:],strict=True):
                    e=bar[(bar>min(v,u))&(bar<max(v,u))]
                    for E in e:
                        r=brentq(lambda zz:chart(zz)-E,a0,b0,xtol=4e-14,rtol=2e-14)
                        edges.add(r)
        if not equilibrium:
            PP,BB=np.broadcast_arrays(p_edges[:,None],b_edges[None,:])
            D=PP+BB-es-k0
            ok=(D>0.)&(D<kcap-k0)
            k=k0+D[ok]
            if len(k):
                S=s*s+2*k*(es-Q*ca)
                r=np.sqrt((S-(M+m)**2)*(S-(M-m)**2))/(2*np.sqrt(S))
                lab=k+es
                P=np.sqrt((k+Q*ca)**2+Q*Q*(1-ca*ca))
                EC=lab*(S+M*M-m*m)/(2*S)
                R=P*r/np.sqrt(S)
                edges.update(D[ok][abs(PP[ok]-EC)<R])
        bounds=np.array(sorted(edges))
        bounds=bounds[np.r_[True,np.diff(bounds)>8*np.spacing(k0+bounds[1:])]]
        whalf=np.diff(bounds)/2
        mid=(bounds[:-1]+bounds[1:])/2
        kn=(k0+mid[:,None]+whalf[:,None]*z).ravel()
        wr=(whalf[:,None]*wk).ravel()
        vec=np.zeros(3)
        if nphi<0 and not equilibrium:
            for start in range(0,len(kn),128):
                vec+=wr[start:start+128]@batched(kn[start:start+128])
        else:
            for k,v in zip(kn,wr,strict=True):
                vec+=v*inner(k)
        rows.append(wm*vec)
    G,L,L0=np.sum(rows,axis=0)
    if L0 < L or min(G,L)<0.:
        raise ArithmeticError('Positive crossed measure violated Pauli monotonicity')
    return dict(gain=float(G),loss=float(L),net=float(G+fs*(G-L)),scalar_f=float(fs),
                loss_if_all_anti_empty=float(L0),
                net_if_all_anti_empty=float(-fs*L0),
                net_change_per_equal_shape_scale=float((1+fs)*G+fs*(L0-L)))


def born(data):
    """Same fixed-q normalizations and one charge assignment as the crossed test."""
    M=data['x_end']; m=M*data['mass_ratio']
    s=data['checkpoints'][-1]['scalar_pole_over_T']
    es,lo,hi,span,_=intervals(M,s,m,Q)
    P=np.hypot(M,data['endpoint_parent_momentum_edges'])
    B=np.hypot(m,data['endpoint_fermion_momentum_edges'])+es
    points=np.unique(np.clip(np.r_[lo,hi,P,B],lo,hi))
    w=np.diff(points)
    E=(points[1:]+points[:-1])/2
    p=np.r_[0.,data['endpoint_parent_occupations'],0.][np.searchsorted(P,E)]
    b=np.r_[0.,data['endpoint_fermion_occupations'],0.][np.searchsorted(B,E)]
    f=data['endpoint_scalar_occupations'][np.searchsorted(data['scalar_momentum_edges'],Q)-1]
    M2=2*((M+m)**2-s*s)
    K=M2/(16*np.pi*Q*es)
    G=K*(w@(p*(1-b)))
    L=K*(w@(b*(1-p)))
    return dict(gain=float(G),loss=float(L),stimulus=float(G-L),net=float(G+f*(G-L)))


def inverse_tail_bound(data,T,kmin):
    M=data['x_end']; m=M*data['mass_ratio']
    s=data['checkpoints'][-1]['scalar_pole_over_T']
    es=np.hypot(Q,s)
    W0=M+m
    eta=1-(s/W0)**2
    amin=eta*np.array([M*M,m*m])/2
    U=2*(W0*W0-s*s)
    C=(M*M+m*m-s*s)/2
    B0=U*(np.array([M*M,m*m])/amin**2).sum()+(
        8*C*C/W0**2+2*s*s)*(1/amin).sum()/eta
    B1=4*(es+Q)*(1/amin).sum()/eta
    thermal=np.exp(-kmin/T)/(-np.expm1(-kmin/T))
    I1=T*(kmin+T)
    I2=T*(kmin*kmin+2*kmin*T+2*T*T)
    return float(thermal*(B0*I1+B1*I2)/(64*np.pi**3*es))


def angle_controls(data):
    M=data['x_end']; m=M*data['mass_ratio']
    s=data['checkpoints'][-1]['scalar_pole_over_T']
    es=np.hypot(Q,s)
    x,w=gl(320)
    z,z_w=gl(32)
    phi=np.cos(np.pi*(z+1.))
    max_angle=max_boost=max_phase=0.
    for k in (1.63,1.9,4.,10.,30.,68.):
        for ca in (-.91,.17,.81):
            S=s*s+2*k*(es-Q*ca)
            W=np.sqrt(S)
            K=k*(es-Q*ca)/W
            r=np.sqrt((S-(M+m)**2)*(S-(M-m)**2))/(2*W)
            lab=np.array([Q*np.sqrt(1-ca*ca),0.,k+Q*ca])
            P=np.linalg.norm(lab)
            sp=Q*np.sqrt(1-ca*ca)*k/(P*K)
            U0=(k+es)/W
            U=-lab/W
            dp=k*U[2]
            V=np.array([0.,0.,k])+U*(k+dp/(U0+1))
            zstar=U0*k+dp
            max_boost=max(max_boost,abs(zstar/K-1.),
                          abs((V@lab)/(P*zstar)-np.sqrt(1-sp*sp)))
            cov_phase=2*k*r/(128*np.pi**3*W*es)
            sigma_unit=r/(16*np.pi*S*K)
            moller=K*W/(es*k)
            flux_phase=sigma_unit*moller*k*k/(4*np.pi**2)
            max_phase=max(max_phase,abs(cov_phase/flux_phase-1.))
            for a,b in ((-1.,-.85),(-.85,.3),(.3,.98),(.98,1.)):
                tt=(a+b)/2+(b-a)/2*x[:,None]
                angle=np.sqrt(1-sp*sp)*tt+sp*np.sqrt((1-tt)*(1+tt))*phi[None,:]
                aa=-K*(np.hypot(r,M)-r*angle)
                bb=K*(np.hypot(r,m)+r*angle)
                by_gl=(b-a)/2*(w@(-invariant_form(aa,bb,(M*M+m*m-S)/2,M,m)@(z_w/2)))
                by_exact=angle_primitive(M,m,s*s,S,K*W,r,K,sp,np.array([a]),np.array([b]))[0]
                max_angle=max(max_angle,abs(by_gl/by_exact-1.))
    if max_angle>8e-10 or max_boost>8e-10 or max_phase>1e-14:
        raise ArithmeticError(('Independent angle/boost control',max_angle,max_boost))
    return dict(independent_GL_angle_max_relative_defect=float(max_angle),
                explicit_Lorentz_projection_defect=float(max_boost),
                unit_amplitude_phase_flux_relative_defect=float(max_phase))


def run():
    results=[]
    for g in ('fine','third'):
        with gzip.open(Path(__file__).parent/'results'/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:
            d=json.load(f)
        r=[]
        for nk,nphi,nmu,nouter in ((3,-1,4,4),(5,-1,6,4),(7,-1,8,4)):
            r.append(dict(mesh=[nk,nphi,nmu,nouter],values=value(d,1.,68.,nk,nphi,nmu,nouter)))
        control=[]
        for T in (.4,1.):
            z=value(d,T,40.,6,0,6,36,True)
            G,L,N,F=(z[k] for k in ('gain','loss','net','scalar_f'))
            ratio=abs(N)/((1+F)*G+F*L)
            if ratio>1e-12:
                raise ArithmeticError(('integrated charged-pair detailed-balance control',ratio))
            control.append(dict(T=T, G=G,L=L,net_over_turnover=ratio))
        tail_short=value(d,1.,40.,5,-1,6,4)
        bound=inverse_tail_bound(d,1.,68.)
        outlim=np.hypot(d['x_end'],d['endpoint_parent_momentum_edges'][-1])+np.hypot(
           d['x_end']*d['mass_ratio'],d['endpoint_fermion_momentum_edges'][-1])-np.hypot(Q,d['checkpoints'][-1]['scalar_pole_over_T'])
        if outlim>=68.:
            raise ArithmeticError('Forward histogram support touches the tail cap')
        fields=('gain','loss','net','net_if_all_anti_empty')
        mesh_error=max(abs(r[-1]['values'][k]/r[-2]['values'][k]-1.) for k in fields)
        if mesh_error>1e-3:
            raise ArithmeticError(('Crossed source refinement did not converge',g,mesh_error))
        base=born(d)
        ratio=(.005**2)*r[-1]['values']['net']/base['net']
        results.append(dict(grid=g,rows=r,short_cap_40=tail_short,equilibrium=control,
                            analytical_inverse_absolute_tail_bound_after_68=bound,
                            forward_exact_maximum_supported_k=float(outlim),
                            angular_controls=angle_controls(d),
                            final_mesh_max_relative_drift=mesh_error,
                            antibath_shape_scale_pauli_max=1./max(d['endpoint_fermion_occupations']),
                            antibath_shape_scale_required_for_neutral_rate=-r[-1]['values']['net_if_all_anti_empty']/r[-1]['values']['net_change_per_equal_shape_scale'],
                            held_archive_Born=base,crossed_to_Born_net_g_point005=ratio))
    return dict(scope='Finite q crossed channel, anti-daughter equals stored daughter as an explicit shape comparator.',results=results)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
