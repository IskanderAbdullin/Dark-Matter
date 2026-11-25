"""Cross-halo consistency screen for a subthreshold virtual mediator.

The same small mixing that suppresses light-light elasticity enhances
LL->HH relative to LL->LL once kinematics open. Test all u/a grid points,
not only the dwarf-optimal ones. Cluster collision optical depth is a
diagnostic for the stipulated dilute two-stream geometry, not a fit to
cluster observations or a transferred elastic SIDM confidence limit.
"""
import argparse,json
from pathlib import Path
import numpy as np
from subthreshold_transport_tradeoff import mb_integral
from soft_mixing_activation import frac_B
from two_state_halo_equilibrium import pl
from direct_activation_consistency import rho_tot,radiation_T4
from late_activation_overlap import critical_density
from resonant_mediator_inverse import M,C,SEC_GYR
from quartic_consistency import lambda_from_sigma

VREF=np.sqrt(2)*44.1e5
SCALE=M*M*(VREF/C)**2
RHO=critical_density(VREF)

def cross_sections(sig0,u,D,delta,a,vkms,x=1.):
    """Exact CM phase space, dominant s-channel only.
    All possible incoming/outgoing L,H species. v is 2p_i/mL (NR convention).
    sigma/m is normalized to reference L-L threshold amplitude at x=1.
    """
    u=np.asarray(u)
    v=vkms*1e5/C
    shift=lambda x:-2*delta**2*x*x/(np.hypot(D,2*delta*x)+D)
    ml2=M*M+shift(x);mh2=ml2+np.hypot(D,2*delta*x)
    ms2=4*(M*M+shift(1))-a*SCALE
    ms=[np.sqrt(ml2),np.sqrt(mh2)]
    mix=frac_B(x,delta,D);f=[mix,1-mix]
    p=np.sqrt(ml2)*v/2
    pairs=[(0,0),(0,1),(1,1)]
    result={}
    for i,j in pairs:
        mi=ms[i]+ms[j]
        dkin=p*p/(np.hypot(ms[i],p)+ms[i])+p*p/(np.hypot(ms[j],p)+ms[j])
        root=mi+dkin;ss=root*root
        pole=(a*SCALE/((mi*mi-ms2)+dkin*(2*mi+dkin)))**2
        norm=sig0*(4*M*M/ss)*pole*f[i]*f[j]/u**4
        for k,l in pairs:
            mf=ms[k]+ms[l];dif=(mi-mf)+dkin
            kf2=np.maximum(dif*(root+mf)*(ss-(ms[k]-ms[l])**2)/(4*ss),0.)
            phase=np.where(dif>=0,np.sqrt(kf2)/p,0.)
            result[f'{i}{j}->{k}{l}']=norm*f[k]*f[l]*phase*(2/(1+(k==l)))
    return result

def scan_one(s,aa,QQ):
    p=s['params'];q=s['screen'];A4=q.get('A4_Gyr',0)*SEC_GYR
    if A4<=0 or 'cosmic_contact_factor' not in q:return dict(feasible=False)
    A2=q['A2_Gyr']*SEC_GYR;A3=q['A3_Gyr']*SEC_GYR
    u=np.geomspace(1e-8,.49,3200)
    delta=p['delta'];D=delta*(1-2*u)/np.sqrt(u*(1-u))
    d=delta/(SCALE*np.sqrt(u*(1-u)))
    R=(1-u)/u
    B=(np.exp(-1e-6)/float(pl(1.5,1e-6)))**2
    counts={k:0 for k in ['loop_pass','loop_and_cosmic','loop_cosmic_conversion','all_including_cluster_1','all_including_cluster_3000_1','all_including_cluster_4500_1']}
    best=None;best_conv_at_cluster=None
    for a,Q in zip(aa,QQ):
        sig0=1/(RHO*VREF*Q*A4)
        g2=float(lambda_from_sigma(sig0))*a*SCALE/u**2
        lam=3*g2/(4*M*M-a*SCALE)
        loop=(g2/(16*np.pi**2)+lam*M*M/(32*np.pi**2))<=D
        loop &= g2/(32*np.pi**2)<=a*SCALE
        loop &= (lam<1.)&(g2/(4*M*M-a*SCALE)<1.)
        counts['loop_pass']+=int(sum(loop))
        cosmic=q['cosmic_contact_factor']/Q
        if cosmic>.01:continue
        counts['loop_and_cosmic']+=int(sum(loop))
        mb=B*(R*mb_integral(d,a)*A3/A4+R**2*mb_integral(2*d,a)*A2/A4)/Q
        good=loop&(mb<=.1)
        counts['loop_cosmic_conversion']+=int(sum(good))
        # Cluster at 1000 km/s, require total cross section <1 cm2/g
        # as generous *diagnostic*, not observational likelihood.
        cv={}
        for v in [1000,3000,4500]:
            mat=cross_sections(sig0,u,D,delta,a,v)
            cv[v]=sum(mat['00->'+f] for f in ['00','01','11'])
        counts['all_including_cluster_1']+=int(sum(good&(cv[1000]<1.)))
        counts['all_including_cluster_3000_1']+=int(sum(good&(cv[3000]<1.)))
        counts['all_including_cluster_4500_1']+=int(sum(good&(cv[4500]<1.)))
        if np.any(good):
            idx=np.where(good)[0][np.argmin(cv[4500][good])]
            v=dict(a=float(a),u=float(u[idx]),D=float(D[idx]),sigma0=float(sig0),gap_over_T=float(d[idx]),
                   MB_initial_H_source_lower=float(mb[idx]),cosmic_proxy=float(cosmic),
                   **{f'LL_total_sigma_m_{k}kms':float(vv[idx]) for k,vv in cv.items()})
            if best is None or v['LL_total_sigma_m_4500kms']<best['LL_total_sigma_m_4500kms']:best=v
        ok=loop&(cv[1000]<1.)
        if np.any(ok):
            idx=np.where(ok)[0][np.argmin(mb[ok])]
            v=dict(a=float(a),u=float(u[idx]),MB_initial_H_source_lower=float(mb[idx]),sigma_1000=float(cv[1000][idx]))
            if best_conv_at_cluster is None or v['MB_initial_H_source_lower']<best_conv_at_cluster['MB_initial_H_source_lower']:best_conv_at_cluster=v
    return dict(feasible=True,counts=counts,min_cluster_sigma_among_dwarf_candidates=best,
                min_MB_conversion_among_cluster_thin=best_conv_at_cluster)

def report_candidate(p):
    out=[]
    for x in [.1,.3,.7,1.,1.2]:
        for v in [30,60,100,150,250,500,1000,2000,3000,4500]:
            s=cross_sections(p['sigma0_for_one_efold'],p['u'],p['D'],p['delta'],p['a'],v,x)
            s={k:float(v) for k,v in s.items()}
            totals={k:sum(s[k+'->'+f] for f in ['00','01','11']) for k in ['00','01','11']}
            out.append(dict(x=x,v_kms=v,channels=s,total_by_incoming=totals,
                       minimum_total_sigma_m=min(totals.values()),tau_0p3g_cm2_min=.3*min(totals.values())))
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--tradeoff',required=True)
    ap.add_argument('--scan',nargs='+',required=True);ap.add_argument('--refined',nargs='*',default=[])
    ap.add_argument('--output',required=True);a=ap.parse_args()
    table=json.loads(Path(a.tradeoff).read_text());aa=np.array(table['aa']);QQ=np.array(table['Q96'])
    scans=sum([json.loads(Path(p).read_text()) for p in a.scan],[])
    rows=[]
    for s in scans:
        v=scan_one(s,aa,QQ)
        rows.append(dict(params=s['params'],scheme=s['scheme'],**v))
        if v['feasible'] and v['counts']['loop_cosmic_conversion']:
            print(s['params']['z_on'],s['params']['dneff'],s['scheme'],v['counts'],v['min_cluster_sigma_among_dwarf_candidates'],flush=True)
    candidates=[]
    for f in a.refined:
        for b in json.loads(Path(f).read_text()):
            candidates.append(dict(index=b['index'],candidate=b['runs'][0]['candidate'],rows=report_candidate(b['runs'][0])))
    counts={k:sum(v.get('counts',{}).get(k,0) for v in rows) for k in rows[0]['counts']}
    print('TOTAL',counts,flush=True)
    Path(a.output).write_text(json.dumps(dict(rows=rows,candidates=candidates,total_counts=counts,
      scope='Discrete grid, small-u roll shortlist, generous full-field cluster optical-depth screen; no inference of observational confidence or universally excluded EFT.'),indent=2,allow_nan=False)+'\n')
