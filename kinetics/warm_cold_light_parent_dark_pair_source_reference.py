"""Tagged positive heavy creation from frozen dark pair channels."""

import argparse
import json

import numpy as np
from pathlib import Path
from scipy.special import roots_legendre
from archived_scalar_pair_reaction_reference import amplitude, gross
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from three_species_decay_onsager_reference import M,m,mp,y,width
from warm_light_occupation_stage_audit_reference import read, digest
from warm_near_rest_light_inverse_decay_reference import E0, window

# Fix tagged H at the lab energy needed to decay to a light particle at rest.
def source(E, occ, gate, d, radial, pol, az, panel, npe=18):
    bmax=2*gate[-1]-E
    if bmax<=M:return 0.
    p=np.sqrt((E-M)*(E+M))
    x,w=roots_legendre(radial)
    edge=np.linspace(M,bmax,npe+1)
    de=np.diff(edge)/2
    f=((edge[:-1]+edge[1:])[:,None]/2+de[:,None]*x).ravel()
    v=(de[:,None]*w).ravel()
    k=np.sqrt((f-M)*(f+M))
    cin,wi=roots_legendre(pol)
    c=np.cos((np.arange(az)+.5)*np.pi/az)
    xo,wo=roots_legendre(panel)
    total=(E+f[:,None])+np.zeros_like(cin)
    S=2*M*M+2*(E*f[:,None]-p*k[:,None]*cin)
    p0=np.sqrt(np.maximum(0.,S/4-M*M))
    p1=np.sqrt(S/4-d*d)
    boost=np.sqrt(np.maximum(total*total-S,0.))/np.sqrt(S)
    dc=boost*p1; mid=total/2.
    u=(gate-mid[:,:,None])/dc[:,:,None]
    gates=np.sort(np.clip(np.concatenate([u,-u,-np.ones_like(u[:,:,:1]),np.ones_like(u[:,:,:1])],axis=2),-1.,1.),axis=2)
    dzo=np.diff(gates,axis=2)/2
    zo=(gates[:,:,1:]+gates[:,:,:-1])/2
    z=(zo[:,:,:,None]+dzo[:,:,:,None]*xo).reshape(*S.shape,-1)
    vw=(dzo[:,:,:,None]*wo).reshape(*S.shape,-1)
    x0=np.clip((E-f[:,None])/(2*np.maximum(boost*p0,1e-280)),-1.,1.)
    mu=np.clip(x0[:,:,None,None]*z[:,:,:,None]+np.sqrt(1-x0*x0)[:,:,None,None]*np.sqrt(1-z*z)[:,:,:,None]*c,-1.,1.)
    if d==mp:
        mat=amplitude(S[:,:,None,None],mu,light=M,heavy=m,scalar=mp).mean(axis=3)
        symmetry=1.
    else:
        tt=M*M+m*m-S[:,:,None,None]/2+2*(p0*p1)[:,:,None,None]*mu
        assert np.max(tt)<mp*mp
        mat=(y**4*(((M+m)**2-tt)/(tt-mp*mp))**2).mean(axis=3)
        symmetry=2.
    eoutp=mid[:,:,None]+dc[:,:,None]*z
    eoutm=mid[:,:,None]-dc[:,:,None]*z
    wf=np.sum(vw*mat*occ(eoutp)*occ(eoutm),axis=2)
    pref=p1/(64*np.pi*E*f[:,None]*np.sqrt(S))
    # H partner has two spins; angular incoming average over cosine.
    return symmetry*np.sum((v*k*f/np.pi**2)[:,None]*(wi/2)*pref*wf)


def check():
    # Integrate the tagged rate over one H charge and compare to the
    # independent scalar-first 2-body event integral, which counts H pairs.
    p=np.linspace(0.,.24,17)
    B=lambda E:np.where(E<np.hypot(mp,p[-1]),.03*np.exp(-E/.08),0.)
    z,w=roots_legendre(35)
    emax=2*np.hypot(p[-1],mp)-M
    pmax=np.sqrt(emax*emax-M*M)
    pp=pmax*(z+1)/2
    Es=np.hypot(pp,M)
    total=np.array([source(E,B,np.hypot(p,mp),mp,6,15,12,4,16) for E in Es])
    tag=np.sum(w*pp*pp*total)*pmax/(2*np.pi**2)
    ref=gross(p,mp,M,B,lambda E:np.zeros_like(E),True,6,15,12,4,exchange_mass=m)[1][1]
    defect=abs(tag/ref-1)
    assert defect<.006,(tag,ref,defect)
    return dict(tagged_integrated_rate=tag, scalar_first_independent_event_rate=ref,
                relative_normalization_defect=defect)


def run(D):
    results=[]
    orders=((3,8,8,3,12),(4,12,10,4,16),(6,18,14,5,20),(8,24,16,6,26))
    for q in (16,32):
        A=read(D/f'warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz')
        for j in (0,3):
            row=A[j]; S=row["snapshot"]
            assert row["snapshot_sha256"]==digest(S)
            B,bg,F,fg=snapshot_functions(S)
            t=S["T_visible_GeV"]
            _,_,_,lo,hi=window(np.asarray(S["lab_momentum_edges_GeV"])[1])
            center=E0+m
            e_all=[lo+m,center,hi+np.hypot(m,S["lab_momentum_edges_GeV"][1])]
            a,b=np.asarray(S["entropy_variable_fermion"])
            maxlight=1/(1+np.exp(min(np.min(a),np.min(a+b))))
            data=[]
            for order in orders:
                vals=np.asarray([[source(E,fn,g,d,*order) for fn,g,d in ((B,bg,mp),(F,fg,m))]
                                 for E in e_all])
                stat=vals*np.asarray(e_all)[:,None]/(width*M*(1-maxlight))
                data.append(dict(order=order,one_way_positive_rates_GeV=vals.tolist(),
                                 divided_by_lower_decay_damping=stat.tolist()))
            if q==16:
                assert max(bg[-1],fg[-1])*2<min(e_all)+M
                assert np.all(vals==0.)
            else:
                assert vals.min()>0.
                last=np.asarray(data[-2]["one_way_positive_rates_GeV"])
                assert np.max(abs(vals/last-1))<.003
            results.append(dict(q=q,j=j,T_GeV=t,snapshot_sha256=row["snapshot_sha256"],
                                tagged_parent_energies_GeV=e_all,orders=data))
    return dict(scope="Frozen represented internal H-pair gain at energies associated "
                      "with near-rest light recoil, with both H Pauli factors omitted. "
                      "Numerical sampled rate, not a continuous bound or a full evolving source.",
                normalization_check=check(),results=results)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(p.parse_args().tapes),indent=2,allow_nan=False))
