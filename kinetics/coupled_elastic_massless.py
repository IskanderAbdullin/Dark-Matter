"""Coupled decay/portal/elastic pilot in the relativistic scalar limit.

Scalar m=0 throughout the kinetic integration (physical 86 eV is used for
the late abundance conversion). Contact quartic coupling corresponding to
a stated vacuum NR cross section. Exact conservation of represented energy
and particle counts; fixed elastic event geometry; no dynamical Hartree
mass, zero mode, number-changing reactions, or resummation. Weak-coupling
diagnostic only, not a relic prediction of the full contact theory.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from parent_angle_kernel import ParentAngle
from laine_portal_rates import PortalTable,make_models
from elastic_implicit_trial import ImplicitElastic
from quartic_consistency import lambda_from_sigma,shell_nodes,gap

_,TableResolved=make_models()

class Coupled(ParentAngle,TableResolved):
    def __init__(self,sigma=0.,**kw):
        super().__init__(**kw)
        assert self.scalar_ratio==0.
        self.sigma=sigma;self.lam=float(lambda_from_sigma(sigma))
        self.elastic=ImplicitElastic(self.edges,0.) if self.lam else None
        self.elastic_coefficient=self.lam**2/(2048*np.pi**5)*1.220890e19/(1.66*np.sqrt(self.g_bath)*self.parent_GeV)
    def dark(self,yd,Hf=1.,jac=False):
        if not self.lam:return super().dark(yd,Hf,jac)
        coef=self.elastic_coefficient*self.geom_x*Hf
        if jac:
            dy,J=super().dark(yd,Hf,True)
            dy[:self.n]+=coef*self.elastic.rates_jac(yd[:self.n],False)
            J[:self.n,:self.n]+=coef*self.elastic.rates_jac(yd[:self.n],True)
            return dy,J
        dy=super().dark(yd,Hf)
        dy[:self.n]+=coef*self.elastic.rates_jac(yd[:self.n],False)
        return dy
    def validate(self):
        r=super().validate()
        if self.elastic is not None:
            assert np.max(abs(self.energy_phi-self.elastic.K)/self.elastic.K)<1e-12
            r['elastic_jac']=self.elastic.validate_jac()
        return r

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--n',type=int,default=64);p.add_argument('--nparent',type=int,default=80)
    p.add_argument('--sigma',type=float,default=1e-7);p.add_argument('--qmin',type=float,default=1e-5)
    p.add_argument('--xend',type=float,default=.5);p.add_argument('--kappa',type=float,default=.00685)
    p.add_argument('--rtol',type=float,default=1e-5);p.add_argument('--wall',type=float,default=480.)
    p.add_argument('--output',required=True);p.add_argument('--max-step',type=float,default=.035)
    a=p.parse_args();t=PortalTable('.runtime/mesh_imPiR.dat')
    model=Coupled(table=t,n=a.n,nparent=a.nparent,kappa=a.kappa,sigma=a.sigma,
        scalar_ratio=0.,fermion_ratio=.9,qmin=a.qmin,g_fermion=4,
        parent_grid='balanced',fermion_grid='balanced',parent_floor=1e-18,record_scalar_history=True)
    r=model.evolve(.1,a.xend,a.rtol,'empty',True,86000.,a.wall,1e-11,'BDF',False,True,a.max_step)
    r.update(sigma_over_m_cm2_g=a.sigma,lambda_quartic=model.lam,
        approximation='massless scalar kinetics including exact discrete elastic and vacuum decay kernels; no dynamic gap or zero mode, real scalar number-changing omitted',
        rate_source_sha256=t.sha,kernel='parent-angle-elastic-massless')
    e=np.array(r['frozen_cell_spectrum']['edges'])
    for h,s in zip(r['history'],r['scalar_history_psd']):
        f=np.array(s['f_phi']);q,w=shell_nodes(e,f,8)
        mass=gap(q,w,1e-12*h['x'],model.lam)
        h.update(hartree_gap_over_peak=float(mass/h['peak_bin_q']),lambda_f_peak=float(model.lam*f.max()),
                 low_boundary_number_fraction=float(np.maximum(f,0)[:2]@model.vol[:2]/h['number']))
    Path(a.output).write_text(json.dumps([r],indent=2,allow_nan=False)+'\n')
    print('FINAL',json.dumps(dict(success=r['success'],message=r['message'],xend=r['xend'],sigma=a.sigma,wall=r['wall_seconds'],final=r['final'])),flush=True)
if __name__=='__main__':main()
