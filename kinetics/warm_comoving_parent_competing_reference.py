"""Selected competing dark reactions for the comoving parent block."""
import numpy as np

from compact_elastic_panel_reference import compile_panels,compact_elastic
from joint_below_pole_elastic_linear_modes_reference import setup
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from scalar_pair_expanding_continuation_reference import MASS,moments,reconstruct
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_late_finite_parent_prompt_reference import FiniteParent
from warm_pair_common_fermion_chart_reference import rate


class ElasticMixin:
    def load_elastic(self,k,T):
        self.panels=[]
        lim=sum(np.hypot(c.p[-1],c.mass) for c in k.c)
        setup(T,k,T,(0.,0.),3,2,2,8,lim/T,.82,cache=self.panels,cache_only=True)
        self.elastic=compile_panels(self.panels,self.n)

    def centered(self,z,ref,R0):
        a,da=self.collision(ref+z)
        b,db=compact_elastic(z,self.t,self.eta,self.H,self.elastic,self.d,True)
        return a+b,lambda v: da(v)+db(v)


class FiniteParentElastic(ElasticMixin,FiniteParent):
    def __init__(self,k,H,T):
        super().__init__(k,H)
        self.load_elastic(k,T)


class StationaryParentElastic(ElasticMixin,LocalAbsolutePrompt):
    def __init__(self,k,H,T):
        super().__init__(k,H,[],True,stream_local=True)
        self.load_elastic(k,T)



def midpoint_pair_kick(p,theta,a,dt,exterior):
    def rhs(th):
        f,_,g,_=snapshot_functions(dict(lab_momentum_edges_GeV=p,
             entropy_variable_scalar=th[0],entropy_variable_fermion=th[1]))
        return rate(p,g,f,(3,7,6,2))
    assert a==exterior.a
    start=moments(p,theta)
    first=rhs(theta)[0]
    mid=reconstruct(p,start+.5*dt*first[:,:,1:-1])
    net,_,_,defect,_=rhs(mid)
    q=dt*net
    result=reconstruct(p,start+q[:,:,1:-1])
    outside=a**3*q[:,0,[0,-1]].sum(axis=1)
    change=(a**3*q[:,0].sum(axis=1))
    for s,m in enumerate(MASS):
        for side,i in enumerate((0,-1)):
            exterior.add(a,m,s,side,p[i],a**3*q[s,0,i],a**4*q[s,1,i])
    return result,change,outside,float(defect)
