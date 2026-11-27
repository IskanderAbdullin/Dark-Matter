"""Finite joined evolution at explicitly chosen elastic packet quadrature orders."""
import argparse
import json
from unittest.mock import patch

from early_pole_archived_cold_screened_elastic_reference import Screened
from early_pole_cold_screened_split_reference import run as evolve
from early_pole_same_flavor_screened_pole_reference import SelfSoft
from early_pole_screened_elastic_channel_sum_reference import ElasticChannels
import early_pole_cold_screened_split_reference as impl


def run(projection,unlike_quad,same_quad,steps,duration,progress_every,checkpoint):
    if min(*unlike_quad,*same_quad)<1:
        raise ValueError('Positive packet quadrature is required.')
    class Unlike(Screened):
        def __init__(self,grid,wide,count,exposure,quadrature,cut,strength,*args,**kwargs):
            super().__init__(grid,wide,count,exposure,tuple(unlike_quad),cut,strength,*args,**kwargs)
    class Both(ElasticChannels):
        def __init__(self,unlike,grid,same_soft_strength):
            self.M=unlike.M
            self.unlike=unlike
            self.soft=SelfSoft(self.M,grid,tuple(same_quad),.08)
            self.soft.strength=same_soft_strength
    # Do not use a checkpoint through this wrapper: graph controls are
    # explicit in the output, while base checkpoint metadata describe
    # the default underlying graph. Interrupted trials may be rerun.
    if checkpoint is not None:
        raise ValueError('Checkpointing is not supported for this quadrature reference.')
    with patch.object(impl,'Screened',Unlike),patch.object(impl,'ElasticChannels',Both):
        z=evolve('fine',duration,steps,True,projection=projection,screening_update='full',
                 same_soft_strength=1.,compton_strength=1.,conversion_strength=1.,snapshots=True,
                 elastic_energy_budget=5e-15,progress_every=progress_every)
    z.update(unlike_packet_quadrature=list(unlike_quad),same_packet_quadrature=list(same_quad))
    return z


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    p.add_argument('--unlike-quad',nargs=3,type=int,default=[2,1,2])
    p.add_argument('--same-quad',nargs=3,type=int,default=[2,1,2])
    p.add_argument('--steps',type=int,default=16)
    p.add_argument('--duration',type=float,default=3.3333333333333335e-7)
    p.add_argument('--progress-every',type=int,default=4)
    p.add_argument('--checkpoint',default=None)
    a=p.parse_args()
    print(json.dumps(run(**vars(a)),indent=2,allow_nan=False))
