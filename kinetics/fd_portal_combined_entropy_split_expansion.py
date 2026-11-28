"""Short expanding FD benchmark with both direct timelike and elastic channels."""
import argparse
import hashlib
import json
from pathlib import Path

from fd_portal_combined_tabulated_bath_entropy import Combined, elastic_heat
from fd_portal_tabulated_entropy_split_expansion import Stage, run as verify


class CombinedExpansion(Stage):
    def closed_bath(self,T,u):
        return Combined(T,u,self.nodes,self.eos)

    def collision(self,Ts,Td,eta):
        f,s=super().collision(Ts,Td,eta)
        mesh=(18,29.) if self.nodes==112 else (25,36.)
        q=elastic_heat(Ts,Td,eta,mesh)
        t=q*(1./Ts-1./Td)
        if t<0.:
            raise ArithmeticError('Invalid elastic source sign in independent nonsplit flow')
        f[1]+=q
        return f,s+t


def run(path,timelike=None):
    z=verify(path,stage=CombinedExpansion)
    z['scope']='Short homogeneous expansion with direct timelike and elastic vacuum charged-lepton FD kernels. Not a full history or final density.'
    if timelike is not None:
        data=timelike.read_bytes()
        b=json.loads(data)
        if z['table_input_sha256']!=b['table_input_sha256'] or len(z['rows'])!=len(b['rows']):
            raise ArithmeticError('Need comparable exact archived bath EOS and trajectories')
        delta=[]
        for u,v in zip(z['rows'],b['rows'],strict=True):
            if u['initial_state']!=v['initial_state']:
                raise ArithmeticError('Different initial thermal states')
            p,q=u['levels'][-1],v['levels'][-1]
            same=bool(u['horizon']==v['horizon'])
            shift=[a-b for a,b in zip(p['state'],q['state'],strict=True)] if same else None
            delta.append(dict(T_SM=u['initial_state'][0],same_duration=same,
                              state_difference_visible_T_eta_log_dark_T=shift,
                              relative_comoving_entropy_difference=(p['relative_comoving_entropy_added']-q['relative_comoving_entropy_added']) if same else None))
        z['timelike_comparison']=dict(reference_sha256=hashlib.sha256(data).hexdigest(),samples=delta)
    return z


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('--timelike',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.timelike),indent=2,allow_nan=False))
