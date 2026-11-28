"""Sum disjoint reciprocal charged packet channels on one variational state."""

from early_pole_same_flavor_screened_pole_reference import SelfSoft


class ElasticChannels:
    def __init__(self,unlike,grid,same_soft_strength):
        if same_soft_strength<0.:
            raise ValueError("negative same-flavor rate")
        self.M=unlike.M
        self.unlike=unlike
        self.soft=SelfSoft(self.M,grid,(1,1,2),.08)
        self.soft.strength=same_soft_strength

    def linear_graphs(self):
        if self.unlike.strength:
            yield self.unlike.graph,self.M.r[2],self.unlike.factor*self.unlike.strength
        if self.soft.strength:
            for k,graph in enumerate(self.soft.graph,1):
                yield graph,self.M.r[k],self.soft.fac*self.soft.strength

    def collision(self,G,F,E):
        D,S,W=self.soft.collision(G,F,E)
        H=0.
        if self.unlike.strength:
            add,s,w,H=self.unlike.collision(G,F,E)
            D+=add
            S.append(s)
            W.append(w)
        return D,sum(S),sum(W),H
