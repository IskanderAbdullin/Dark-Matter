"""Finite-displacement check for a filtered fixed-energy tabulated FD bath."""
import argparse
import hashlib
import json
from pathlib import Path

from fd_portal_finite_bath_discrete_entropy import run as trajectories
from fd_portal_tabulated_exponential_filter import Filtered
from sm_dark_entropy_exchange_reference import StandardModel


def run(path):
    eos=StandardModel(path)
    a=trajectories(lambda T,u,n:Filtered(T,u,n,eos))
    return dict(rows=a['rows'],table_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                scope='Fixed-energy positive-filter trajectories compared against the unmodified timelike moment ODE, not expanding or abundance estimates.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    print(json.dumps(run(p.parse_args().eos),indent=2,allow_nan=False))
