"""Verify bitwise restoration at an accepted SciPy 1.17 DOP853 step."""
import hashlib
import inspect
import sys
import json
from pathlib import Path
import tempfile

import numpy as np
import scipy
from scipy.integrate import DOP853, solve_ivp


def encode(d):
    return json.dumps(d,sort_keys=True,separators=(',',':'),allow_nan=False).encode()


def fingerprint():
    return hashlib.sha256(inspect.getsource(DOP853._step_impl).encode()+
                          inspect.getsource(DOP853._dense_output_impl).encode()).hexdigest()


def evolve(rhs,span,initial,times,max_step,rtol,atol,tag,path,start=None,stop=None,save_model_state=None,load_model_state=None,save_every=None):
    if save_every is not None and save_every<1:
        raise ValueError('Accepted-step save cadence must be positive')
    a,b=span
    times=np.asarray(times,dtype=float)
    if (times.ndim!=1 or not len(times) or not np.all(np.isfinite(times))
            or not np.all(np.diff(times)>0.) or not a<=times[0]<=times[-1]<=b):
        raise ValueError('Use forward integration with strictly ordered report times')
    config=dict(span=list(span),initial=np.asarray(initial).tolist(),times=np.asarray(times).tolist(),
                max_step=max_step,rtol=rtol,atol=np.asarray(atol).tolist(),
                scipy=scipy.__version__,numpy=np.__version__,python=sys.version.split()[0],solver_sources=fingerprint(),caller=tag,external_state=bool(save_model_state or load_model_state),
                schema=3 if (save_model_state or load_model_state) else 2)
    steps=[];history=[];widths=[];i=0
    if start is not None:
        package=json.loads(start.read_text())
        saved=package['body']
        if hashlib.sha256(encode(saved)).hexdigest()!=package['sha256'] or saved['spec']!=config:
            raise ValueError('Damaged or mismatched accepted-step checkpoint')
        if load_model_state is not None:
            load_model_state(saved['model_state'])
        t,state,f,H,nev=(saved[k] for k in ('t','y','f','next_step','calls'))
        state,f=np.asarray(state),np.asarray(f)
        steps=list(saved['steps'])
        history=list(saved['samples'])
        widths=list(saved['widths'])
        i=len(history)
        if sum(widths)!=i or any(w<=0 for w in widths):
            raise ValueError('Invalid dense-evaluation block grouping')
        if not a<t<b or H<=0. or not np.isfinite(H) or not np.all(np.isfinite(state)&np.isfinite(f)):
            raise ValueError('Bad step state')
        if state.shape!=np.asarray(initial).shape or f.shape!=state.shape or steps[-1]!=t:
            raise ValueError('Bad step dimensions or inconsistent time')
        flag=[True]
        def wrapped(T,Y):
            if flag[0]:
                flag[0]=False
                if T!=t or not np.array_equal(Y,state):
                    raise ValueError('SciPy constructor called unexpected initial point')
                return f
            return rhs(T,Y)
        s=DOP853(wrapped,t,state,b,max_step=max_step,rtol=rtol,atol=atol,
                 first_step=min(H,max_step,b-t))
        s.h_abs=H
        s.nfev=nev
    else:
        s=DOP853(rhs,a,initial,b,max_step=max_step,rtol=rtol,atol=atol)
    while s.status=='running':
        s.step()
        if s.status=='failed':raise ArithmeticError('Step rejection hit the minimum')
        steps.append(float(s.t))
        j=np.searchsorted(times,s.t,'right')
        if j>i:
            D=s.dense_output()
            history.extend(D(times[i:j]).T.tolist())
            widths.append(int(j-i))
        i=j
        suspended=stop is not None and len(steps)==stop
        periodic=save_every is not None and len(steps)%save_every==0
        if (suspended or periodic) and s.status=='running':
            data=dict(spec=config,t=float(s.t),y=s.y.tolist(),f=s.f.tolist(),next_step=float(s.h_abs),
                      calls=s.nfev,samples=history,steps=steps,widths=widths,
                      model_state=None if save_model_state is None else save_model_state())
            temp=path.with_suffix('.tmp')
            import os
            with temp.open('wb') as f:
                f.write(encode(dict(body=data,sha256=hashlib.sha256(encode(data)).hexdigest())))
                f.flush()
                os.fsync(f.fileno())
            temp.replace(path)
            parent_fd=os.open(path.parent,os.O_RDONLY)
            try:
                os.fsync(parent_fd)
            finally:
                os.close(parent_fd)
            if suspended:
                return None
    return dict(steps=steps,samples=history,calls=s.nfev,widths=widths)


def dense_matrix(result):
    y=np.asarray(result['samples'])
    ends=np.r_[0,np.cumsum(result['widths'])]
    return np.hstack([y[a:b].T for a,b in zip(ends[:-1],ends[1:],strict=True)])


def rhs(t,y):
    a,b,c=y
    return np.array([-.3*a+np.sin(51*t)*b, .05*a-.9*b+.2*c, .3*a-.1*c+np.cos(25*t)*b])


def parameters():
    t=np.linspace(0,1.3,89)
    return [0.,1.3],[1.1,.15,.6],t,.034,2e-10,3e-13,'fixture'


def run():
    with tempfile.TemporaryDirectory() as temp:
        f=Path(temp)/'checkpoint.json'
        c=parameters()
        t=c[2]
        a=evolve(rhs,*c,f)
        ref=solve_ivp(rhs,c[0],c[1],method='DOP853',t_eval=t,max_step=c[3],rtol=c[4],atol=c[5])
        np.testing.assert_array_equal(dense_matrix(a),ref.y)
        if dense_matrix(a).strides!=ref.y.strides:
            raise ArithmeticError('Restored solve_ivp layout differs')
        checks=[]
        for j in (1,3,8,29,45):
            if evolve(rhs,*c,f,stop=j) is not None:
                raise ValueError('Missing test suspension')
            b=evolve(rhs,*c,f,start=f)
            if a!=b:raise ArithmeticError(('Nonbitwise adaptive continuation',j))
            checks.append(j)
        live=evolve(rhs,*c,f,save_every=9)
        if live!=a or evolve(rhs,*c,f,start=f)!=a:
            raise ArithmeticError('Periodic writes or their restart changed the result')
        bad=json.loads(f.read_text())
        bad['body']['next_step']*=1.01
        f.write_bytes(encode(bad))
        try:
            evolve(rhs,*c,f,start=f)
        except ValueError:
            pass
        else:
            raise ArithmeticError('Corrupt checkpoint accepted')
        return dict(scipy=scipy.__version__,adaptive_steps=len(a['steps']),dense_samples=len(t),
                    nfev=a['calls'],exact_restarted_steps=checks,bitwise_solver_and_resume=True,
                    corrupted_checkpoint_rejected=True,periodic_replay_is_bitwise=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    x=p.add_mutually_exclusive_group()
    x.add_argument('--write-checkpoint',type=Path)
    x.add_argument('--restore-checkpoint',type=Path)
    a=p.parse_args()
    if a.write_checkpoint:
        evolve(rhs,*parameters(),a.write_checkpoint,stop=8)
    elif a.restore_checkpoint:
        r=evolve(rhs,*parameters(),a.restore_checkpoint,start=a.restore_checkpoint)
        print(hashlib.sha256(encode(r)).hexdigest())
    else:
        print(json.dumps(run(),indent=2,allow_nan=False))
