"""Positive marginal balancing for sampled reversible quartet weights."""
import numpy as np
def balance(w,ids,targets,N,maxit=50,tol=1e-8):
    raw=w/w.sum();x=raw.copy();missing=False
    for j in range(4):
        c=np.bincount(ids[j],weights=x,minlength=N);missing|=bool(np.any((targets[j]>0)&(c==0)))
    if missing:return raw,1.,0,False
    for it in range(maxit):
        for j in range(4):
            c=np.bincount(ids[j],weights=x,minlength=N);ratio=np.divide(targets[j],c,out=np.zeros(N),where=c>0);x*=ratio[ids[j]]
        if it%5==4 or it==maxit-1:
            err=max(float(np.sum(abs(np.bincount(ids[j],weights=x,minlength=N)-targets[j]))) for j in range(4))
            if err<tol:return x,err,it+1,True
    return x,err,maxit,True
