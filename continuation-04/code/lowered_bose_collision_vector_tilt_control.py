"""Synthetic vector-moment tilt control; not an actual halo C22 solve."""
import numpy as np,json
rng=np.random.default_rng(25);n=40;k=3;N=1000
f=np.exp(rng.uniform(-5,2,n));o=rng.normal(size=(n,k));ids=rng.integers(0,n,(4,N));w=np.exp(rng.uniform(-2,2,N));i,j,l,m=ids;a=f[i]*f[j]*(1+f[l])*(1+f[m]);b=(1+f[i])*(1+f[j])*f[l]*f[m];nu=o[l]+o[m]-o[i]-o[j]
def eval(x):
    u=nu@x/2;A=w*a*np.exp(-u);B=w*b*np.exp(u);J=A-B
    return 2*(A+B).sum(),-nu.T@J,nu.T@((A+B)[:,None]*nu)/2,J
x=np.zeros(k)
for it in range(30):
    val,g,H,J=eval(x)
    if np.linalg.norm(g)<1e-9:break
    step=np.linalg.solve(H,g);t=1.
    while eval(x-t*step)[0]>val+1e-12:t/=2
    x-=t*step
v,g,H,J=eval(x);raw=w*(a-b);S=J@np.log(a/b);St=J@(np.log(a/b)-nu@x)
print(json.dumps(dict(events=N,constraints=k,iterations=it,raw_moment_norm=float(np.linalg.norm(nu.T@raw)),final_moment_norm=float(np.linalg.norm(g)),entropy=float(S),tilted_entropy=float(St),entropy_difference=float(S-St),min_hessian_eigenvalue=float(np.linalg.eigvalsh(H)[0]),relative_flux_change=float(abs(J-raw).sum()/abs(raw).sum()))))
