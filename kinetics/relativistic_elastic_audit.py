"""Deterministic relativistic on-shell finite-T event audit; NOT an integrated Boltzmann kernel."""
import json, math, random
def dot(a,b): return a[0]*b[0]-sum(x*y for x,y in zip(a[1:],b[1:]))
def boost(v,b):
 b2=sum(x*x for x in b)
 if b2<1e-28:return v
 g=1/math.sqrt(1-b2);bp=sum(x*y for x,y in zip(b,v[1:]));c=(g-1)*bp/b2+g*v[0]
 return (g*(v[0]+bp),*(v[i+1]+c*b[i] for i in range(3)))
def fd(E,T,mu=0):
 x=(E-mu)/T
 return 1/(1+math.exp(x)) if x<700 else 0.
def event(M,p,k,theta,phi):
 E=math.hypot(M,p);P=(E,0.,0.,p);K=(k,k*math.sin(theta),0.,k*math.cos(theta))
 total=tuple(a+b for a,b in zip(P,K));s=dot(total,total)
 beta=tuple(x/total[0] for x in total[1:]);star=boost(P,tuple(-x for x in beta))
 pmag=math.sqrt(max(0,star[0]**2-M*M))
 n=[x/math.sqrt(sum(t*t for t in star[1:])) for x in star[1:]] if pmag>1e-14 else [0,0,1]
 axis=[0.,1.,0.] if abs(n[1])<.9 else [1.,0.,0.]
 proj=sum(x*y for x,y in zip(n,axis));a=[axis[i]-proj*n[i] for i in range(3)]
 an=math.sqrt(sum(x*x for x in a));a=[x/an for x in a]
 c=math.cos(phi);t=math.sin(phi);direction=[c*n[i]+t*a[i] for i in range(3)]
 P2=boost((math.hypot(M,pmag),*(pmag*x for x in direction)),beta)
 K2=tuple(total[i]-P2[i] for i in range(4))
 return P,K,P2,K2,s
def run(n=4000,seed=21112026):
 rng=random.Random(seed);shell=balance=eq=0.;entropy_min=1e9;entropy_weighted_heat=0.;count=0
 for _ in range(n):
  M=10**rng.uniform(-1,1);T=10**rng.uniform(-1,1);p=rng.uniform(.05,5)*T;k=rng.uniform(.05,5)*T
  P,K,P2,K2,s=event(M,p,k,rng.uniform(.02,3.12),rng.uniform(.02,3.12))
  shell=max(shell,max(abs(dot(P,P)-M*M),abs(dot(P2,P2)-M*M),abs(dot(K,K)),abs(dot(K2,K2)))/max(1,M*M,s))
  balance=max(balance,max(abs(P[i]+K[i]-P2[i]-K2[i]) for i in range(4))/max(1,P[0]+K[0]))
  muc=rng.uniform(-.2,.2)*T;mue=rng.uniform(-.2,.2)*T
  f,g,f2,g2=fd(P[0],T,muc),fd(K[0],T,mue),fd(P2[0],T,muc),fd(K2[0],T,mue)
  F=f*g*(1-f2)*(1-g2);R=f2*g2*(1-f)*(1-g)
  eq=max(eq,abs(F-R)/max(F,R,1e-290))
  Tc=T*rng.uniform(.7,1.3)
  f,g,f2,g2=fd(P[0],Tc,muc),fd(K[0],T,mue),fd(P2[0],Tc,muc),fd(K2[0],T,mue)
  F=f*g*(1-f2)*(1-g2);R=f2*g2*(1-f)*(1-g)
  if min(F,R)>1e-290:
   entropy_min=min(entropy_min,(F-R)*math.log(F/R))
   entropy_weighted_heat+=(P2[0]-P[0])*(F-R)*(1/Tc-1/T);count+=1
 assert shell<1e-12 and balance<1e-12 and eq<1e-11 and entropy_min>=-1e-15 and entropy_weighted_heat>=-1e-13
 return dict(seed=seed,n=n,max_relative_shell=shell,max_relative_4momentum=balance,max_relative_equilibrium_gain_loss=eq,min_event_entropy=entropy_min,weighted_heat_entropy_sign=entropy_weighted_heat,nonequilibrium_samples=count,scope="random event checks, NOT integrated collision operator")
if __name__=="__main__":print(json.dumps(run(),indent=2))
