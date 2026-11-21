/* Preliminary SPARC fixed-m=86eV halo comparison, 2026-11-21.
Node >=18, no npm dependencies:
  node sparc_pilot_compare.js /path/to/MassModels_Lelli2016c.mrt > pilot-results.json
SOURCE: https://astroweb.cwru.edu/SPARC/MassModels_Lelli2016c.mrt
NOT a hierarchical fit; fixed M/L 0.5,0.7, diagonal random errors only.
*/
"use strict";
const fs=require("node:fs");
const ZETA=2.612375348685488, K86=41.87085711141446, PI=Math.PI;
function readData(filename) {
 const groups={};
 for(const line of fs.readFileSync(filename,"utf8").split(/\r?\n/)) {
  const a=line.trim().split(/\s+/);
  if(a.length!==10||!/^[A-Za-z]/.test(a[0]))continue;
  const n=a.slice(1).map(Number);
  if(n.some(x=>!Number.isFinite(x)))continue;
  (groups[a[0]]??=[]).push({r:n[1],v:n[2],dv:n[3],vg:n[4],vd:n[5],vb:n[6]});
 }
 for(const g of Object.values(groups))g.sort((a,b)=>a.r-b.r);
 const count=Object.values(groups).reduce((s,g)=>s+g.length,0);
 if(Object.keys(groups).length!==175||count!==3391)
  throw Error("Expected exact SPARC Lelli+2016c 175 galaxies and 3391 rows");
 return groups;
}
function dm2(o){return o.v**2-o.vg*Math.abs(o.vg)-.5*o.vd**2-.7*o.vb**2;}
function choose(groups) {
 return Object.entries(groups).filter(([id,g])=>{
  const f=g.map(x=>dm2(x)/x.v**2).sort((a,b)=>a-b);
  return g.length>=15&&f[Math.floor(f.length/2)]>.7&&g.every(x=>dm2(x)>=0);
 }).sort((a,b)=>a[0].localeCompare(b[0])).slice(0,20);
}
function li(t) {
 if(t<=0)return ZETA;
 if(t<.35)return ZETA-2*Math.sqrt(PI*t)+1.4603545088095868*t-.1039431124886773*t*t
   +.0042475336483055*t*t*t+.0003548720324104*t*t*t*t;
 let sum=0;
 for(let k=1;k<250;k++){const add=Math.exp(-k*t)/k**1.5;sum+=add;if(add<1e-13)break;}
 return sum;
}
function invert(n) {
 if(n>=ZETA)return 0;
 let l=0,r=60;
 for(let i=0;i<50;i++){const m=(l+r)/2;if(li(m)>n)l=m;else r=m;}
 return(l+r)/2;
}
const S0=Math.log(1e-4), S_MAX=Math.log(1000), N=440, H=(S_MAX-S0)/N;
function profile(n) {
 const t=invert(n),p=new Float64Array(N+1);
 let s=S0,y=Math.exp(2*S0)/6,v=Math.exp(2*S0)/3;
 const eq=(s,y,v)=>[v,Math.exp(2*s)*li(t+y)/n-v];
 for(let i=0;i<=N;i++) {
  p[i]=v;if(i===N)break;
  const k1=eq(s,y,v),k2=eq(s+H/2,y+H*k1[0]/2,v+H*k1[1]/2),
        k3=eq(s+H/2,y+H*k2[0]/2,v+H*k2[1]/2),k4=eq(s+H,y+H*k3[0],v+H*k3[1]);
  y+=H*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6;
  v+=H*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;s+=H;
 }
 return p;
}
const NUGRID=[
 ...Array.from({length:12},(_,i)=>.0001*(.05/.0001)**(i/11)),
 ...Array.from({length:90},(_,i)=>.05+(ZETA-.05)*(i+1)/90)
];
const PROFILES=NUGRID.map(profile);
function radial(p,x) {
 if(x<1e-4)return x*x/3;
 const a=(Math.log(x)-S0)/H;
 if(a>=N)return p[N];
 const i=Math.max(0,Math.floor(a)),f=a-i;
 return p[i]*(1-f)+p[i+1]*f;
}
function u2(nu,x) {
 if(nu<=NUGRID[0])return radial(PROFILES[0],x);
 let l=0,r=NUGRID.length-1;
 while(r-l>1){const m=(l+r)>>1;if(NUGRID[m]<nu)l=m;else r=m;}
 const w=(nu-NUGRID[l])/(NUGRID[r]-NUGRID[l]);
 return(1-w)*radial(PROFILES[l],x)+w*radial(PROFILES[r],x);
}
function haloV2(kind,r,R,A) {
 const x=r/R;
 if(kind==="Bose86") {
  const nu=K86/(R*R*A);
  return nu>ZETA?NaN:A*A*u2(nu,x);
 }
 if(kind==="Burkert") {
  const f=x<.001?4*x*x/3-x*x*x:(2*Math.log1p(x)+Math.log1p(x*x)-2*Math.atan(x))/x;
  return A*A*Math.max(0,f);
 }
 const f=x<.001?x*x/3:1-Math.atan(x)/x;
 return A*A*Math.max(0,f);
}
function chi2(g,kind,z,held) {
 const R=Math.exp(z[0]),A=Math.exp(z[1]);
 if(R<.08||R>120||A<3||A>650)return 1e14;
 let s=0;
 for(let i=0;i<g.length;i++){
  const x=g[i];
  const tot=x.vg*Math.abs(x.vg)+.5*x.vd*x.vd+.7*x.vb*x.vb+haloV2(kind,x.r,R,A);
  // Require positivity on both train and heldout radii without reading heldout velocities.
  if(!Number.isFinite(tot)||tot<0)return 1e12;
  if((i%4===0)!==held)continue;
  s+=((Math.sqrt(tot)-x.v)/x.dv)**2;
 }
 return s;
}
function nm(g,k,start) {
 let P=[start,[start[0]+.19,start[1]],[start[0],start[1]+.19]]
   .map(x=>({x:[...x],y:chi2(g,k,x,false)}));
 const ex=x=>({x,y:chi2(g,k,x,false)});
 for(let it=0;it<105;it++){
  P.sort((a,b)=>a.y-b.y);
  const [a,b,c]=P,center=[(a.x[0]+b.x[0])/2,(a.x[1]+b.x[1])/2];
  const refl=ex([2*center[0]-c.x[0],2*center[1]-c.x[1]]);
  if(refl.y<a.y){
   const expand=ex([3*center[0]-2*c.x[0],3*center[1]-2*c.x[1]]);
   P[2]=expand.y<refl.y?expand:refl;
  }else if(refl.y<b.y)P[2]=refl;
  else{
   const tar=refl.y<c.y?refl:c;
   const con=ex([center[0]+.5*(tar.x[0]-center[0]),center[1]+.5*(tar.x[1]-center[1])]);
   if(con.y<c.y)P[2]=con;
   else{
    P[1]=ex([a.x[0]+.5*(b.x[0]-a.x[0]),a.x[1]+.5*(b.x[1]-a.x[1])]);
    P[2]=ex([a.x[0]+.5*(c.x[0]-a.x[0]),a.x[1]+.5*(c.x[1]-a.x[1])]);
   }
  }
  if(Math.abs(P[0].y-P[1].y)<1e-9&&Math.abs(P[1].y-P[2].y)<1e-9
      &&Math.hypot(P[0].x[0]-P[2].x[0],P[0].x[1]-P[2].x[1])<1e-6)break;
 }
 P.sort((a,b)=>a.y-b.y);return P[0];
}
function fit(g,kind){
 const candidates=[];
 for(let a=0;a<=10;a++)for(let b=0;b<=10;b++){
  const z=[Math.log(.2)+a*Math.log(350)/10,Math.log(7)+b*Math.log(450/7)/10];
  candidates.push({z,c:chi2(g,kind,z,false)});
 }
 candidates.sort((a,b)=>a.c-b.c);
 const starts=[candidates[0].z];
 for(const p of candidates) {
  if(starts.length>=3)break;
  if(starts.every(s=>Math.hypot(s[0]-p.z[0],s[1]-p.z[1])>.3))starts.push(p.z);
 }
 const best=starts.map(s=>nm(g,kind,s)).sort((a,b)=>a.y-b.y)[0];
 const R=Math.exp(best.x[0]),A=Math.exp(best.x[1]);
 return{model:kind,R_kpc:R,amp_km_s:A,nu:kind==="Bose86"?K86/(R*R*A):null,
   chi2_train:best.y,chi2_holdout:chi2(g,kind,best.x,true)};
}
function main(){
 const groups=readData(process.argv[2]),sel=choose(groups),models=["Bose86","Burkert","PIS"];
 if(sel.length!==20)throw Error("Subset mismatch");
 if(Math.abs(Math.sqrt(u2(ZETA,10))-1.46422)>.01)throw Error("ODE regression check failed");
 const rows=sel.map(([id,g])=>({id,n:g.length,n_train:g.filter((_,i)=>i%4!==0).length,
   n_holdout:g.filter((_,i)=>i%4===0).length,fits:models.map(k=>fit(g,k))}));
 const sums={};
 for(let j=0;j<3;j++){
  const name=models[j];
  sums[name]={train_chi2:rows.reduce((s,x)=>s+x.fits[j].chi2_train,0),
   holdout_chi2:rows.reduce((s,x)=>s+x.fits[j].chi2_holdout,0),
   best_holdout:rows.filter(x=>x.fits[j].chi2_holdout
     ===Math.min(...x.fits.map(v=>v.chi2_holdout))).length};
 }
 const result={data:"SPARC Lelli2016c mass-model rows",mass_eV_fixed:86,
  M_L_disk:.5,M_L_bulge:.7,nu_templates:NUGRID.length,
  selection:"first 20 alphabetic, >=15 points, median nominal DM2 fraction>.7, no negative nominal DM2",
  train_points:rows.reduce((s,x)=>s+x.n_train,0),
  holdout_points:rows.reduce((s,x)=>s+x.n_holdout,0),
  comparisons:sums,per_galaxy:rows,
  disclaimer:"Preliminary, selection-biased; fixed M/L, random errors only, no distance/inclination nuisance, no rigorous global inference."};
 console.log(JSON.stringify(result,null,2));
}
main();
