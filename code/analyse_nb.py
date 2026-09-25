import pickle, glob, math, cmath, numpy as np
from vform2 import r4
def fourier(V,n,sh,R):
    O=n+1; acc={q:0j for q in [(1,1),(1,-1),(-1,1),(-1,-1)]}; cnt=0
    for i in range(-n,n):
        for j in range(-n,n):
            if abs(2*i+1)+abs(2*j+1)<=2*n and math.hypot(i+.5,j+.5)<=R:
                cnt+=1; x,y=i+sh[0],j+sh[1]
                for q in acc: acc[q]+=V[i+O,j+O]*cmath.exp(-1j*math.pi/2*(q[0]*x+q[1]*y))
    return {q:a/cnt for q,a in acc.items()}
def amp(F): return math.sqrt(sum(abs(a)**2 for a in F.values()))
rows=[]
for f in sorted(glob.glob('nb/a_*.pkl')):
    D=pickle.load(open(f,'rb')); n=D['n']; sh=D['shift']; R=0.25*n/math.sqrt(2)
    F=fourier(D['vmean'],n,sh,R); A=amp(F)
    Ab=[amp(fourier(Vb,n,sh,R)) for Vb in D['vblocks']]; err=np.std(Ab,ddof=1)/math.sqrt(len(Ab))
    # H/V split
    O=n+1; s={0:[],1:[]}
    for i in range(-n,n):
        for j in range(-n,n):
            if abs(2*i+1)+abs(2*j+1)<=2*n and math.hypot(i+.5,j+.5)<=R: s[(i+sh[0]+j+sh[1])%2].append(D['vmean'][i+O,j+O])
    AH=math.sqrt(max(np.mean(np.square(s[0]))-1/D['N'],0)); AV=math.sqrt(max(np.mean(np.square(s[1]))-1/D['N'],0))
    # mean central height relative to its residue: (h - r)/4 has mean m; phase variable
    r=r4((sh[0],sh[1])); hc=D['hcenter']; dev=(hc-r)/4.0
    rows.append(dict(n=n,P=D['placement'],A=A,err=err,alpha=F[(1,1)],beta=F[(1,-1)],AH=AH,AV=AV,asym=(AH**2-AV**2)/(AH**2+AV**2) if AH+AV>0 else 0,
                     varh=hc.var(),hbar=hc.mean(),frac=np.mean(dev-np.floor(dev)),ec=D['emean'][O,O],V=sorted(set(D['tot'].tolist())),N=D['N']))
rows.sort(key=lambda r:(r['P'],r['n']))
pickle.dump(rows,open('nb/rows.pkl','wb'))
for r in rows:
    print('%-2s n=%3d A=%.3f±%.3f asym=%+.2f  arg(alpha)=%+6.1f° arg(beta)=%+6.1f°  <h_c>=%.2f  Var h=%.2f  V=%s'%(r['P'],r['n'],r['A'],r['err'],r['asym'],math.degrees(cmath.phase(r['alpha'])),math.degrees(cmath.phase(r['beta'])),r['hbar'],r['varh'],r['V']))
for P in ('I',):
    R=[r for r in rows if r['P']==P and r['n']>=16]
    x=np.log([r['n'] for r in R]); y=np.log([r['A'] for r in R]); w=1/np.array([r['err']/r['A'] for r in R])**2
    p,cov=np.polyfit(x,y,1,w=np.sqrt(w),cov=True); print(P,'weighted fit exponent %.3f ± %.3f'%(-p[0],math.sqrt(cov[0,0])))
    for k in range(4):
        Rk=[r for r in R if r['n']%4==k]; xk=np.log([r['n'] for r in Rk]); yk=np.log([r['A'] for r in Rk])
        pk,ck=np.polyfit(xk,yk,1,cov=True); print('   n≡%d exponent %.3f ± %.3f'%(k,-pk[0],math.sqrt(ck[0,0])))
    xv=np.log([r['n'] for r in R]); 
    for k in (0,3,1,2):
        Rk=[r for r in R if r['n']%4==k]; pv=np.polyfit(np.log([r['n'] for r in Rk]),[r['varh'] for r in Rk],1); print('   n≡%d dVar/dln n = %.3f (theory 8/pi^2 = %.3f)'%(k,pv[0],8/math.pi**2))
