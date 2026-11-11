"""Exact finite-temperature edge selection for block-aligned straight-edged Shakti samples (order n Aztec II + 2(n+1)
edge plaquettes). Free edge plaquettes are monomers on the outer face with fugacity z=exp(J/2T); the partition function
is the Pfaffian Pf(K + z^2 C), C_ab = +1 for edge cells a<b in counterclockwise order starting from the west-side group
(validated against brute-force subset sums for n=5,7: pf_test7.py). <eps(p)> relative to the bottom vertex of the Aztec
diamond follows from the seam-flipped Pfaffian. Usage: python edge_exact.py n x1 x2 ..."""
import sys, math, cmath, pickle, time, numpy as np
from pfapack.pfaffian import pfaffian
from pf_edge import geometry, skew
from exact_v import Region, start_vertex, black
from vform2 import r4
n=int(sys.argv[1]); XS=[float(a) for a in sys.argv[2:]]
A,box,E=geometry(n); cells,idx,K=skew(box)
cx=np.mean([c[0] for c in box]); cy=np.mean([c[1] for c in box])
Eo=sorted(E,key=lambda c: math.atan2(c[1]+.5-cy,c[0]+.5-cx)); Eo=[Eo[-1]]+Eo[:-1]; ie=[idx[e] for e in Eo]
C=np.zeros_like(K)
for a in range(len(ie)):
    for b in range(a+1,len(ie)): C[ie[a],ie[b]]=1; C[ie[b],ie[a]]=-1
RA=Region(A); v0,hs=start_vertex(A)
acx=np.mean([c[0]+.5 for c in A]); acy=np.mean([c[1]+.5 for c in A])
cen=[c for c in A if math.hypot(c[0]+.5-acx,c[1]+.5-acy)<=4]
verts=sorted({(c[0]+a,c[1]+b) for c in cen for a in (0,1) for b in (0,1)})
seams={}
for p in verts:
    pts=RA.path(v0,p); h0=hs; fl=[]
    for a,b in zip(pts,pts[1:]):
        L,R=RA.edge_cells(a,b); h0+= 1 if black(L) else -1
        if L in box and R in box: fl.append((idx[L],idx[R]))
    seams[p]=(h0,fl)
out={}
for x in XS:
    t=time.time(); z=math.exp(x/2); M=K+z*z*C
    s,l=np.linalg.slogdet(M); lnZ=0.5*l
    dz=1e-4; s2,l2=np.linalg.slogdet(K+(z*math.exp(dz))**2*C); nfree=(0.5*l2-lnZ)/dz   # = z dlnZ/dz
    Minv=np.linalg.inv(M); eps={}
    Np=0.0
    for e in E:
        i=idx[e]
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            q=(e[0]+d[0],e[1]+d[1])
            if q in idx: j=idx[q]; Np+=M[i,j]*Minv[j,i]
    nfree=len(E)-Np
    for p,(h0,fl) in seams.items():
        # Pf(M + U C U^T)/Pf(M) = Pf(C^{-1} + U^T M^{-1} U) / Pf(C^{-1})  (low-rank seam update)
        k=len(fl); e0=1 if ((h0-r4(p))//4)%2==0 else -1
        if k==0: eps[p]=e0; continue
        cols=[]; Ci=np.zeros((2*k,2*k))
        for m,(i,j) in enumerate(fl):
            cols+= [i,j]; v=M[i,j]; Ci[2*m,2*m+1]=1/(2*v); Ci[2*m+1,2*m]=-1/(2*v)   # C block [[0,-2v],[2v,0]] -> inverse
        S=Ci+Minv[np.ix_(cols,cols)]; S=(S-S.T)/2
        eps[p]=e0*pfaffian(S)/pfaffian(Ci)
    psi={c:0.5*sum(cmath.exp(1j*math.pi*r4((c[0]+a,c[1]+b))/4)*eps[(c[0]+a,c[1]+b)] for a in (0,1) for b in (0,1)) for c in cen}
    amp=np.mean([abs(v) for v in psi.values()])/1.3066
    out[x]=dict(lnZ=lnZ,paired=1-nfree/len(E),npaired=len(E)-nfree,amp=amp,psi=psi)
    print(f'n={n} x=J/T={x:6.2f}  paired fraction {1-nfree/len(E):.4f}  (paired plaquettes {len(E)-nfree:.2f} of {len(E)})  amp {amp:.4f}  {time.time()-t:.0f}s',flush=True)
pickle.dump(dict(n=n,E=len(E),out=out),open(f'exact/edge_exact_n{n}.pkl','wb'))
