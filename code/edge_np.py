"""<N_p> (paired edge plaquettes) from Pfaffian edge probabilities P(i~j) = M_ij (M^{-1})_ji (no cancellation)."""
import sys, math, numpy as np
from pf_edge import geometry, skew
n=int(sys.argv[1]); XS=[float(a) for a in sys.argv[2:]]
A,box,E=geometry(n); cells,idx,K=skew(box)
cx=np.mean([c[0] for c in box]); cy=np.mean([c[1] for c in box])
Eo=sorted(E,key=lambda c: math.atan2(c[1]+.5-cy,c[0]+.5-cx)); Eo=[Eo[-1]]+Eo[:-1]; ie=[idx[e] for e in Eo]
C=np.zeros_like(K)
for a in range(len(ie)):
    for b in range(a+1,len(ie)): C[ie[a],ie[b]]=1; C[ie[b],ie[a]]=-1
for x in XS:
    M=K+math.exp(x)*C; Mi=np.linalg.inv(M); Np=0.0
    for e in E:
        i=idx[e]
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            q=(e[0]+d[0],e[1]+d[1])
            if q in idx: j=idx[q]; Np+=M[i,j]*Mi[j,i]
    print(n,x,Np,math.log(Np)+x-(n+2)*math.log(2) if Np>0 else None,flush=True)
