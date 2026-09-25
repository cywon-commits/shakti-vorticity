import io, contextlib, numpy as np
from exact_defect import *
_p=open('placement2.py').read(); exec(_p[_p.index('def enum_cells'):_p.index('for n in range(1,6)')])
for n,pl in [(5,'I'),(5,'II'),(6,'I')]:
    cells,sh=aztec_cells(n,pl); c1=(sh[0],0) if black((sh[0],0)) else (sh[0]-1,0)
    for dxy in [(1,0),(3,0),(2,1),(3,2),(1,2),(4,1)]:
        c2=(c1[0]+dxy[0],c1[1]+dxy[1])
        O1,O2=O_at(c1),O_at(c2)
        if set(O1)&set(O2) or not set(O2)<=cells: continue
        rest=cells-set(O1)-set(O2)
        N=sum(1 for _ in enum_cells(rest)) if n<6 else None
        B=sorted(c for c in rest if black(c)); W=sorted(c for c in rest if not black(c)); wi={c:k for k,c in enumerate(W)}
        K=np.zeros((len(B),len(W)),complex)
        for i,b in enumerate(B):
            for di,dj in ((1,0),(-1,0),(0,1),(0,-1)):
                w=(b[0]+di,b[1]+dj)
                if w in wi: K[i,wi[w]]=1 if dj==0 else 1j
        d=np.linalg.det(K)
        print(n,pl,dxy,'brute',N,'|det K|',round(abs(d),6))
