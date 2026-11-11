import math, itertools, numpy as np, sys
from pf_edge import *
from exact_v import black, Region
n=int(sys.argv[1]); A,box,E=geometry(n); cells,idx,K=skew(box)
cx=np.mean([c[0] for c in box]); cy=np.mean([c[1] for c in box])
Eo=sorted(E,key=lambda c: math.atan2(c[1]+.5-cy,c[0]+.5-cx)); Eo=[Eo[-1]]+Eo[:-1]   # start at the west-side group
ie=[idx[e] for e in Eo]
def Zpf(z):
    Mx=K.copy()
    for a in range(len(ie)):
        for b in range(a+1,len(ie)): Mx[ie[a],ie[b]]+=z*z; Mx[ie[b],ie[a]]-=z*z
    s,l=np.linalg.slogdet(Mx); return 0.5*l
Eb=[e for e in E if black(e)]; Ew=[e for e in E if not black(e)]; terms=[]
for kb in range(len(Eb)+1):
    for Sb in itertools.combinations(Eb,kb):
        for Sw in itertools.combinations(Ew,kb):
            try: terms.append((len(E)-2*kb,Region(A|set(Sb)|set(Sw)).ld0))
            except Exception: pass
for z in (0.3,1.0,2.5):
    ref=np.log(sum(math.exp(l)*z**f for f,l in terms))
    print(n,z,'brute lnZ',ref,'pfaffian',Zpf(z))
