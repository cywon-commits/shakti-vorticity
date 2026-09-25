import sys, random, math
from mono_ilp import *
P=int(sys.argv[1]); M=Model(P); T=M.T
cen=lambda i:(sum(p[0] for p in T.plaqs[i])/2+.5, sum(p[1] for p in T.plaqs[i])/2+.5)
def tdist(a,b):
    dx=abs(a[0]-b[0]); dy=abs(a[1]-b[1]); dx=min(dx,P-dx); dy=min(dy,P-dy); return math.hypot(dx,dy)
A=next(i for i in range(len(T.plaqs)) if T.horiz[i])
rng=random.Random(1)
Bs=sorted([i for i in range(len(T.plaqs)) if not T.horiz[i]],key=lambda b:tdist(cen(A),cen(b)))
seen=set()
for B in Bs:
    r=round(tdist(cen(A),cen(B)),2)
    if r in seen: continue
    seen.add(r)
    for trial in range(2):
        d=M.tiling(removed=(A,B),rng=random.Random(trial))
        E,sp,labs,el=M.solve(d)
        extra=[(p,l) for p,l in labs.items() if l not in ('z2ok','z3happy','z4I') and not (l=='z3unhappy' and p in el)]
        print(f'r={r:5.2f} trial {trial}: excess E = {E-len(d):.0f}; defects:',Counter(l for p,l in extra),flush=True)
