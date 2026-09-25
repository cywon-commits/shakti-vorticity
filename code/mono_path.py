"""Forced-elbow ILP: monomer pair created by sliding along an alternating path from a ground-state tiling (trivial sector)."""
import sys, math, pickle
from mono_ilp import *
P=int(sys.argv[1]); D=float(sys.argv[2]) if len(sys.argv)>2 else 3.0
COST=dict(COST0); [COST.__setitem__(k,D) for k in ('z4II','z4q2','z3q3','z2q')]
M=Model(P); T=M.T
cen=lambda i:(sum(p[0] for p in T.plaqs[i])/2+.5, sum(p[1] for p in T.plaqs[i])/2+.5)
def td(a,b):
    dx=abs(a[0]-b[0]); dy=abs(a[1]-b[1]); return math.hypot(min(dx,P-dx),min(dy,P-dy))
doms=M.tiling(); E0,sp,labs,el=M.solve(doms,cost=COST); print('reference E',E0,'tiles',len(doms),flush=True)
part={}
for c in doms: part[c[0]]=c; part[c[1]]=c
nbr=defaultdict(set)
for c in T.cands: nbr[c[0]].add((c[1],c)); nbr[c[1]].add((c[0],c))
c0=doms[0]; A,m=c0[0],c0[1]; cur=[c for c in doms if c is not c0]
res=[]
for step in range(int(sys.argv[3]) if len(sys.argv)>3 else 12):
    E,sp,labs,el=M.solve(cur,cost=COST)
    extra=Counter(l for p,l in labs.items() if l not in ('z2ok','z3happy','z4I') and not (l=='z3unhappy' and p in el))
    r=td(cen(A),cen(m)); res.append((r,E-len(cur),extra))
    print(f'step {step}: r={r:5.2f} excess E={E-len(cur):4.1f} extra {dict(extra)}',flush=True)
    pickle.dump(res,open(f'mono/path{P}_D{D:g}.pkl','wb'))
    # move monomer m two plaquettes away from A
    pc={}; 
    for c in cur: pc[c[0]]=c; pc[c[1]]=c
    best=None
    for X,cx in nbr[m]:
        if X==A or X not in pc: continue
        dXY=pc[X]; Y=dXY[1] if dXY[0]==X else dXY[0]
        if Y==A: continue
        dd=td(cen(A),cen(Y))
        if best is None or dd>best[0]: best=(dd,X,Y,dXY,cx)
    dd,X,Y,dXY,cx=best
    cur=[c for c in cur if c is not dXY]+[cx]; m=Y
