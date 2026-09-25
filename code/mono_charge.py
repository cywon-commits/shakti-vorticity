import sys, math
from mono_ilp import *
P=20; M=Model(P); T=M.T
cen=lambda i:(sum(p[0] for p in T.plaqs[i])/2+.5, sum(p[1] for p in T.plaqs[i])/2+.5)
def dv(a,b):
    dx=(b[0]-a[0]+P/2)%P-P/2; dy=(b[1]-a[1]+P/2)%P-P/2; return math.hypot(dx,dy)
doms=M.tiling(); E0,sp0,_,_=M.solve(doms)
q0={p:T.vstate(p,sp0)[1] for p in T.verts}
nbr=defaultdict(set)
for c in T.cands: nbr[c[0]].add((c[1],c)); nbr[c[1]].add((c[0],c))
c0=doms[0]; A,m=c0[0],c0[1]; cur=[c for c in doms if c is not c0]
for step in range(6):
    pc={}
    for c in cur: pc[c[0]]=c; pc[c[1]]=c
    best=None
    for X,cx in nbr[m]:
        if X==A or X not in pc: continue
        d=pc[X]; Y=d[1] if d[0]==X else d[0]
        if Y==A: continue
        dd=dv(cen(A),cen(Y))
        if best is None or dd>best[0]: best=(dd,X,Y,d,cx)
    dd,X,Y,d,cx=best; cur=[c for c in cur if c is not d]+[cx]; m=Y
    if step in (1,3,5):
        E,sp,labs,el=M.solve(cur)
        q={p:T.vstate(p,sp)[1] for p in T.verts}
        dq={p:q[p]-q0[p] for p in T.verts if q[p]!=q0[p]}
        cA,cm=cen(A),cen(m)
        def disk(c,R): return sum(v for p,v in dq.items() if dv(c,p)<=R)
        defects=[(p,labs[p],q[p]) for p in T.verts if labs[p] not in ('z2ok','z3happy','z4I') and not (labs[p]=='z3unhappy' and p in el)]
        print(f'r={dv(cA,cm):.2f} E-N={E-len(cur)}  net Δq total={sum(dq.values())}  Δq within R=2,3 of monomer A: {disk(cA,2)},{disk(cA,3)}  of monomer B: {disk(cm,2)},{disk(cm,3)}')
        print('   defects:',[(p,l,qq, round(dv(cA,p),1), round(dv(cm,p),1)) for p,l,qq in defects])
        vor=sum(v for p,v in dq.items() if dv(cA,p)<dv(cm,p)); print('   Voronoi charge: A side',vor,' B side',sum(dq.values())-vor)
