import math
exec(open('find.py').read().split('# candidate')[0])
keys=list(rs)
cen={k:(sum(p[0] for p in r)/2+.5, sum(p[1] for p in r)/2+.5) for k,r in rs.items()}
idx={}
for k,r in rs.items():
    for p in r: idx[p]=k
def nbrs(k):
    out=set()
    for p in rs[k]:
        for dx,dy in [(1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1),(2,0),(0,2),(-2,0),(0,-2)]:
            q=(p[0]+dx,p[1]+dy)
            if q in idx and idx[q]!=k and isL(rs[k],rs[idx[q]]): out.add(idx[q])
    return out
def ordered(k):
    c=cen[k]; ns=list(nbrs(k))
    return sorted(ns,key=lambda n: math.atan2(cen[n][1]-c[1],cen[n][0]-c[0]))
root=(0,0,0)
coord={root:(0,0)}
dirs=[(1,0),(0,1),(-1,0),(0,-1)]
from collections import deque
Q=deque([root]); bad=0
while Q:
    k=Q.popleft(); a,b,_=k
    if max(abs(a),abs(b))>R-2: continue
    ns=ordered(k)
    assert len(ns)==4,(k,ns)
    # find a neighbor with known coord to anchor
    if k==root: start=0; d0=0
    else:
        start=next(i for i,n in enumerate(ns) if n in coord)
        dv=(coord[ns[start]][0]-coord[k][0],coord[ns[start]][1]-coord[k][1]); d0=dirs.index(dv)
    for i in range(4):
        n=ns[(start+i)%4]; d=dirs[(d0+i)%4]
        c=(coord[k][0]+d[0],coord[k][1]+d[1])
        if n in coord:
            if coord[n]!=c: bad+=1
        else: coord[n]=c; Q.append(n)
print('inconsistencies',bad)
for a in range(-2,3):
  for b in range(-2,3):
    print((a,b),'H' if (a+b)%2==0 else 'V',[coord.get((a,b,k)) for k in (0,1)])
