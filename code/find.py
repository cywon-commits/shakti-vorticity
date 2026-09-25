import itertools
# Shakti (basketweave): block (a,b) = [2a,2a+2]x[2b,2b+2]; H if (a+b)%2==0
def rects(a,b):
    x,y=2*a,2*b
    if (a+b)%2==0: return [frozenset({(x,y),(x+1,y)}), frozenset({(x,y+1),(x+1,y+1)})]
    return [frozenset({(x,y),(x,y+1)}), frozenset({(x+1,y),(x+1,y+1)})]
def isL(r1,r2):
    u=r1|r2
    if len(u)!=4: return False
    # connected & not straight & not square
    xs={p[0] for p in u}; ys={p[1] for p in u}
    if sorted((len(xs),len(ys)))!=[2,3]: return False
    adj=sum(1 for p in u for q in u if abs(p[0]-q[0])+abs(p[1]-q[1])==1)//2
    return adj==3 and (r1 & r2)==frozenset()
R=10
rs={}
for a in range(-R,R):
  for b in range(-R,R):
    for k,r in enumerate(rects(a,b)): rs[(a,b,k)]=r
# candidate: block (a,b) -> base cell (p,q) = (a*m1+b*m2, a*m3+b*m4); k -> offset 0 or diagonal d
best=None
for m in itertools.product(range(-2,3),repeat=4):
  for d in [(1,1),(1,-1)]:
    for flipH in [0,1]:
     for flipV in [0,1]:
      cell={}
      for (a,b,k),r in rs.items():
        p=m[0]*a+m[1]*b; q=m[2]*a+m[3]*b
        f = flipH if (a+b)%2==0 else flipV
        kk=k^f
        cell[(a,b,k)]=(p+kk*d[0], q+kk*d[1])
      if len(set(cell.values()))!=len(cell): continue
      inv={v:kk for kk,v in cell.items()}
      ok=True
      for key,c in cell.items():
        a,b,k=key
        if abs(a)>R-3 or abs(b)>R-3: continue
        for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)]:
          n=(c[0]+dx,c[1]+dy)
          if n not in inv or not isL(rs[key],rs[inv[n]]): ok=False;break
        if not ok: break
      if ok: print(m,d,flipH,flipV)
