"""Local proof checks for v(d) = s(d) * (-1)^floor(h(p_d)/4)."""
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from vort_rule import *
from vform2 import S_TAB, r4
def corners(d): return {(x+dx,y+dy) for (x,y) in d for dx in (0,1) for dy in (0,1)}
def perim_edges(d):
    E=Counter()
    for (x,y) in d:
        for e in [((x,y),(x+1,y)),((x,y+1),(x+1,y+1)),((x,y),(x,y+1)),((x+1,y),(x+1,y+1))]: E[key(e)]+=1
    return {e for e,c in E.items() if c==1}
def black_left(p,q):
    d=(q[0]-p[0],q[1]-p[1])
    left={(1,0):(p[0],p[1]),(-1,0):(q[0],q[1]-1),(0,1):(p[0]-1,p[1]),(0,-1):(q[0],q[1])}[d]
    return 1 if black(left) else -1
def heights_on(ds, p0):
    """heights on the perimeters of the dominoes ds (connected through shared vertices), relative h(p0)=r4(p0)"""
    E=set().union(*[perim_edges(d) for d in ds]); adj=defaultdict(list)
    for a,b in E: adj[a].append(b); adj[b].append(a)
    h={p0:r4(p0)}; Q=[p0]
    for p in Q:
        for q in adj[p]:
            val=h[p]+black_left(p,q)
            if q not in h: h[q]=val; Q.append(q)
            elif h[q]!=val: raise Exception('inconsistent')
    return h
def vform(d,h):
    d=tuple(sorted(d)); p=d[0]; typ='H' if d[1][0]!=d[0][0] else 'V'
    assert h[p]%4==r4(p); return S_TAB[(typ,p[0]%4,p[1]%4)]*(-1)**(h[p]//4)
res=Counter(); nocontact=0
base=[((i,j),(i+dx,j+dy)) for i in range(0,4) for j in range(0,4) for dx,dy in [(1,0),(0,1)]]   # all positions mod 4
for d1 in base:
    A=info2(d1)
    for i in range(-5,9):
        for j in range(-5,9):
            for dx,dy in [(1,0),(0,1)]:
                d2=((i,j),(i+dx,j+dy))
                if set(d2)&set(d1): continue
                B=info2(d2); js,common,kinds=coupling(A,B)
                if not common: continue
                J=js.pop()
                if not (corners(d1)&corners(d2)): nocontact+=1; continue
                h=heights_on([d1,d2],tuple(sorted(d1))[0])
                ok = vform(d1,h)*vform(d2,h)==J
                res[ok]+=1
print('coupled pairs checked:',dict(res),' coupled pairs without a common chess vertex:',nocontact)
