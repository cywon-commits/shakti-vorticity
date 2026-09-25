"""Chessboard <-> Shakti map (placement I), domino shuffling, short-arm graph with frozen boundary ring."""
import math, random
from collections import Counter
_g={}
exec(open('embed.py').read().split("print('inconsistencies'")[0],_g)
_inv0={c:k for k,c in _g['coord'].items()}
INV={(-c[0]-1,-c[1]-1):k for c,k in _inv0.items()}          # same convention as the figures
RECTS=_g['rects']
# periodicity of the map: cell + (2,2) and cell + (2,-2) translate the block index
def _T(v):
    ts={ (INV[(i+v[0],j+v[1])][0]-INV[(i,j)][0], INV[(i+v[0],j+v[1])][1]-INV[(i,j)][1]) for i in range(-3,3) for j in range(-3,3)}
    assert len(ts)==1; return ts.pop()
T1,T2=_T((2,2)),_T((2,-2))
def inv(c):
    i,j=c
    # c = r + p(2,2) + q(2,-2) with r in a small window
    p=math.floor((i+j)/4); q=math.floor((i-j)/4)
    r=(i-2*p-2*q, j-2*p+2*q)
    a,b,k=INV[r]
    return (a+p*T1[0]+q*T2[0], b+p*T1[1]+q*T2[1], k)
def rect(c):
    a,b,k=inv(c); return RECTS(a,b)[k]
def black(c): return (c[0]+c[1])%2==0
def partner(c):
    a,b,_=inv(c)
    return next(c2 for c2 in [(c[0]+1,c[1]+1),(c[0]-1,c[1]-1),(c[0]+1,c[1]-1),(c[0]-1,c[1]+1)] if inv(c2)[:2]==(a,b))
def short_cell(x,y):
    """short arm of the L for domino {x,y}: y if y is orthogonally adjacent to partner(x), else x"""
    p=partner(x); return y if abs(y[0]-p[0])+abs(y[1]-p[1])==1 else x
def pcen(c):
    r=rect(c); return (sum(p[0] for p in r)/2+.5, sum(p[1] for p in r)/2+.5)
def aztec(n, shift=None):
    if shift is None: shift=((n+1)%2,0)       # placement I
    return {(i+shift[0],j+shift[1]) for i in range(-n,n) for j in range(-n,n) if abs(i+.5)+abs(j+.5)<=n}, shift
# ---------- uniform random tiling by domino shuffling (unshifted coordinates) ----------
def shuffle(n, rng=random):
    dom={}   # key: (type, anchor cell) ; horizontal anchor=left cell, vertical anchor=bottom cell
    for k in range(n):
        # 1. destruction of colliding pairs
        S={(t,a) for t,a in dom}
        kill=set()
        for t,a in S:
            if t=='S' and ('N',(a[0],a[1]-1)) in S: kill|={('S',a),('N',(a[0],a[1]-1))}
            if t=='E' and ('W',(a[0]+1,a[1])) in S: kill|={('E',a),('W',(a[0]+1,a[1]))}
        for x in kill: dom.pop(x)
        # 2. sliding
        mv={'N':(0,1),'S':(0,-1),'E':(1,0),'W':(-1,0)}
        dom={(t,(a[0]+mv[t][0],a[1]+mv[t][1])):1 for t,a in dom}
        # 3. creation in empty 2x2 blocks of the order-(k+1) diamond
        m=k+1; occ=set()
        for t,a in dom:
            occ.add(a); occ.add((a[0]+1,a[1]) if t in 'NS' else (a[0],a[1]+1))
        cells=[(i,j) for j in range(m-1,-m-1,-1) for i in range(-m,m) if abs(i+.5)+abs(j+.5)<=m]
        for (i,j) in cells:
            if (i,j) in occ: continue
            blk=[(i,j),(i+1,j),(i,j-1),(i+1,j-1)]
            assert all(b not in occ for b in blk),'block not empty'
            if rng.random()<.5: dom[('N',(i,j))]=1; dom[('S',(i,j-1))]=1
            else: dom[('W',(i,j-1))]=1; dom[('E',(i+1,j-1))]=1
            occ|=set(blk)
    out=[]
    for t,a in dom:
        out.append((a,(a[0]+1,a[1])) if t in 'NS' else (a,(a[0],a[1]+1)))
    return out
# ---------- graph ----------
class DSU:
    def __init__(s): s.p={}
    def f(s,x):
        s.p.setdefault(x,x)
        while s.p[x]!=x: s.p[x]=s.p[s.p[x]]; x=s.p[x]
        return x
    def u(s,a,b): s.p[s.f(a)]=s.f(b)
def analyse(doms, cells):
    """returns dict cell->spin (+1/-1) for ring-connected short arms, set of free short arms, conflicts"""
    short={short_cell(*d) for d in doms}
    def touches(c):
        return any(nb not in cells for nb in [(c[0]+1,c[1]),(c[0]-1,c[1]),(c[0],c[1]+1),(c[0],c[1]-1),partner(c)])
    allsq=[p for c in cells for p in rect(c)]
    mx=sum(p[0]+.5 for p in allsq)/len(allsq); my=sum(p[1]+.5 for p in allsq)/len(allsq)
    ring=sorted([c for c in short if touches(c)], key=lambda c: math.atan2(pcen(c)[1]-my,pcen(c)[0]-mx))
    adj={c:[] for c in short}
    for c in short:
        for d in [(1,0),(-1,0),(0,1),(0,-1)]:
            nb=(c[0]+d[0],c[1]+d[1])
            if nb in short: adj[c].append(nb)
    # ring spins alternate along the ring; reference: first ring node after angle -pi (fixed, since ring is frozen)
    spin={}
    for idx,c in enumerate(ring): spin[c]=1 if idx%2==0 else -1
    for a,b in zip(ring,ring[1:]+ring[:1]):
        adj[a].append(b); adj[b].append(a)
    Q=list(ring); conflicts=0
    for u in Q:
        for v in adj[u]:
            if v not in spin: spin[v]=-spin[u]; Q.append(v)
            elif spin[v]==spin[u]: conflicts+=1
    free=short-set(spin)
    return spin, free, ring, conflicts//2
