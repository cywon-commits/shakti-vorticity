"""Consistency of the O-tetromino (diagonal defect) template with the L-tetromino templates."""
import sys, io, contextlib
from collections import Counter, defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    from vort import tet_of, template_for, key, sgn, is_island, RECTS
from core import rect, partner, black, aztec, inv
# O templates in block coordinates [0,2]^2, state o=+1 (centre two-in, charge +2), read off the figure
O_H=[((1,2),(0,2)),((1,2),(2,2)),((0,2),(0,1)),((0,0),(0,1)),((0,1),(1,1)),((2,1),(1,1)),((2,2),(2,1)),((2,0),(2,1)),((1,0),(0,0)),((1,0),(2,0))]
O_V=[((0,2),(1,2)),((2,2),(1,2)),((0,1),(0,2)),((0,1),(0,0)),((1,2),(1,1)),((1,0),(1,1)),((2,1),(2,2)),((2,1),(2,0)),((0,0),(1,0)),((2,0),(1,0))]
# external islands (read off the second figure), state o=+1
O_H_EXT=[((0,2),(0,3)),((-1,2),(0,2)),((1,3),(1,2)),((2,2),(2,3)),((3,2),(2,2)),((-1,0),(0,0)),((0,0),(0,-1)),((1,-1),(1,0)),((3,0),(2,0)),((2,0),(2,-1))]
O_V_EXT=[((0,3),(0,2)),((0,2),(-1,2)),((2,3),(2,2)),((2,2),(3,2)),((-1,1),(0,1)),((3,1),(2,1)),((0,0),(-1,0)),((0,-1),(0,0)),((2,0),(3,0)),((2,-1),(2,0))]
USE_EXT=True
def o_template(c1,c2):
    sq=set(rect(c1))|set(rect(c2)); x0=min(p[0] for p in sq); y0=min(p[1] for p in sq)
    assert len(sq)==4 and max(p[0] for p in sq)==x0+1 and max(p[1] for p in sq)==y0+1, 'not a 2x2 block'
    horiz = len({p[1] for p in rect(c1)})==1
    T=O_H if horiz else O_V
    out={}
    for a,b in T:
        e=((a[0]+x0,a[1]+y0),(b[0]+x0,b[1]+y0)); assert is_island(e),('missing island',e)
        out[key(e)]=(sgn(e),'per')
    if USE_EXT:
        for a,b in (O_H_EXT if horiz else O_V_EXT):
            e=((a[0]+x0,a[1]+y0),(b[0]+x0,b[1]+y0)); assert is_island(e),('missing external island',e)
            out[key(e)]=(sgn(e),'ext')
    return out
_cache={}
def tile_template(t):
    if t in _cache: return _cache[t]
    kind,c1,c2=t
    if kind=='D':
        rs_,rl,s=tet_of((c1,c2)); tp,fr=template_for(rs_,rl)
    else:
        tp=o_template(c1,c2)
    _cache[t]=tp; return tp
def analyse(tiles):
    owner=defaultdict(list)
    for k,t in enumerate(tiles):
        for e,(b,kind) in tile_template(t).items(): owner[e].append((k,b,kind))
    adj=defaultdict(dict); conflicts=[]
    for e,lst in owner.items():
        for x in range(len(lst)):
            for y in range(x+1,len(lst)):
                (i,bi,ki),(j,bj,kj)=lst[x],lst[y]
                if i==j: continue
                J=bi*bj
                if j in adj[i] and adj[i][j]!=J: conflicts.append((tiles[i],tiles[j]))
                adj[i][j]=J; adj[j][i]=J
    v={}; fr=0; comps=0
    for s in range(len(tiles)):
        if s in v: continue
        comps+=1; v[s]=1; Q=[s]
        for u in Q:
            for w,J in adj[u].items():
                if w not in v: v[w]=v[u]*J; Q.append(w)
                elif v[w]!=v[u]*J: fr+=1
    return conflicts, fr//2, comps, v
def enum_with_defects(cells, maxdef=None):
    cl=sorted(cells,key=lambda c:(c[1],c[0]))
    def rec(used,tiles,nd):
        free=next((c for c in cl if c not in used),None)
        if free is None: yield list(tiles); return
        for dd in [(1,0),(0,1)]:
            o=(free[0]+dd[0],free[1]+dd[1])
            if o in cells and o not in used: yield from rec(used|{free,o},tiles+[('D',)+tuple(sorted((free,o)))],nd)
        p=partner(free)
        if p in cells and p not in used and (maxdef is None or nd<maxdef):
            yield from rec(used|{free,p},tiles+[('O',)+tuple(sorted((free,p)))],nd+1)
    return rec(frozenset(),[],0)
if __name__=='__main__':
    for name,cells in [('Aztec n=2 (I)',aztec(2)[0]),('Aztec n=3 (I)',aztec(3)[0]),('Aztec n=3 (II)',aztec(3,(1,0))[0]),('Aztec n=4 (I)',aztec(4)[0]),
                       ('4x4 board',{(i,j) for i in range(4) for j in range(4)}),('4x6 board',{(i,j) for i in range(4) for j in range(6)}),('6x6 board',{(i,j) for i in range(6) for j in range(6)})]:
        C=Counter(); bad=Counter(); examples=[]
        for tiles in enum_with_defects(cells):
            nO=sum(t[0]=='O' for t in tiles)
            conf,fr,comps,v=analyse(tiles)
            status='conflict' if conf else ('frustrated' if fr else 'ok')
            C[(nO//2,status)]+=1
            if status!='ok' and len(examples)<3: examples.append((tiles,conf[:2]))
        print(name, dict(sorted(C.items())),flush=True)
        for ex in examples[:1]: print('   example of a problem:', [t for t in ex[0] if t[0]=='O'], ex[1][:1])
