import random, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect
exec(open('pf_fig1.py').read().split('GAP=0.13')[0].split('from pubstyle import *')[0])
import itertools
n=5
base={(i,j) for i in range(-n,n) for j in range(-n,n) if abs(i+.5)+abs(j+.5)<=n}
src=open('pf_fig1.py').read(); s=src.index('dom={}'); e=src.index('typecol=')
exec(src[s:e])
def orient(r):
    xs={p[0] for p in r}; ys={p[1] for p in r}; return 'H' if len(xs)==2 else 'V'
for pl,sh in (('I',((n+1)%2,0)),('II',(n%2,0))):
    bad_shade=0; bad_perp=0; bad_adj=0; occ={}
    for d in doms0:
        dd=[(c[0]+sh[0],c[1]+sh[1]) for c in d]
        rs=[list(rect(c)) for c in dd]
        for c,r in zip(dd,rs):
            if (orient(r)=='H')!=black(c): bad_shade+=1
            for q in r: occ[q]=occ.get(q,0)+1
        if orient(rs[0])==orient(rs[1]): bad_perp+=1
        a,b=set(map(tuple,rs[0])),set(map(tuple,rs[1]))
        if not any((x+dx,y+dy) in b for (x,y) in a for dx,dy in ((1,0),(-1,0),(0,1),(0,-1))): bad_adj+=1
    print(pl,'plaquettes whose shade != orientation:',bad_shade,' non-perpendicular pairs:',bad_perp,' non-adjacent pairs:',bad_adj,' overlaps:',sum(v>1 for v in occ.values()))
    # black cell of (a) -> orientation in this placement
    cnt={}
    for d in doms0:
        for c in d:
            r=list(rect((c[0]+sh[0],c[1]+sh[1]))); k=('black(a)' if black(c) else 'white(a)',orient(r)); cnt[k]=cnt.get(k,0)+1
    print('   (a)-colour -> plaquette orientation:',cnt)
