import io, contextlib
from collections import Counter, defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    from otet import *
from vform2 import S_TAB, r4
def height_mod8(tiles,cells):
    part={}
    for t in tiles:
        if t[0]=='D': part[t[1]]=t[2]; part[t[2]]=t[1]
    V=set()
    for (x,y) in cells: V|={(x,y),(x+1,y),(x,y+1),(x+1,y+1)}
    for t in tiles:     # the corner shared by the two cells of a diagonal tile lies inside the defect
        if t[0]=='O':
            c1,c2=t[1],t[2]; V.discard((max(c1[0],c2[0]),max(c1[1],c2[1])))
    start=min(V); h={start:r4(start)}; Q=[start]; bad=0
    for p in Q:
        for d in [(1,0),(-1,0),(0,1),(0,-1)]:
            q=(p[0]+d[0],p[1]+d[1])
            if q not in V: continue
            if d==(1,0): L,R=(p[0],p[1]),(p[0],p[1]-1)
            elif d==(-1,0): L,R=(q[0],q[1]-1),(q[0],q[1])
            elif d==(0,1): L,R=(p[0]-1,p[1]),(p[0],p[1])
            else: L,R=(q[0],q[1]),(q[0]-1,q[1])
            if L not in cells and R not in cells: continue
            if L in part and part[L]==R: continue
            val=h[p]+(1 if black(L) else -1)
            if q not in h: h[q]=val; Q.append(q)
            elif (h[q]-val)%8!=0: bad+=1
    return h,bad
def corner_height(d,h):
    # walk the domino perimeter counter-clockwise from p; the domino lies on the left of every step
    cells_=set(d); p=d[0]
    if d[1][0]!=d[0][0]: cyc=[p,(p[0]+1,p[1]),(p[0]+2,p[1]),(p[0]+2,p[1]+1),(p[0]+1,p[1]+1),(p[0],p[1]+1)]
    else: cyc=[p,(p[0]+1,p[1]),(p[0]+1,p[1]+1),(p[0]+1,p[1]+2),(p[0],p[1]+2),(p[0],p[1]+1)]
    inc=[]
    for a,b in zip(cyc,cyc[1:]+cyc[:1]):
        dx,dy=b[0]-a[0],b[1]-a[1]
        left={(1,0):(a[0],a[1]),(-1,0):(b[0],b[1]-1),(0,1):(a[0]-1,a[1]),(0,-1):(b[0],b[1])}[(dx,dy)]
        inc.append(1 if black(left) else -1)
    for i,q in enumerate(cyc):
        if q in h:
            # go backwards from q to p along the cycle
            val=h[q]
            for j in range(i-1,-1,-1): val-=inc[j]
            return val
    raise KeyError('no corner')
def eps(hv,p): return 1 if ((hv-r4(p))//4)%2==0 else -1
res=Counter(); tabO=defaultdict(Counter); badh=0
for name,cells in [('Aztec n=4 (I)',aztec(4)[0]),('6x6',{(i,j) for i in range(6) for j in range(6)}),('Aztec n=3 (II)',aztec(3,(1,0))[0])]:
    for tiles in enum_with_defects(cells):
        if not any(t[0]=='O' for t in tiles): continue
        conf,fr,comps,v=analyse(tiles); h,bad=height_mod8(tiles,cells); badh+=bad
        # dominoes: v = s(d) eps(h(p)) up to a global sign
        prods=set(); skip=False
        for k,t in enumerate(tiles):
            if t[0]=='D':
                d=(t[1],t[2]); p=d[0]; typ='H' if d[1][0]!=d[0][0] else 'V'
                try: hp=corner_height(d,h)
                except KeyError: skip=True; break
                prods.add(v[k]*S_TAB[(typ,p[0]%4,p[1]%4)]*eps(hp,p))
        if skip: res[(name,'enclosed region, skipped')]+=1; continue
        res[(name,len(prods),comps)]+=1
        if len(prods)==1:
            g=prods.pop()
            for k,t in enumerate(tiles):
                if t[0]=='O':
                    c1,c2=t[1],t[2]; p=(min(c1[0],c2[0]),min(c1[1],c2[1]))   # lower-left corner of the diagonal pair's bounding box
                    horiz=len({q[1] for q in rect(c1)})==1
                    for off in [(0,0),(1,0),(2,0),(2,1),(2,2),(1,2),(0,2),(0,1)]:
                        q=(p[0]+off[0],p[1]+off[1])
                        if q in h: tabO[('H' if horiz else 'V',p[0]%4,p[1]%4,off)][v[k]*g*eps(h[q],q)]+=1
print('height mod-8 inconsistencies:',badh)
print('domino formula holds in coverings with defects (prods==1):',dict(res))
from collections import defaultdict as dd
byoff=dd(lambda:[0,0])
for k,c in tabO.items():
    byoff[(k[0],k[3])][0]+=1; byoff[(k[0],k[3])][1]+= (len(c)>1)
print('reference vertex (orientation, offset in the bounding box): [classes, ambiguous]')
for k in sorted(byoff): print('  ',k,byoff[k])
good=[k for k in byoff if byoff[k][1]==0]
for g in good[:2]:
    print('table for',g,{kk[1:3]:list(c.keys())[0] for kk,c in tabO.items() if (kk[0],kk[3])==g})
