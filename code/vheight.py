import sys
from vort_tiling import *
from stats_run import enum
def height(doms,cells):
    part={}
    for a,b in doms: part[a]=b; part[b]=a
    V=set()
    for (x,y) in cells: V|={(x,y),(x+1,y),(x,y+1),(x+1,y+1)}
    start=min(V); h={start:0}; Q=[start]
    for p in Q:
        for d in [(1,0),(-1,0),(0,1),(0,-1)]:
            q=(p[0]+d[0],p[1]+d[1])
            if q not in V: continue
            # squares on the two sides of edge p->q
            if d==(1,0): left,right=(p[0],p[1]),(p[0],p[1]-1)
            elif d==(-1,0): left,right=(q[0],q[1]-1),(q[0],q[1])
            elif d==(0,1): left,right=(p[0]-1,p[1]),(p[0],p[1])
            else: left,right=(q[0],q[1]),(q[0]-1,q[1])
            if left not in cells and right not in cells: continue
            if left in cells and right in cells and part[left]==right: continue   # crosses a domino
            dh=1 if black(left) else -1
            if q not in h: h[q]=h[p]+dh; Q.append(q)
            else: assert h[q]==h[p]+dh
    return h
if __name__=='__main__':
    n=int(sys.argv[1]); it,cells=enum(n)
    data=defaultdict(lambda: defaultdict(Counter))
    for doms in it:
        I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
        spin,free,ring,_=analyse(doms,cells)
        sc={short_cell(*doms[k]):k for k in range(len(doms))}; ref=v[sc[ring[0]]]
        h=height(doms,cells)
        for k,(a,b) in enumerate(doms):
            d=tuple(sorted((a,b)))
            verts=sorted({(x+dx,y+dy) for (x,y) in d for dx in (0,1) for dy in (0,1)})
            for vi,p in enumerate(verts): data[d][vi][(h[p]%8, v[k]*ref)]+=1
    # for each placement, is there a vertex whose h mod 8 determines v?
    ok=0; tot=0; which=Counter()
    for d,dv in data.items():
        tot+=1; good=[]
        for vi,c in dv.items():
            m=defaultdict(set)
            for (hm,vv),cnt in c.items(): m[hm].add(vv)
            if all(len(s)==1 for s in m.values()): good.append(vi)
        ok+=bool(good); which[tuple(good)]+=1
    print('n',n,'placements',tot,'determined by h mod 8 at some vertex:',ok, which.most_common(5))
