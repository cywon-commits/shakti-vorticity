import sys, math
from collections import Counter
exec(open('viz_shortarm.py').read().split("fig,(A,B)=plt.subplots")[0].replace("n=5","n=int(sys.argv[1])").replace("range(20000)","range(0)"))
import networkx as nx
SX,SY=int(sys.argv[2]),int(sys.argv[3])
cells={(c[0]+SX,c[1]+SY) for c in cells}
def partner(c):
    return next(c2 for c2 in [(c[0]+1,c[1]+1),(c[0]-1,c[1]-1),(c[0]+1,c[1]-1),(c[0]-1,c[1]+1)] if inv[c2][:2]==inv[c][:2])
def pcen(c):
    r=rs[inv[c]]; return (sum(p[0] for p in r)/2+.5, sum(p[1] for p in r)/2+.5)
allsq={p for c in cells for p in rs[inv[c]]}
mx=sum(p[0]+.5 for p in allsq)/len(allsq); my=sum(p[1]+.5 for p in allsq)/len(allsq)
def touches_perimeter(c):              # Shakti definition: plaquette shares an edge with the outside
    return any((p[0]+dx,p[1]+dy) not in allsq for p in rs[inv[c]] for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)])
def outer_ring(c):                     # chessboard definition: cell has an orthogonal neighbour outside
    return any((c[0]+dx,c[1]+dy) not in cells for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)])
def graph(short, bdef):
    G=nx.Graph(); G.add_nodes_from(short)
    for c in short:
        for dx,dy in [(1,0),(0,1)]:
            if (c[0]+dx,c[1]+dy) in short: G.add_edge(c,(c[0]+dx,c[1]+dy),kind='int')
    bd=[c for c in short if bdef(c)]
    bd.sort(key=lambda c: math.atan2(pcen(c)[1]-my,pcen(c)[0]-mx))
    ring=[]
    if len(bd)>=2:
        for a,b in zip(bd,bd[1:]+bd[:1]):
            if a!=b and not G.has_edge(a,b): G.add_edge(a,b,kind='ring'); ring.append((a,b))
    return G,bd,ring
if __name__=='__main__':
    cl=sorted(cells,key=lambda c:(c[1],c[0]))
    stats={k:Counter() for k in ('shakti','chess')}
    def rec(used,ds):
        free=next((c for c in cl if c not in used),None)
        if free is None:
            short={short_cell(tuple(sorted(d))) for d in ds}
            for name,bdef in (('shakti',touches_perimeter),('chess',outer_ring)):
                G,bd,ring=graph(short,bdef)
                same=sum(black(a)==black(b) for a,b in ring)
                stats[name][(nx.is_bipartite(G), same%2)]+=1
            return
        for d in [(1,0),(0,1)]:
            o=(free[0]+d[0],free[1]+d[1])
            if o in cells and o not in used: rec(used|{free,o},ds+[(free,o)])
    rec(frozenset(),[])
    print('n=',n,'shift',(SX,SY), {k:dict(v) for k,v in stats.items()})
