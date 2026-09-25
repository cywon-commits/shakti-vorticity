import sys
exec(open('viz_shortarm.py').read().split("fig,(A,B)=plt.subplots")[0].replace("n=5","n=int(sys.argv[1])").replace("range(20000)","range(0)"))
import networkx as nx
cl=sorted(cells,key=lambda c:(c[1],c[0]))
res=[];cnt=0
def rec(used,ds):
    global cnt
    free=next((c for c in cl if c not in used),None)
    if free is None:
        cnt+=1
        short={short_cell(tuple(sorted(d))) for d in ds}
        G=nx.Graph();G.add_nodes_from(short)
        for c in short:
            for dx,dy in [(1,0),(0,1)]:
                if (c[0]+dx,c[1]+dy) in short: G.add_edge(c,(c[0]+dx,c[1]+dy))
        res.append((nx.is_bipartite(G),tuple(sorted(len(c) for c in nx.cycle_basis(G)))));return
    for d in [(1,0),(0,1)]:
        o=(free[0]+d[0],free[1]+d[1])
        if o in cells and o not in used: rec(used|{free,o},ds+[(free,o)])
rec(frozenset(),[])
from collections import Counter; print(cnt,'all bipartite:',all(r[0] for r in res),'cycle lengths:',Counter(r[1] for r in res))
