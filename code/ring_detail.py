import sys
exec(open('ring_shift.py').read().split("if __name__")[0])
cl=sorted(cells,key=lambda c:(c[1],c[0]))
C=Counter(); oddlen=Counter()
def rec(used,ds):
    free=next((c for c in cl if c not in used),None)
    if free is None:
        short={short_cell(tuple(sorted(d))) for d in ds}
        G,bd,ring=graph(short,touches_perimeter)
        same=sum(black(a)==black(b) for a,b in ring)
        C[(nx.is_bipartite(G),same)]+=1
        if not nx.is_bipartite(G):
            cyc=nx.find_cycle(G.subgraph(max(nx.connected_components(G),key=len)))
        return
    for d in [(1,0),(0,1)]:
        o=(free[0]+d[0],free[1]+d[1])
        if o in cells and o not in used: rec(used|{free,o},ds+[(free,o)])
rec(frozenset(),[])
print(sys.argv[1:], sorted(C.items()))
