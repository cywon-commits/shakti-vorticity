"""Tiling language: coupling-graph frustration with monomers (uncovered chess cells) vs O defects."""
import sys, random, math, io, contextlib, pulp, networkx as nx
from collections import Counter, defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    import otet
from otet import tile_template, partner, aztec
def couplings(tiles):
    owner=defaultdict(list)
    for k,t in enumerate(tiles):
        for e,(b,kind) in tile_template(t).items(): owner[e].append((k,b))
    J={}; conf=0
    for e,lst in owner.items():
        for x in range(len(lst)):
            for y in range(x+1,len(lst)):
                (i,bi),(j,bj)=lst[x],lst[y]
                if i==j: continue
                k=(min(i,j),max(i,j))
                if k in J and J[k]!=bi*bj: conf+=1
                J[k]=bi*bj
    return J,conf
def min_frustration(n,J):
    prob=pulp.LpProblem('f',pulp.LpMinimize)
    u=[pulp.LpVariable(f'u{i}',cat='Binary') for i in range(n)]
    w={k:pulp.LpVariable(f'w{k[0]}_{k[1]}',lowBound=0) for k in J}
    for (i,j),s in J.items():
        if s==1: prob+= w[i,j]>=u[i]-u[j]; prob+= w[i,j]>=u[j]-u[i]
        else: prob+= w[i,j]>=u[i]+u[j]-1; prob+= w[i,j]>=1-u[i]-u[j]
    prob+= pulp.lpSum(w.values()); prob.solve(pulp.PULP_CBC_CMD(msg=0))
    viol=[k for k in J if w[k].value()>.5]
    return round(pulp.value(prob.objective)), viol, [int(x.value()>.5) for x in u]
def random_tiling(cells,rng,allow_O=()):
    G=nx.Graph(); cl=sorted(cells); rng.shuffle(cl); G.add_nodes_from(cl)
    for c in cl:
        for d in [(1,0),(0,1),(-1,0),(0,-1)]:
            o=(c[0]+d[0],c[1]+d[1])
            if o in cells: G.add_edge(c,o,w=rng.random())
    M=nx.max_weight_matching(G,maxcardinality=True,weight='w')
    assert 2*len(M)==len(cells)
    return [('D',)+tuple(sorted(e)) for e in M]
if __name__=='__main__':
    L=int(sys.argv[1]); rng=random.Random(0)
    board={(i,j) for i in range(L) for j in range(L)}
    a=(L//2-1,L//2-1)
    print('full board: min frustration',min_frustration(len(random_tiling(board,rng)),couplings(random_tiling(board,rng))[0])[0])
    rows=[]
    for dx in range(1,L//2+1,1):
        for dy in (0,):
            b=(a[0]+dx,a[1]+dy)
            if (sum(a)+sum(b))%2==0: b=(b[0],b[1]+1)
            if b not in board: continue
            fs=[]
            for t in range(6):
                tiles=random_tiling(board-{a,b},rng); J,conf=couplings(tiles)
                f,viol,u=min_frustration(len(tiles),J); fs.append(f)
            d=abs(a[0]-b[0])+abs(a[1]-b[1])
            print(f'monomers {a},{b} (chess distance {d}): min # violated couplings over 6 tilings = {fs}',flush=True)
