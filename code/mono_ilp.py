"""Shakti vertex-model ILP with the unhappy-vertex (dimer) positions fixed by a tiling that may contain monomers."""
import itertools, random, pulp, networkx as nx
from shakti_torus import *
COST0={'z2ok':0,'z3happy':0,'z3unhappy':1,'z4I':0,'z4II':3,'z4q2':3,'z3q3':3,'z2q':3,'z4q4':20}
class Model:
    def __init__(s,P):
        s.T=T=Torus(P)
        s.adj=defaultdict(list)
        for c in T.cands: s.adj[c[0]].append(c); s.adj[c[1]].append(c)
        s.elb={}
        for c in T.cands:
            t=T.tmpl(c[2],c[3],ext=True); u=[]
            for p in T.verts:
                ks=[k for k,_,_,_ in T.inc[p]]
                if all(k in t for k in ks) and classify(*T.vstate(p,t))=='z3unhappy': u.append(p)
            assert len(u)==1,(u,c[:2]); s.elb[(c[0],c[1])]=u[0]
        # vertex state tables
        s.tab={}
        for p in T.verts:
            L=T.inc[p]; rows=[]
            for bits in itertools.product((1,-1),repeat=len(L)):
                sp={k:b for (k,_,_,_),b in zip(L,bits)}
                rows.append((bits,classify(*T.vstate(p,sp))))
            s.tab[p]=rows
    def tiling(s,removed=(),rng=None):
        T=s.T; G=nx.Graph(); keep=[i for i in range(len(T.plaqs)) if i not in removed]
        G.add_nodes_from(keep)
        cl=list(T.cands); 
        if rng: rng.shuffle(cl)
        for c in cl:
            if c[0] in removed or c[1] in removed: continue
            G.add_edge(c[0],c[1],c=c)
        top=[i for i in keep if T.horiz[i]]
        M=nx.bipartite.maximum_matching(G,top_nodes=top)
        doms=[G[i][j]['c'] for i,j in M.items() if i in top]
        assert 2*len(doms)==len(keep),'no perfect matching'
        return doms
    def solve(s,doms,cost=COST0,forced_unhappy=True,fix=None):
        T=s.T; prob=pulp.LpProblem('shakti',pulp.LpMinimize)
        y={}; obj=[]
        el={s.elb[(c[0],c[1])] for c in doms}
        for p,rows in s.tab.items():
            for a,(bits,lab) in enumerate(rows):
                y[p,a]=pulp.LpVariable(f'y_{p[0]}_{p[1]}_{a}',cat='Binary'); obj.append(cost[lab]*y[p,a])
            prob+= pulp.lpSum(y[p,a] for a in range(len(rows)))==1
            if forced_unhappy and p in el:
                prob+= pulp.lpSum(y[p,a] for a,(b,l) in enumerate(rows) if l=='z3unhappy')==1
        # island consistency: s_k in {0,1}
        sv={k:pulp.LpVariable(f's_{k}',cat='Binary') for k in range(len(T.isl))}
        for p,rows in s.tab.items():
            for j,(k,_,_,_) in enumerate(T.inc[p]):
                prob+= sv[k]==pulp.lpSum(y[p,a] for a,(bits,l) in enumerate(rows) if bits[j]==1)
        if fix:
            for p,labs in fix.items():
                prob+= pulp.lpSum(y[p,a] for a,(b,l) in enumerate(s.tab[p]) if l in labs)==1
        prob+= pulp.lpSum(obj)
        prob.solve(pulp.PULP_CBC_CMD(msg=0,timeLimit=300))
        E=pulp.value(prob.objective)
        spins={k:(1 if sv[k].value()>.5 else -1) for k in sv}
        labs={p:classify(*T.vstate(p,spins)) for p in T.verts}
        return E,spins,labs,el
if __name__=='__main__':
    import sys
    P=int(sys.argv[1]) if len(sys.argv)>1 else 8
    M=Model(P); d=M.tiling()
    E,sp,labs,el=M.solve(d)
    print('P',P,'tiles',len(d),'E',E, Counter(labs.values()))
