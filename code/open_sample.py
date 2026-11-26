# Copied from the companion repository github.com/cywon-commits/shakti-arctic-design (code/gen_ilp.py, same authors): open Shakti samples and vertex energies.
"""General open Shakti sample (any polygonal cut): exact vertex-model ground state (ILP) vs effective boundary-monomer dimer model (LP).
Run from inside code/ (core.py reads embed.py from the working directory)."""
import io, contextlib, itertools, math, sys
from collections import defaultdict, Counter
import numpy as np, pulp
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
with contextlib.redirect_stdout(io.StringIO()):
    from vort import is_island
    from core import rect, black
from shakti_torus import classify
COST={'z1':0,'z2ok':0,'z3happy':0,'z3unhappy':1,'z4I':0,'z4II':3,'z4q2':3,'z3q3':3,'z2q':3,'z4q4':20}
_rc={}
def R(c):
    if c not in _rc: _rc[c]=tuple(sorted(rect(c)))
    return _rc[c]
def ppts(c): return {(p[0]+a,p[1]+b) for p in R(c) for a in (0,1) for b in (0,1)}
def perim(c):
    r=R(c); xs=[p[0] for p in r]; ys=[p[1] for p in r]; X0,X1,Y0,Y1=min(xs),max(xs)+1,min(ys),max(ys)+1
    E=[]
    for x in range(X0,X1): E+= [((x,Y0),(x+1,Y0)),((x,Y1),(x+1,Y1))]
    for y in range(Y0,Y1): E+= [((X0,y),(X0,y+1)),((X1,y),(X1,y+1))]
    return E
def stem(c):
    """the island that makes the long-side midpoint of plaquette c a z=3 vertex (exactly one exists in the bulk)"""
    r=R(c); xs=sorted({p[0] for p in r}); ys=sorted({p[1] for p in r})
    if len(xs)==2:   # horizontal
        x=xs[1]; y=ys[0]; cand=[((x,y+1),(x,y+2)),((x,y-1),(x,y))]
    else:
        x=xs[0]; y=ys[1]; cand=[((x+1,y),(x+2,y)),((x-1,y),(x,y))]
    s=[e for e in cand if is_island(e)]; assert len(s)==1,(c,s); return s[0]
class Sample:
    def __init__(s, inside, W=40, mode='A'):
        """inside(point)->bool on lattice points. mode A: every island with both ends inside. mode B: only perimeters of complete plaquettes."""
        s.cells={(i,j) for i in range(-W,W) for j in range(-W,W) if all(inside(p) for p in ppts((i,j)))}
        if mode=='A':
            P=[(x,y) for x in range(-2*W,2*W) for y in range(-2*W,2*W) if inside((x,y))]; Ps=set(P)
            s.isl=[e for (x,y) in P for e in (((x,y),(x+1,y)),((x,y),(x,y+1))) if e[1] in Ps and is_island(e)]
        else:
            s.isl=sorted({e for c in s.cells for e in perim(c)})
        s.islset=set(s.isl)
        s.free={c for c in s.cells if stem(c) not in s.islset}
    def ilp(s, tl=900):
        inc=defaultdict(list)
        for k,(a,b) in enumerate(s.isl): inc[a].append((k,'start',a,b)); inc[b].append((k,'end',a,b))
        prob=pulp.LpProblem('open',pulp.LpMinimize); sv={k:pulp.LpVariable(f's{k}',cat='Binary') for k in range(len(s.isl))}
        obj=[]
        def lab(L,bits):
            ins=[];dirs=[]
            for (k,role,a,b),bt in zip(L,bits):
                ins.append((bt==1 and role=='end') or (bt==-1 and role=='start')); dirs.append((b[0]-a[0],b[1]-a[1]))
            z=len(L); q=sum(1 if i else -1 for i in ins)
            return 'z1' if z==1 else classify(z,q,ins,dirs)
        n=0
        for p,L in inc.items():
            rows=[(bits,lab(L,bits)) for bits in itertools.product((1,-1),repeat=len(L))]
            ys=[pulp.LpVariable(f'y{n}_{a}',cat='Binary') for a in range(len(rows))]; n+=1
            prob+=pulp.lpSum(ys)==1; obj+=[COST[l]*v for (b,l),v in zip(rows,ys)]
            for j,(k,_,_,_) in enumerate(L): prob+= sv[k]==pulp.lpSum(v for (b,l),v in zip(rows,ys) if b[j]==1)
        prob+=pulp.lpSum(obj); prob.solve(pulp.PULP_CBC_CMD(msg=0,timeLimit=tl))
        spins={k:(1 if sv[k].value()>.5 else -1) for k in sv}
        labs=Counter(lab(inc[p],[spins[k] for k,_,_,_ in inc[p]]) for p in inc)
        return round(pulp.value(prob.objective)), {l:c for l,c in labs.items() if COST[l]>0}, pulp.LpStatus[prob.status]
    def edges(s):
        return [(c,n) for c in s.cells if black(c) for n in [(c[0]+1,c[1]),(c[0]-1,c[1]),(c[0],c[1]+1),(c[0],c[1]-1)] if n in s.cells]
    def lp(s, free=None, maximize_free=True):
        """min number of dimers s.t. every non-freeable cell is matched, freeable cells matched at most once. returns (#dimers, #free cells used) or None"""
        free=s.free if free is None else free
        E=s.edges(); cl=sorted(s.cells); idx={c:k for k,c in enumerate(cl)}
        if not cl: return (0,0,[])
        A=lil_matrix((len(cl),len(E)))
        for k,(a,b) in enumerate(E): A[idx[a],k]=1; A[idx[b],k]=1
        A=A.tocsr(); eq=[k for k,c in enumerate(cl) if c not in free]; le=[k for k,c in enumerate(cl) if c in free]
        res=linprog(np.ones(len(E)),A_ub=A[le] if le else None,b_ub=np.ones(len(le)) if le else None,A_eq=A[eq] if eq else None,b_eq=np.ones(len(eq)) if eq else None,bounds=(0,1),method='highs')
        if res.status!=0: return None
        x=res.x; assert np.allclose(x,np.round(x),atol=1e-6),'fractional LP'
        dim=[E[k] for k in range(len(E)) if x[k]>.5]
        return (len(dim), len(cl)-2*len(dim), dim)
def halfplanes(*hp):
    """each hp = (a,b,lo,hi): lo <= a*x+b*y <= hi"""
    return lambda p: all(lo<=a*p[0]+b*p[1]<=hi for a,b,lo,hi in hp)
