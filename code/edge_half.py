"""Boundary excitations of straight-edged Shakti samples at finite temperature.
Two kinds of boundary monomer, both on the outer face (so the model stays Pfaffian):
  * freeable edge plaquettes (block-aligned edges): a free plaquette closes into a complete vortex; weight y = e^{J/2T}
    relative to a paired one (stem-island rule);
  * half-dimer sites: a boundary z=3 vertex that touches a single plaquette (or two non-adjacent ones) can be unhappy while
    covering only that plaquette.  Cost J/2 relative to a dimer, weight n_c e^{-J/2T}, with n_c the number of such vertices
    of the plaquette counted with their spin multiplicity (2 for a boundary point of coordination 3 in the full lattice, 1 for a
    corner point; for a point shared by two non-adjacent plaquettes, 1 for each non-freeable one).
The weights were fixed by complete enumeration of the first excited level of the open vertex model and are checked by
edge_vertex_check.py (exact density of states, transfer matrix).  Normalisation: dimers have weight 1."""
import io, contextlib, math, numpy as np
from collections import defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    from vort import is_island
    from open_sample import Sample, halfplanes, ppts
from bmono import BMono
def fulldeg(p):
    x,y=p; return sum(is_island(e) for e in (((x,y),(x+1,y)),((x-1,y),(x,y)),((x,y),(x,y+1)),((x,y-1),(x,y))))
def halfsites(S):
    deg=defaultdict(int)
    for a,b in S.isl: deg[a]+=1; deg[b]+=1
    touch=defaultdict(list)
    for c in S.cells:
        for p in ppts(c): touch[p].append(c)
    n=defaultdict(int)
    for p,d in deg.items():
        if d!=3: continue
        t=touch[p]
        if len(t)==1: n[t[0]]+=2 if fulldeg(p)==3 else 1
        elif len(t)==2 and abs(t[0][0]-t[1][0])+abs(t[0][1]-t[1][1])!=1:
            for c in t:
                if c not in S.free: n[c]+=1
    return dict(n)
def square(L,kind):
    """X x X Shakti square, X = L: kind 'aligned' (edges along block boundaries, Aztec II of order L/2-1 plus 2(n+1) edge
    plaquettes) or 'through' (edges through the blocks, fringed: Aztec I of order L/2-1)."""
    hp=halfplanes((1,0,0,L),(0,1,0,L)) if kind=='aligned' else halfplanes((1,0,1,L+1),(0,1,1,L+1))
    S=Sample(hp,W=L+10,mode='A'); return set(S.cells),set(S.free),halfsites(S)
class EdgeModel:
    def __init__(s,cells,free,nh):
        s.cells=cells; s.free=free; s.nh=nh; s.sites=set(free)|set(nh)
        s.B=BMono(cells,s.sites); assert not s.B.odd
        s.idx=s.B.idx; s.m=s.B.m; s.order=s.B.order
        s.isfree=np.array([c in free for c in s.order[:s.m]],float); s.nc=np.array([nh.get(c,0) for c in s.order[:s.m]],float)
    def matrix(s,x):
        """x = J/T"""
        ell=s.isfree*math.exp(x/2)+s.nc*math.exp(-x/2); m=s.m; A=s.B.K.copy(); i,j=np.triu_indices(m,1)
        v=s.B.g[i]*s.B.g[j]*ell[i]*ell[j]; A[i,j]+=v; A[j,i]-=v; return A,ell
    def stats(s,x):
        """returns M, M^{-1}, monomer probability of each site split into (free-vortex part, half-dimer part)"""
        A,ell=s.matrix(x); Ai=np.linalg.inv(A); m=s.m
        mono=np.einsum('ij,ji->i',A[:m,:m]-s.B.K[:m,:m],Ai[:m,:m])
        fr=np.where(ell>0,s.isfree*math.exp(x/2)/np.where(ell>0,ell,1),0)
        return A,Ai,mono*fr,mono*(1-fr)
