"""Boundary-monomer dimer model at finite fugacity: exact partition function and local statistics by a Pfaffian.
Finite-temperature effective model of an open Shakti sample: every freeable edge plaquette may stay unpaired with weight
y = exp(J/2T) (stem-island rule).  Monomers restricted to the outer face keep the model Pfaffian
(Giuliani, Jauslin, Lieb, J. Stat. Phys. 163, 211 (2016)):
    Z(y) = Pf(K + M),   M_ij = G_i G_j y^2  (i < j, freeable cells labelled in cyclic order along the boundary),
where G_i G_j is the Kasteleyn sign of a virtual edge drawn outside the sample between cells i and j."""
import io, contextlib, math, sys, numpy as np
from functools import lru_cache
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect, black
def centre(c):
    r=rect(c); return (sum(p[0] for p in r)/2+0.5, sum(p[1] for p in r)/2+0.5)
def outer_walk(cells):
    """closed walk along the outer face of the cell-adjacency graph (counterclockwise), as a list of cells"""
    cells=set(cells); start=min(cells,key=lambda c:(c[1],c[0])); dirs=[(1,0),(0,1),(-1,0),(0,-1)]
    def step(c,din):          # arriving at c moving in direction index din: take the first existing neighbour, scanning counterclockwise from the right-hand side
        for t in (-1,0,1,2):
            d=(din+t)%4; n=(c[0]+dirs[d][0],c[1]+dirs[d][1])
            if n in cells: return n,d
        raise ValueError('isolated cell')
    walk=[start]; c,d=step(start,0); first=(c,d)
    while True:
        walk.append(c); c,d=step(c,d)
        if (c,d)==first and walk[-1]==start: break
    return walk[:-1]
class BMono:
    def __init__(s,cells,free,c_sign=1):
        cen={c:centre(c) for c in cells}; walk=outer_walk(cells); seen={}
        for k,c in enumerate(walk):
            if c in free and c not in seen: seen[c]=k     # a cell visited twice (corner) keeps its first position
        assert len(seen)==len(free),'a freeable cell is not on the outer face'
        F=sorted(free,key=seen.get)
        s.order=F+sorted(set(cells)-set(free)); s.idx={c:k for k,c in enumerate(s.order)}; s.m=len(F); s.N=len(s.order); s.cen=cen
        K=np.zeros((s.N,s.N))
        for c in cells:
            for d,w in (((1,0),1.0),((0,1),-1.0 if c[0]%2 else 1.0)):
                n=(c[0]+d[0],c[1]+d[1])
                if n in s.idx: K[s.idx[c],s.idx[n]]=w; K[s.idx[n],s.idx[c]]=-w
        s.K=K; s.col=np.array([1 if black(c) else -1 for c in s.order])
        # G_i = product of the Kasteleyn signs along the boundary walk from freeable cell 0 to freeable cell i
        G=np.ones(s.m); sg=1.0; k0=seen[F[0]] if F else 0; j=0
        for k in range(k0,len(walk)-1):
            sg*=np.sign(K[s.idx[walk[k]],s.idx[walk[k+1]]])
            if walk[k+1] in seen and seen[walk[k+1]]==k+1: j+=1; G[j]=sg
        s.g=G; s.c=c_sign
        s.odd=s.N%2==1      # odd number of cells: one auxiliary boundary cell that is always a monomer (weight 1), last in the cyclic order
    def matrix(s,y):
        m=s.m; A=s.K.copy(); i,j=np.triu_indices(m,1)
        v=s.c*s.g[i]*s.g[j]*y*y; A[i,j]+=v; A[j,i]-=v
        return A
    def lnZ(s,y):
        sg,ld=np.linalg.slogdet(s.matrix(y)); assert sg>0; return ld/2
    def stats(s,y):
        """monomer probability of each freeable cell, pairing probabilities (dict), and per-cell max pairing probability"""
        A=s.matrix(y); Ai=np.linalg.inv(A); m=s.m
        M=A[:m,:m]-s.K[:m,:m]; mono=np.einsum('ij,ji->i',M,Ai[:m,:m])
        P=s.K*Ai.T; 
        return mono,P
    def exact(s,y):
        """Exact rational arithmetic (python-flint) for integer y: returns ln Z, monomer probabilities (freeable cells, in s.order),
        and the dict {(i,j): pairing probability} over adjacent cells (indices in s.order).  Handles an odd number of cells."""
        import flint
        m,N=s.m,s.N; n=N+1 if s.odd else N; A=[[0]*n for _ in range(n)]
        ii,jj=np.nonzero(np.triu(s.K,1))
        for i,j in zip(ii,jj): v=int(s.K[i,j]); A[i][j]=v; A[j][i]=-v
        lab=list(range(m))+([N] if s.odd else []); ell=[y]*m+([1] if s.odd else []); G=[int(g) for g in s.g]+([1] if s.odd else [])
        for a in range(len(lab)):
            for b in range(a+1,len(lab)):
                v=G[a]*G[b]*ell[a]*ell[b]; A[lab[a]][lab[b]]+=v; A[lab[b]][lab[a]]-=v
        M=flint.fmpz_mat(n,n,[v for r in A for v in r]); d=M.det()
        if d==0: return None
        assert d>0; lnZ=math.log(int(d))/2; Ai=flint.fmpq_mat(M).inv()
        f=lambda q: int(q.p)/int(q.q)
        pair={(int(i),int(j)): int(s.K[i,j])*f(Ai[int(j),int(i)]) for i,j in zip(ii,jj)}
        mono=np.zeros(m)
        for a in range(len(lab)):
            for b in range(a+1,len(lab)):
                p=G[a]*G[b]*ell[a]*ell[b]*f(Ai[lab[b],lab[a]])
                if a<m: mono[a]+=p
                if b<m: mono[b]+=p
        return lnZ,mono,pair
def brute(cells,free):
    """Z as a polynomial in y (coefficient list by number of monomers), by exhaustive recursion.  Small samples only."""
    cl=sorted(cells); idx={c:k for k,c in enumerate(cl)}; nb=[[idx[(c[0]+a,c[1]+b)] for a,b in ((1,0),(-1,0),(0,1),(0,-1)) if (c[0]+a,c[1]+b) in idx] for c in cl]
    fr=[c in free for c in cl]; n=len(cl)
    @lru_cache(None)
    def rec(mask):
        if mask==(1<<n)-1: return (1,)
        k=0
        while mask>>k&1: k+=1
        tot=[0]*(n+1)
        for j in nb[k]:
            if not mask>>j&1:
                for a,v in enumerate(rec(mask|1<<k|1<<j)): tot[a]+=v
        if fr[k]:
            for a,v in enumerate(rec(mask|1<<k)): tot[a+1]+=v
        while len(tot)>1 and tot[-1]==0: tot.pop()
        return tuple(tot)
    return rec(0)
if __name__=='__main__':
    from gen_ilp import Sample, halfplanes
    tests=[('block-aligned box 8',halfplanes((1,0,0,8),(0,1,0,8)),'A'),('block-aligned box 8x12',halfplanes((1,0,0,8),(0,1,0,12)),'A'),
           ('mixed box X block-aligned 8, Y through-block 6',halfplanes((1,0,0,8),(0,1,1,7)),'A'),('mixed box 12 x 6',halfplanes((1,0,0,12),(0,1,1,7)),'A'),
           ('through-block box 8, clean',halfplanes((1,0,1,9),(0,1,1,9)),'B'),('45-degree diamond D=4, clean',halfplanes((1,1,0,8),(1,-1,0,8)),'B'),
           ('mixed 8 x 9',halfplanes((1,0,0,8),(0,1,0,9)),'A'),('mixed 12 x 9 (odd number of cells)',halfplanes((1,0,0,12),(0,1,0,9)),'A'),('clean 9 x 9',halfplanes((1,0,1,10),(0,1,0,9)),'B'),
           ('plus n=6 c=2, clean',lambda p:(abs(p[0])<=6 and abs(p[1])<=2) or (abs(p[0])<=2 and abs(p[1])<=6),'B')]
    for name,ins,mode in tests:
        S=Sample(ins,W=14,mode=mode)
        if len(S.cells)>44: print(name,'too large for brute force',len(S.cells)); continue
        poly=brute(S.cells,S.free)
        print(f'{name}: cells {len(S.cells)}, freeable {len(S.free)}, Z(y) coefficients {poly}')
        if sum(poly)==0: continue
        for c in (1,-1):
            try:
                B=BMono(S.cells,S.free,c)
                if not B.odd:
                    d=[B.lnZ(y)-math.log(sum(v*y**a for a,v in enumerate(poly))) for y in (0.3,1.0,2.5)]
                    print(f'   c={c:+d}: ln Z(Pfaffian) - ln Z(brute) at y = 0.3, 1, 2.5: {np.round(d,10)}')
                if c==1:
                    for y in (1,3):
                        r=B.exact(y); zb=sum(v*y**a for a,v in enumerate(poly)); nb=sum(a*v*y**a for a,v in enumerate(poly))/zb
                        print(f'   exact, y={y}: Z {math.exp(r[0]):.6f} (brute {zb}), <monomers> {r[1].sum():.6f} (brute {nb:.6f})')
            except AssertionError as e: print(f'   c={c:+d}: failed ({e})')
