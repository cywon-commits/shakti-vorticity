"""Check of the boundary-excitation model (edge_half.py) against the full open vertex model of straight-edged squares.
Vertex model: exact density of states g(E) by a transfer matrix (vertex_dos.py; all spin states, unhappy z=3 vertex J,
type-II and charged vertices 3J).  Model: Z = Pf(...) as an exact polynomial in t = e^{-J/T} (integer evaluation, python-flint).
Usage: python edge_vertex_check.py L1 L2 ...  (L = X, the side of the square; Aztec order n = L/2 - 1)"""
import sys, math, numpy as np, flint
from edge_half import square, EdgeModel
from vertex_dos import dos
from open_sample import Sample, halfplanes
def poly(E):
    """coefficients of Z in powers of t = e^{-J/T}, with dimer weight t, free-vortex weight 1, half-dimer weight n_c t"""
    B=E.B; N=B.N; k=3*N+80; t=1<<k
    A=[[0]*N for _ in range(N)]
    ii,jj=np.nonzero(np.triu(B.K,1))
    for i,j in zip(ii,jj): v=int(B.K[i,j])*t; A[i][j]=v; A[j][i]=-v
    ell=[int(E.isfree[a])+int(E.nc[a])*t for a in range(E.m)]; G=[int(g) for g in B.g]
    for a in range(E.m):
        for b in range(a+1,E.m):
            v=G[a]*G[b]*ell[a]*ell[b]; A[a][b]+=v; A[b][a]-=v
    d=int(flint.fmpz_mat(N,N,[v for r in A for v in r]).det()); z=math.isqrt(d); assert z*z==d
    c=[]
    while z: c.append(z&(t-1)); z>>=k
    return c
for L in [int(v) for v in sys.argv[1:]]:
    for kind in ('aligned','through'):
        cells,free,nh=square(L,kind); E=EdgeModel(cells,free,nh); c=poly(E)
        hp=halfplanes((1,0,0,L),(0,1,0,L)) if kind=='aligned' else halfplanes((1,0,1,L+1),(0,1,1,L+1))
        S=Sample(hp,W=L+10,mode='A'); E0=next(i for i,v in enumerate(c) if v); g,_,_=dos(S,E0+2)
        row=' | '.join(f'E0+{k}: vertex {g[E0+k]/2:.0f} model {c[E0+k] if E0+k<len(c) else 0}' for k in range(3))
        print(f'X={L} {kind:8s} n={L//2-1}: half-dimer sites {sum(nh.values())} | {row}',flush=True)
