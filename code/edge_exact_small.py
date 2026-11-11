"""Exact edge-pairing statistics for small n by summing over all balanced subsets S of edge plaquettes."""
import sys, math, itertools, io, contextlib, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    import open_ilp as O
from exact_v import aztec_cells, black, Region
n=int(sys.argv[1])
A,_=aztec_cells(n,'II'); pts=[p for c in A for p in O.CELLS[c]]
x0=min(p[0] for p in pts); x1=max(p[0] for p in pts)+1; y0=min(p[1] for p in pts); y1=max(p[1] for p in pts)+1
E=sorted(O.shakti_rect(x0,y0,x1-x0,y1-y0)-A); Eb=[e for e in E if black(e)]; Ew=[e for e in E if not black(e)]
L=[]   # (k=|S|, lnM)
for kb in range(len(Eb)+1):
    for Sb in itertools.combinations(Eb,kb):
        for Sw in itertools.combinations(Ew,kb):
            try: L.append((2*kb,Region(A|set(Sb)|set(Sw)).ld0))
            except Exception: pass
k=np.array([a for a,b in L]); lm=np.array([b for a,b in L])
for x in map(float,sys.argv[2:]):
    lw=lm-x*k/2; lw-=lw.max(); p=np.exp(lw); p/=p.sum()
    print(f'n={n} x={x:5.2f} exact paired fraction {np.sum(p*k)/len(E):.3f}')
