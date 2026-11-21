"""Pfaffian for monomers restricted to the edge plaquettes (all on the outer face)."""
import sys, math, itertools, io, contextlib, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    import open_ilp as O
from exact_v import aztec_cells, black, Region
def geometry(n):
    A,_=aztec_cells(n,'II'); pts=[p for c in A for p in O.CELLS[c]]
    x0=min(p[0] for p in pts); x1=max(p[0] for p in pts)+1; y0=min(p[1] for p in pts); y1=max(p[1] for p in pts)+1
    box=O.shakti_rect(x0,y0,x1-x0,y1-y0); return A,box,sorted(box-A)
def skew(cells):
    cells=sorted(cells); idx={c:k for k,c in enumerate(cells)}; N=len(cells); M=np.zeros((N,N))
    for (x,y),k in idx.items():
        q=(x+1,y)
        if q in idx: M[k,idx[q]]=1; M[idx[q],k]=-1
        q=(x,y+1)
        if q in idx:
            s=1 if x%2==0 else -1; M[k,idx[q]]=s; M[idx[q],k]=-s
    return cells,idx,M
if __name__=='__main__':
    n=int(sys.argv[1]); A,box,E=geometry(n)
    cells,idx,K=skew(box)
    cx=np.mean([c[0] for c in box]); cy=np.mean([c[1] for c in box])
    Eo=sorted(E,key=lambda c: math.atan2(c[1]+.5-cy,c[0]+.5-cx))   # counterclockwise
    # brute force reference Z(z) = sum_S z^{#free} M(A u S)
    Eb=[e for e in E if black(e)]; Ew=[e for e in E if not black(e)]
    terms=[]
    for kb in range(len(Eb)+1):
        for Sb in itertools.combinations(Eb,kb):
            for Sw in itertools.combinations(Ew,kb):
                try: R=Region(A|set(Sb)|set(Sw)); terms.append((len(E)-2*kb,R.ld0))
                except Exception: pass
    def Zref(z): return sum(math.exp(l)*z**f for f,l in terms)
    zt=1.7; ref=Zref(zt)
    ie=[idx[e] for e in Eo]; best=None
    for sig in itertools.product((1,-1),repeat=len(E)-1):
        sg=(1,)+sig; Mx=K.copy()
        for a in range(len(ie)):
            for b in range(a+1,len(ie)):
                v=sg[a]*sg[b]*zt*zt; Mx[ie[a],ie[b]]+=v; Mx[ie[b],ie[a]]-=v
        val=math.sqrt(abs(np.linalg.det(Mx)))
        if abs(val/ref-1)<1e-9: best=sg; print('match',sg); break
    print('ref',ref,'best',best)
