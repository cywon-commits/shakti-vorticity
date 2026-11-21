"""Open vertex model: cost of pairing a subset S of edge plaquettes (n=5 block-aligned box). Expect E = n(n+1) + |S|/2."""
import io, contextlib, random, sys
with contextlib.redirect_stdout(io.StringIO()):
    import open_ilp as O
from exact_v import aztec_cells, black, Region
n=int(sys.argv[1]) if len(sys.argv)>1 else 5
A,_=aztec_cells(n,'II'); pts=[p for c in A for p in O.CELLS[c]]
x0=min(p[0] for p in pts); x1=max(p[0] for p in pts)+1; y0=min(p[1] for p in pts); y1=max(p[1] for p in pts)+1
box=O.shakti_rect(x0,y0,x1-x0,y1-y0); E=sorted(box-A); Eb=[e for e in E if black(e)]; Ew=[e for e in E if not black(e)]
random.seed(1); tests=[]
for k in (1,2,3,4,5):
    for rep in range(2):
        S=set(random.sample(Eb,k))|set(random.sample(Ew,k)); tests.append(S)
for S in tests:
    try: Region(A|S)
    except Exception: print('|S|=',len(S),'not tileable'); continue
    En,labs=O.solve(x0,y0,x1-x0,y1-y0,force=A|S)
    print(f'|S|={len(S):2d}  E_ILP={En:.0f}  n(n+1)+|S|/2={n*(n+1)+len(S)//2}',flush=True)
