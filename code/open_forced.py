"""Energy cost of forcing every plaquette of a block-aligned box into a tetromino (vs free edge vortices)."""
import sys, pickle
from collections import Counter
import open_ilp as O
from exact_v import aztec_cells
out={}
for n in map(int,sys.argv[1:]):
    A,_=aztec_cells(n,'II'); pts=[p for c in A for p in O.CELLS[c]]
    x0=min(p[0] for p in pts); x1=max(p[0] for p in pts)+1; y0=min(p[1] for p in pts); y1=max(p[1] for p in pts)+1
    box=O.shakti_rect(x0,y0,x1-x0,y1-y0)
    E0,l0=O.solve(x0,y0,x1-x0,y1-y0); g0=dict(O.solve.last)
    E1,l1=O.solve(x0,y0,x1-x0,y1-y0,force=box); g1=dict(O.solve.last)
    c1=Counter(l for l in l1.values() if O.COST[l]>0)
    print(f'n={n}: free edges E={E0:.0f}; all plaquettes forced E={E1:.0f} {dict(c1)}; difference {E1-E0:.0f} (n+1={n+1}); box plaquettes {len(box)}',flush=True)
    out[n]=dict(box=box,A=A,bbox=(x0,y0,x1,y1),free=dict(isl=g0['isl'],spins=g0['spins'],labs=l0,E=E0),forced=dict(isl=g1['isl'],spins=g1['spins'],labs=l1,E=E1))
pickle.dump(out,open('exact/open_forced.pkl','wb'))
