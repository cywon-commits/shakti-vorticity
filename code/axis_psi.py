"""Straight-edged X x X Shakti samples. Edges along block boundaries: the sample is the placement-II Aztec diamond of
order n = X/2-1 plus 2(n+1) edge plaquettes that close their flux on their own (open_vortex.py); edges through the
blocks: exactly the placement-I Aztec diamond. Exact Psi over the tiled region; key "forced" = block-aligned box
with all plaquettes in tetrominoes (comparison). Usage: python axis_psi.py 48"""
import sys, pickle, time
from axis_region import CELLS, shakti_rect
from exact_psi import psi_map
from exact_v import Region, start_vertex, aztec_cells
X=int(sys.argv[1]); n=X//2-1; out={}
for pl in ('II','I'):
    t=time.time(); A,_=aztec_cells(n,pl); pts=[p for c in A for p in CELLS[c]]
    x0=min(p[0] for p in pts); y0=min(p[1] for p in pts); box=shakti_rect(x0,y0,X,X)
    R=Region(A); s,hs=start_vertex(A); P,E=psi_map(R,s,hs,sorted(A))
    out[pl]=dict(cells=A,box=box,P=P,origin=(x0,y0),X=X,n=n)
    print(pl,'n',n,'Aztec cells',len(A),'complete plaquettes in box',len(box),'%.0fs'%(time.time()-t),flush=True)
# comparison: the same block-aligned box with every plaquette forced into a tetromino (no free edge vortices)
B=out['II']['box']; R=Region(B); s,hs=start_vertex(B); P,E=psi_map(R,s,hs,sorted(B))
out['forced']=dict(cells=B,box=B,P=P,origin=out['II']['origin'],X=X,n=n)
print('forced: tiled plaquettes',len(B),flush=True)
pickle.dump(out,open(f'exact/axis_psi_{X}.pkl','wb'))
