import math, itertools
from exact_v import *
from vform2 import height_norm, v_formula
from placement2 import enum_cells
import io, contextlib
n=4; cells,sh=aztec_cells(n); R=Region(cells)
print('ln det K =',R.ld0,' expected ln 2^(n(n+1)/2) =',n*(n+1)/2*math.log(2))
start,hs=start_vertex(cells)
# brute force
tot=0; E={}; Vc={}
with contextlib.redirect_stdout(io.StringIO()):
    tilings=list(enum_cells(cells))
for doms in tilings:
    h=height_norm(doms,cells); c0=h[start]-hs
    for p in [(sh[0],0),(sh[0]+1,1),(sh[0]-1,2),(sh[0]+2,-1)]:
        e=1 if (((h[p]-c0)-r4(p))//4)%2==0 else -1; E[p]=E.get(p,0)+e
    for d in doms:
        v=v_formula(d,{q:x-c0 for q,x in h.items()})
        for c in d: Vc[c]=Vc.get(c,0)+v
N=len(tilings)
for p in E: print('vertex',p,' brute <eps> = %+.5f   exact = %+.5f'%(E[p]/N,R.expect_eps(start,hs,p)))
Vx=v_map(R,start,hs,[(sh[0]+i,j) for i in (-1,0,1) for j in (-1,0,1)])
for c,val in Vx.items(): print('cell',c,' brute <v> = %+.5f   exact = %+.5f'%(Vc[c]/N,val))
