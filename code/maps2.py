# Exact <v> in a window around the centre of the Aztec diamond (n=36-39) for placement I or II.
# Usage (from code/):  python maps2.py I ; python maps2.py II   -> exact/maps_I.pkl, exact/maps_II.pkl (read by pf_fig4.py)
import sys, pickle, time
from exact_v import *
pl=sys.argv[1]; res={}
for n in [36,37,38,39]:
    t=time.time(); cells,sh=aztec_cells(n,pl); R=Region(cells); s,hs=start_vertex(cells); cx,cy=sh[0],0
    win=[c for c in cells if abs(c[0]+.5-cx)<=12.5 and abs(c[1]+.5-cy)<=12.5]
    res[n]=dict(n=n,sh=sh,V=v_map(R,s,hs,win),centre=(cx,cy),cells=cells)
    print(pl,n,len(win),'cells  %.0fs'%(time.time()-t),flush=True)
pickle.dump(res,open(f'exact/maps_{pl}.pkl','wb'))
