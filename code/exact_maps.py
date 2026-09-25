import sys, math, pickle, time
from exact_v import *
res={}
for n,mode in [(36,'full'),(37,'full'),(38,'full'),(39,'full'),(24,'full')]:
    t=time.time(); cells,sh=aztec_cells(n); R=Region(cells); s,hs=start_vertex(cells); cx,cy=sh[0],0
    win=[c for c in cells if (mode=='full') or (abs(c[0]+.5-cx)<=8 and abs(c[1]+.5-cy)<=8)]
    res[n]=dict(n=n,sh=sh,V=v_map(R,s,hs,win),centre=(cx,cy),cells=cells)
    print(n,mode,len(win),'cells  %.0fs'%(time.time()-t),flush=True)
pickle.dump(res,open('exact/maps.pkl','wb'))
