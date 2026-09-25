import sys, pickle, time
from exact_mono import *
W=8; DX=7
pl=sys.argv[1]; ns=[int(x) for x in sys.argv[2:]]
out={}
for n in ns:
    cells,sh=aztec_cells(n,pl); cx,cy=sh[0],0; s,hs=start_vertex(cells)
    win=[(cx+i,cy+k) for i in range(-W,W) for k in range(-W,W)]
    A=(cx-4,cy); B=(A[0]+DX,cy)
    res={}
    t=time.time()
    R=Region(cells); res['none']=v_map(R,s,hs,[c for c in win if c in R.cells])
    Os=[O_at(A),O_at(B)]; assert not set(Os[0])&set(Os[1])
    R=DefRegion(cells,Os); res['O']=v_map(R,s,hs,[c for c in win if c in R.cells]); res['Os']=Os
    R=MonoRegion(cells,A,B); res['mono']=v_map(R,s,hs,[c for c in win if c in R.cells]); res['AB']=(A,B)
    res.update(n=n,pl=pl,centre=(cx,cy))
    out[n]=res; print(pl,n,'%.0fs'%(time.time()-t),flush=True)
    pickle.dump(out,open(f'exact/defect_maps_{pl}.pkl','wb'))
