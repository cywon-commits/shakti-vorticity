import sys, math, cmath, pickle, time, numpy as np
from exact_v import *
out=[]
for n in range(int(sys.argv[1]),int(sys.argv[2])+1):
    t=time.time(); cells,sh=aztec_cells(n); R=Region(cells); s,hs=start_vertex(cells)
    err=abs(R.ld0-n*(n+1)/2*math.log(2))
    cx,cy=sh[0],0      # centre vertex
    ec=R.expect_eps(s,hs,(cx,cy))
    win=[c for c in cells if math.hypot(c[0]+.5-cx,c[1]+.5-cy)<=4]
    V=v_map(R,s,hs,win)
    a=sum(V[c]*cmath.exp(-1j*math.pi/2*(c[0]+c[1])) for c in win)/len(win)
    b=sum(V[c]*cmath.exp(-1j*math.pi/2*(c[0]-c[1])) for c in win)/len(win)
    A=math.sqrt(2*(abs(a)**2+abs(b)**2))
    out.append(dict(n=n,eps_c=ec,A=A,alpha=a,beta=b,lnerr=err,V=V))
    print('n=%3d  <eps(centre)>=%+.6f  A_loc(r<=4)=%.5f  alpha=%+.4f%+.4fi beta=%+.4f%+.4fi  lndet err %.1e  (%.0fs)'%(n,ec,A,a.real,a.imag,b.real,b.imag,err,time.time()-t),flush=True)
pickle.dump(out,open(f'exact/scan_{sys.argv[1]}_{sys.argv[2]}.pkl','wb'))
