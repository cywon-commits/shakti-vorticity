"""Exact mean height <h(p)> (Thurston, normalised h(start)=r4(start)) from edge probabilities K(b,w)K^{-1}(w,b)."""
import sys, math, cmath, pickle, numpy as np
from collections import deque
from exact_v import *
def mean_height(R):
    Kd=R.K.toarray(); Ki=np.linalg.inv(Kd)          # W x B
    s,hs=start_vertex(R.cells); H={s:float(hs)}; Q=deque([s])
    while Q:
        u=Q.popleft()
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            v=(u[0]+d[0],u[1]+d[1])
            if v in H or v not in R.V: continue
            L,Rc=R.edge_cells(u,v)
            if L not in R.cells and Rc not in R.cells: continue
            sg=1 if black(L) else -1
            P=0.0
            if L in R.cells and Rc in R.cells:
                b,w=(L,Rc) if black(L) else (Rc,L); P=(Kd[R.bi[b],R.wi[w]]*Ki[R.wi[w],R.bi[b]]).real
            H[v]=H[u]+sg*(1-4*P); Q.append(v)
    return H
if __name__=='__main__':
    for pl in ('I','II'):
        for n in range(8,41):
            cells,sh=aztec_cells(n,pl); R=Region(cells); H=mean_height(R)
            cx,cy=sh[0],0
            hc=np.mean([H[(cx+a,cy+b)] for a in (0,1) for b in (0,1)])
            print(pl,n,n%8,'<h> at centre vertex %.6f   cell-centre mean %.4f'%(H[(cx,cy)],hc),flush=True)
