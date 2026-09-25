"""Exact local XY order parameter  Psi(cell) = 1/2 sum_{corners p} e^{i pi r(p)/4} <eps(p)>,  <eps(p)> = Re[e^{-i pi r/4} Psi]."""
import sys, pickle, time, cmath
from exact_mono import *
W=9; DX=7
def psi_map(R,s,hs,cells_win):
    V={}
    need={(c[0]+a,c[1]+b) for c in cells_win for a in (0,1) for b in (0,1)}
    E={}
    for p in need:
        if p in R.V:
            try: E[p]=R.expect_eps(s,hs,p)
            except (KeyError,RuntimeError): pass
    P={}
    for c in cells_win:
        cs=[(c[0]+a,c[1]+b) for a in (0,1) for b in (0,1)]
        if all(p in E for p in cs): P[c]=0.5*sum(cmath.exp(1j*math.pi*r4(p)/4)*E[p] for p in cs)
    return P,E
if __name__=='__main__':
    pl=sys.argv[1]; ns=[int(x) for x in sys.argv[2:]]; out={}
    for n in ns:
        cells,sh=aztec_cells(n,pl); cx,cy=sh[0],0; s,hs=start_vertex(cells)
        win=[(cx+i,cy+k) for i in range(-W,W) for k in range(-W,W)]
        A=(cx-4,cy); B=(A[0]+DX,cy); t=time.time(); res=dict(n=n,pl=pl,centre=(cx,cy),AB=(A,B),Os=[O_at(A),O_at(B)])
        R=Region(cells); res['none']=psi_map(R,s,hs,[c for c in win if c in R.cells])
        R=DefRegion(cells,res['Os']); res['O']=psi_map(R,s,hs,[c for c in win if c in R.cells])
        R=MonoRegion(cells,A,B); res['mono']=psi_map(R,s,hs,[c for c in win if c in R.cells])
        out[n]=res; print(pl,n,'%.0fs'%(time.time()-t),flush=True)
        pickle.dump(out,open(f'exact/psi_{pl}.pkl','wb'))
