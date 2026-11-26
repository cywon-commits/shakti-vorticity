"""Exact finite-temperature statistics of straight-edged Shakti squares with all boundary excitations (edge_half.py):
freeable edge plaquettes (block-aligned edges) and half-dimer sites (all edges).  Replaces edge_exact.py, which had only the
first kind.  Observables: number of boundary defects N_d (paired edge plaquettes + occupied half-dimer sites; each shifts the
boundary height by 4), fraction of paired edge plaquettes, and the central amplitude |Psi|/Psi_sat (cells within distance 4
of the centre; sign referred to the bottom vertex of the Aztec diamond, seam-flipped Pfaffian as in edge_exact.py).
Usage: python edge_exact2.py kind n t1 t2 ...   kind = aligned | through, t = T/T_on, T_on = J/[(n+c) ln 2], c = 4 (aligned), 2 (through)"""
import sys, math, cmath, pickle, time, numpy as np
from pfapack.pfaffian import pfaffian
from edge_half import square, EdgeModel
from exact_v import Region, start_vertex, black
from vform2 import r4
kind=sys.argv[1]; n=int(sys.argv[2]); TS=[float(a) for a in sys.argv[3:]]; L=2*n+2; cT=4 if kind=='aligned' else 2
cells,free,nh=square(L,kind); EM=EdgeModel(cells,free,nh); idx=EM.idx
A=cells-free; RA=Region(A); v0,hs=start_vertex(A)
acx=np.mean([c[0]+.5 for c in A]); acy=np.mean([c[1]+.5 for c in A])
cen=[c for c in A if math.hypot(c[0]+.5-acx,c[1]+.5-acy)<=4]
verts=sorted({(c[0]+a,c[1]+b) for c in cen for a in (0,1) for b in (0,1)})
seams={}
for p in verts:
    pts=RA.path(v0,p); h0=hs; fl=[]
    for a,b in zip(pts,pts[1:]):
        Lc,Rc=RA.edge_cells(a,b); h0+= 1 if black(Lc) else -1
        if Lc in cells and Rc in cells: fl.append((idx[Lc],idx[Rc]))
    seams[p]=(h0,fl)
out={}
for tt in TS:
    t0=time.time(); x=(n+cT)*math.log(2)/tt
    M,Mi,pf,ph=EM.stats(x)
    nfree=pf.sum(); nhalf=ph.sum(); npair=len(free)-nfree; Nd=npair+nhalf
    eps={}
    for p,(h0,fl) in seams.items():
        k=len(fl); e0=1 if ((h0-r4(p))//4)%2==0 else -1
        if k==0: eps[p]=e0; continue
        cols=[]; Ci=np.zeros((2*k,2*k))
        for m,(i,j) in enumerate(fl):
            cols+=[i,j]; v=M[i,j]; Ci[2*m,2*m+1]=1/(2*v); Ci[2*m+1,2*m]=-1/(2*v)
        Sm=Ci+Mi[np.ix_(cols,cols)]; Sm=(Sm-Sm.T)/2
        eps[p]=e0*pfaffian(Sm)/pfaffian(Ci)
    psi={c:0.5*sum(cmath.exp(1j*math.pi*r4((c[0]+a,c[1]+b))/4)*eps[(c[0]+a,c[1]+b)] for a in (0,1) for b in (0,1)) for c in cen}
    amp=float(np.mean([abs(v) for v in psi.values()])/1.3066)
    out[tt]=dict(x=x,Nd=float(Nd),npair=float(npair),nhalf=float(nhalf),paired=float(npair/len(free)) if free else float('nan'),amp=amp)
    print(f'{kind} n={n} T/T_on={tt:7.3f} (J/T={x:6.2f}): N_d {Nd:9.3f} (paired {npair:.3f} of {len(free)}, half-dimer {nhalf:.3f} on {len(nh)} cells)  amp {amp:.4f}  {time.time()-t0:.0f}s',flush=True)
fn=f'exact/edge_half_{kind}_n{n}.pkl'
try: old=pickle.load(open(fn,'rb')); old['out'].update(out); out=old['out']
except FileNotFoundError: pass
pickle.dump(dict(n=n,kind=kind,E=len(free),H=len(nh),out=out),open(fn,'wb'))
