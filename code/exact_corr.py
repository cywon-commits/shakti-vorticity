"""Exact phase correlation G(p,q) = < exp(i pi (h(q)-h(p))/4) > = e^{i pi dh0/4} det K_gamma / det K."""
import sys, math, cmath, pickle, time
from exact_v import *
def G(R,p,q):
    K,h0=R.seam(p,0,q)
    try: la,ph=R.logdet(K)
    except RuntimeError: return 0j
    return cmath.exp(1j*math.pi*h0/4)*(ph/R.ph0*math.exp(la-R.ld0)).real
if __name__=='__main__':
    if sys.argv[1]=='test':
        # brute force on a 6x6 board
        from vform2 import height_norm
        _p=open('placement2.py').read(); exec(_p[_p.index('def enum_cells'):_p.index('for n in range(1,6)')])
        cells={(i,j) for i in range(6) for j in range(6)}; R=Region(cells)
        tl=list(enum_cells(cells)); p=(3,3)
        for q in [(4,3),(5,3),(5,5),(3,6),(6,6),(1,2)]:
            b=sum(cmath.exp(1j*math.pi*(h[q]-h[p])/4) for h in (height_norm(d,cells) for d in tl))/len(tl)
            print(q,'brute',b,'exact',G(R,p,q))
        sys.exit()
    kind=sys.argv[1]; L=int(sys.argv[2]); Rw=int(sys.argv[3])
    if kind=='square': cells={(i,j) for i in range(L) for j in range(L)}; p0=(L//2,L//2)
    else:
        cells,sh=aztec_cells(L,kind); p0=(sh[0],0)
    R=Region(cells); t=time.time(); out={}
    for a in (0,1):
        for b in (0,1):
            p=(p0[0]+a,p0[1]+b)
            for dx in range(-Rw,Rw+1):
                for dy in range(-Rw,Rw+1):
                    q=(p[0]+dx,p[1]+dy)
                    if q in R.V: out[(p,(dx,dy))]=G(R,p,q)
            print(kind,L,'origin',p,'%.0fs'%(time.time()-t),flush=True)
            pickle.dump(dict(kind=kind,L=L,Rw=Rw,G=out,p0=p0),open(f'exact/corr_{kind}_{L}.pkl','wb'))
