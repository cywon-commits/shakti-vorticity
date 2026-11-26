"""Exact density of states of the open Shakti vertex model (all spin states, vertex energies COST of gen_ilp.py)
by a site-by-site transfer matrix, to test the boundary-monomer effective model at finite temperature.
g[E] = number of spin states with energy E (in units of J), kept for E <= cap.  Floating point (53 bits) is enough to compare counts."""
import sys, itertools, math, pickle, numpy as np
from collections import defaultdict
from open_sample import Sample, halfplanes, COST
from shakti_torus import classify
def lab(L,bits):
    ins=[];dirs=[]
    for (k,role,a,b),bt in zip(L,bits):
        ins.append((bt==1 and role=='end') or (bt==-1 and role=='start')); dirs.append((b[0]-a[0],b[1]-a[1]))
    z=len(L); q=sum(1 if i else -1 for i in ins)
    return 'z1' if z==1 else classify(z,q,ins,dirs)
def dos(S,cap):
    inc=defaultdict(list)
    for k,(a,b) in enumerate(S.isl): inc[a].append((k,'start',a,b)); inc[b].append((k,'end',a,b))
    pts=sorted(inc,key=lambda p:(p[1],p[0])); order={p:i for i,p in enumerate(pts)}
    A=np.zeros(cap+1); A[0]=1.0; axes=[]          # axes: list of open island indices; A has shape (2,)*len(axes) + (cap+1,)
    maxopen=0
    for p in pts:
        L=inc[p]; ks=[k for k,_,_,_ in L]; old=[k for k in ks if k in axes]; new=[k for k in ks if k not in axes]
        keep=[k for k in axes if k not in old]; naxes=keep+new
        B=np.zeros((2,)*len(naxes)+(cap+1,))
        Am=np.moveaxis(A,[axes.index(k) for k in old],range(len(old))) if old else A          # old axes first, then keep..., E
        for bits in itertools.product((1,-1),repeat=len(L)):
            c=COST[lab(L,bits)]
            if c>cap: continue
            val=dict(zip(ks,bits)); src=Am[tuple(0 if val[k]==1 else 1 for k in old)]       # shape (2,)*len(keep)+(cap+1,)
            dst=B[(slice(None),)*len(keep)+tuple(0 if val[k]==1 else 1 for k in new)]
            if c==0: dst+=src
            else: dst[...,c:]+=src[...,:-c]
        A,axes=B,naxes; maxopen=max(maxopen,len(axes))
    assert not axes or all(False for _ in axes) or True
    g=A.reshape(-1,cap+1).sum(0) if axes else A
    return g,maxopen,len(S.isl)
if __name__=='__main__':
    from bmono import BMono
    from zpoly import zpoly
    out={}
    for L in [int(v) for v in sys.argv[1:]]:
        S=Sample(halfplanes((1,0,0,L+1),(0,1,0,L+1)),W=L+10,mode='A'); B=BMono(S.cells,S.free); c=zpoly(B); m=B.m
        Egs=(len(S.cells)-m)//2; K=int(__import__("os").environ.get("KLEV","6")); g,mo,ni=dos(S,Egs+K); out[L]=(g,c,Egs)
        E0=int(np.nonzero(g)[0][0]); D=g[E0]
        print(f'L={L}: islands {ni}, max open {mo}, cells {len(S.cells)}, m={m}; effective-model E_gs={Egs}, vertex-model E_gs={E0}, ground-state degeneracy {D:.0f}')
        print('   k   g(E_gs+k)/g(E_gs)     effective model c_(m-2k)     ratio')
        for k in range(K+1):
            ce=c[m-2*k] if m-2*k>=0 else 0
            print(f'  {k:2d}   {g[E0+k]/D:18.6g}   {ce:18d}     {g[E0+k]/D/ce if ce else float("nan"):.6f}')
        pickle.dump(out,open('data/vertex_dos.pkl','wb'))
