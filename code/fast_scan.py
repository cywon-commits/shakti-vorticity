"""Fast sampling using the exact formula v(d) = s(d) (-1)^floor(h(p_d)/4)."""
import sys, time, pickle, random, math
import numpy as np
from core import aztec, shuffle, black, short_cell, analyse
from vform2 import S_TAB, r4
def height_fast(doms, cells, n):
    # vertices in a (2n+3)^2 array, BFS along edges not crossing a domino
    part={}
    for a,b in doms: part[a]=b; part[b]=a
    xs=[c[0] for c in cells]; ys=[c[1] for c in cells]; x0,y0=min(xs),min(ys)
    W=max(xs)-x0+2; Hh=max(ys)-y0+2
    h=np.full((W,Hh),np.iinfo(np.int32).min,dtype=np.int64)
    start=None
    for c in cells:
        if start is None or c<start: start=c
    p=(start[0],start[1]); h[p[0]-x0,p[1]-y0]=r4(p); st=[p]
    while st:
        p=st.pop(); hp=h[p[0]-x0,p[1]-y0]
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            q=(p[0]+d[0],p[1]+d[1])
            if d==(1,0): L,R=(p[0],p[1]),(p[0],p[1]-1)
            elif d==(-1,0): L,R=(q[0],q[1]-1),(q[0],q[1])
            elif d==(0,1): L,R=(p[0]-1,p[1]),(p[0],p[1])
            else: L,R=(q[0],q[1]),(q[0]-1,q[1])
            iL=L in part; iR=R in part
            if not iL and not iR: continue
            if iL and iR and part[L]==R: continue
            val=hp+(1 if (L[0]+L[1])%2==0 else -1)
            if h[q[0]-x0,q[1]-y0]==np.iinfo(np.int32).min: h[q[0]-x0,q[1]-y0]=val; st.append(q)
    return h,x0,y0
def run(n,N,seed,placement='I'):
    cells,sh=aztec(n, None if placement=='I' else (n%2,0)); rng=random.Random(seed)
    xs=[c[0] for c in cells]; ys=[c[1] for c in cells]; x0,y0=min(xs),min(ys)
    W=max(xs)-x0+2; Hh=max(ys)-y0+2
    vsum=np.zeros((W,Hh)); esum=np.zeros((W,Hh)); hs=[]; ref=None; tot=[]
    cx=(min(xs)+max(xs)+1)//2; cy=(min(ys)+max(ys)+1)//2   # central vertex
    for t in range(N):
        doms=[((a[0]+sh[0],a[1]+sh[1]),(b[0]+sh[0],b[1]+sh[1])) for a,b in shuffle(n,rng)]
        h,_,_=height_fast(doms,cells,n)
        if ref is None and placement=='I':
            spin,free,ring,_=analyse(doms,cells); refcell=ring[0]
        elif placement!='I': refcell=None
        S=0; sgn_ref=None
        for d in doms:
            d=tuple(sorted(d)); p=d[0]; typ='H' if d[1][0]!=d[0][0] else 'V'
            hv=h[p[0]-x0,p[1]-y0]; v=S_TAB[(typ,p[0]%4,p[1]%4)]*(1 if (hv//4)%2==0 else -1)
            if refcell in d: sgn_ref=v
            for c in d: vsum[c[0]-x0,c[1]-y0]+=v   # accumulate raw, fix sign after
            S+=v
        if placement!='I': ref=1
        else:
            if ref is None: ref=sgn_ref
            assert sgn_ref==ref, 'boundary value changed'
        valid=h>np.iinfo(np.int32).min
        R=np.zeros_like(h)
        for i in range(W):
            for j in range(Hh): R[i,j]=r4((i+x0,j+y0))
        eps=np.where(valid, np.where(((h-R)//4)%2==0,1,-1),0)
        esum+=eps
        hs.append(h[cx-x0,cy-y0]); tot.append(S)
    return dict(n=n,N=N,cells=cells,x0=x0,y0=y0,vmean=vsum*ref/N,emean=esum*ref/N,hcenter=np.array(hs),tot=np.array(tot)*ref,center=(cx,cy))
if __name__=='__main__':
    n=int(sys.argv[1]); N=int(sys.argv[2]); P=sys.argv[3] if len(sys.argv)>3 else 'I'; t=time.time()
    D=run(n,N,1000+n,P); D['placement']=P
    pickle.dump(D,open(f'fast/f_{n}_{N}.pkl' if P=='I' else f'fast2/f_{n}_{N}.pkl','wb'))
    print(n,N,'time %.0fs'%(time.time()-t),'sum v',set(D['tot'].tolist()),'Var h(center) %.3f'%D['hcenter'].var(),'<eps(center)> %+.3f'%D['emean'][D['center'][0]-D['x0'],D['center'][1]-D['y0']])
