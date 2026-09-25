"""Square L x L board with a central square hole (size H), height-function Monte Carlo with a loop move
that shifts the whole inner boundary by +-4 (mode 'free') or +-8 (mode 'fixed': inner vorticity pinned).
mode 'nohole' = ordinary board.  Vorticity from v(d) = s(d)(-1)^floor(h/4)."""
import sys, time, pickle, numpy as np
from vform2 import S_TAB, r4
from placement2 import height_norm
SH=np.zeros((4,4)); SV=np.zeros((4,4))
for x in range(4):
    for y in range(4): SH[x,y]=S_TAB[('H',x,y)]; SV[x,y]=S_TAB[('V',x,y)]
def build(L,Hh):
    cell=np.ones((L,L),bool); o=(L-Hh)//2
    if Hh>0: cell[o:o+Hh,o:o+Hh]=False
    return cell,o
def run(L,Hh,mode,burn,nmeas,gap,seed,target_parity=0):
    rng=np.random.default_rng(seed)
    cell,o=build(L,Hh if mode!='nohole' else 0)
    cells={(i,j) for i in range(L) for j in range(L) if cell[i,j]}
    doms=[]
    for j in range(L):
        i=0
        while i<L:
            if cell[i,j]: doms.append(((i,j),(i+1,j))); i+=2
            else: i+=1
    hd=height_norm(doms,cells)
    X,Y=np.meshgrid(np.arange(L+1),np.arange(L+1),indexing='ij')
    h=np.zeros((L+1,L+1),np.int64); vdef=np.zeros((L+1,L+1),bool)
    for (p,val) in hd.items(): h[p]=val; vdef[p]=True
    R=np.vectorize(lambda a,b:r4((a,b)))(X,Y)
    blk=lambda i,j: ((i+j)%2==0)
    sgn={'px':np.where(blk(X,Y),1,-1),'mx':np.where(blk(X-1,Y-1),1,-1),'py':np.where(blk(X-1,Y),1,-1),'my':np.where(blk(X,Y-1),1,-1)}
    C=np.zeros((L+2,L+2),bool); C[1:-1,1:-1]=cell       # padded cell mask, C[i+1,j+1]=cell(i,j)
    # vertex (i,j) touches cells (i-1,j-1),(i,j-1),(i-1,j),(i,j) -> padded indices
    t00=C[0:L+1,0:L+1]; t10=C[1:L+2,0:L+1]; t01=C[0:L+1,1:L+2]; t11=C[1:L+2,1:L+2]
    interior=t00&t10&t01&t11
    touch=t00|t10|t01|t11
    inner=np.zeros_like(interior)
    if mode!='nohole':
        hole_touch=np.zeros_like(interior); Hc=np.zeros((L+2,L+2),bool); Hc[1+o:1+o+Hh,1+o:1+o+Hh]=True
        hole_touch=Hc[0:L+1,0:L+1]|Hc[1:L+2,0:L+1]|Hc[0:L+1,1:L+2]|Hc[1:L+2,1:L+2]
        inner=hole_touch&touch
    # edges to check for the loop move: horizontal/vertical edges with at least one adjacent region cell
    ehx=np.zeros((L,L+1),bool)   # edge (i,j)-(i+1,j): cells (i,j-1),(i,j)
    ehx=C[1:L+1,0:L+1]|C[1:L+1,1:L+2]
    evy=C[0:L+1,1:L+1]|C[1:L+2,1:L+1]   # edge (i,j)-(i,j+1): cells (i-1,j),(i,j)
    def valid_all(hh):
        d=(hh[1:,:]-hh[:-1,:])*sgn['px'][:-1,:]; ok=np.all(((d==1)|(d==-3))|~ehx)
        d=(hh[:,1:]-hh[:,:-1])*sgn['py'][:,:-1]; return ok and np.all(((d==1)|(d==-3))|~evy)
    assert valid_all(h)
    ref=tuple(np.argwhere(inner)[0]) if inner.any() else None
    h0=h[ref] if ref else 0
    step=4 if mode=='free' else 8
    def sweep():
        for a in (0,1):
            for b in (0,1):
                sl=(slice(a,L+1,2),slice(b,L+1,2)); I,J=np.meshgrid(np.arange(a,L+1,2),np.arange(b,L+1,2),indexing='ij')
                m=interior[sl]; hp=h[sl]; prop=hp+4*rng.choice([-1,1],size=hp.shape); ok=m.copy()
                for (di,dj,key) in ((1,0,'px'),(-1,0,'mx'),(0,1,'py'),(0,-1,'my')):
                    ii=np.clip(I+di,0,L); jj=np.clip(J+dj,0,L)
                    d=(h[ii,jj]-prop)*sgn[key][I,J]; ok&=(d==1)|(d==-3)
                hp[ok]=prop[ok]; h[sl]=hp
        if mode in ('free','fixed'):
            dlt=step*rng.choice([-1,1]); hn=h.copy(); hn[inner]+=dlt
            if valid_all(hn): h[:]=hn
    for _ in range(burn): sweep()
    if mode=='fixed':   # reach the requested parity of the inner offset using free moves, then keep it
        tries=0
        while ((h[ref]-h0)//4)%2!=target_parity and tries<200000:
            hn=h.copy(); hn[inner]+=4*rng.choice([-1,1])
            if valid_all(hn): h[:]=hn
            else:
                for _ in range(1): 
                    step_save=step; sweep()
            tries+=1
        assert ((h[ref]-h0)//4)%2==target_parity
        for _ in range(burn//4): sweep()
    vs=np.zeros((L,L)); tots=[]; ks=[]
    for mm in range(nmeas):
        for _ in range(gap): sweep()
        dV=(h[:,1:]-h[:,:-1])*sgn['py'][:,:-1]; dH=(h[1:,:]-h[:-1,:])*sgn['px'][:-1,:]
        Hd=np.zeros((L,L),bool); Hd[:-1,:]=(dV[1:L,:L]==-3)&cell[:-1,:]&cell[1:,:]
        Vd=np.zeros((L,L),bool); Vd[:,:-1]=(dH[:L,1:L]==-3)&cell[:,:-1]&cell[:,1:]
        eps=np.where(((h-(h%4))//4)%2==0,1,-1)
        vH=SH[X[:L,:L]%4,Y[:L,:L]%4]*eps[:L,:L]; vV=SV[X[:L,:L]%4,Y[:L,:L]%4]*eps[:L,:L]
        v=np.zeros((L,L)); v[Hd]=vH[Hd]; v[1:,:][Hd[:-1,:]]=vH[:-1,:][Hd[:-1,:]]
        v[Vd]=vV[Vd]; v[:,1:][Vd[:,:-1]]=vV[:,:-1][Vd[:,:-1]]
        assert np.all((v!=0)==cell)
        vs+=v; tots.append(v.sum()/2); ks.append((h[ref]-h0)//4 if ref else 0)
    return dict(L=L,H=Hh,mode=mode,parity=target_parity,N=nmeas,vmean=vs/nmeas,cell=cell,tot=np.array(tots),k=np.array(ks))
if __name__=='__main__':
    L,Hh=int(sys.argv[1]),int(sys.argv[2]); mode=sys.argv[3]; par=int(sys.argv[4]) if len(sys.argv)>4 else 0
    burn=20*L*L; nm=int(sys.argv[5]) if len(sys.argv)>5 else 3000; gap=L*L//16
    t=time.time(); D=run(L,Hh,mode,burn,nm,gap,L*1000+Hh*10+par+hash(mode)%97,par)
    tag=f'{mode}{par if mode=="fixed" else ""}'
    pickle.dump(D,open(f'hole/h_{L}_{Hh}_{tag}.pkl','wb'))
    from collections import Counter
    print(L,Hh,tag,'time %.0fs'%(time.time()-t),'Σv values',dict(Counter(D['tot'].astype(int).tolist())),'inner offset k:',dict(Counter(D['k'].tolist())),flush=True)
