"""Uniform domino tilings of an L x L board (L even) by height-function Metropolis flips (vectorised),
vorticity from v(d) = s(d)(-1)^floor(h(p_d)/4).  Board cells [sx, sx+L) x [sy, sy+L)."""
import sys, time, pickle
import numpy as np
from vform2 import S_TAB, r4
S_H=np.zeros((4,4)); S_V=np.zeros((4,4))
for x in range(4):
    for y in range(4): S_H[x,y]=S_TAB[('H',x,y)]; S_V[x,y]=S_TAB[('V',x,y)]
def setup(W,H,sx,sy):
    X,Y=np.meshgrid(np.arange(sx,sx+W+1),np.arange(sy,sy+H+1),indexing='ij')   # vertex coords
    blk=lambda i,j: ((i+j)%2==0)
    # orientation sign of edge p->q with black square on the left: +1 if black on left of p->q
    # directions: +x (left square (x,y)), -x (left square (x-1,y-1)), +y (left square (x-1,y)), -y (left square (x,y-1))
    sgn={'px':np.where(blk(X,Y),1,-1),'mx':np.where(blk(X-1,Y-1),1,-1),'py':np.where(blk(X-1,Y),1,-1),'my':np.where(blk(X,Y-1),1,-1)}
    R=np.vectorize(lambda a,b:r4((a,b)))(X,Y)
    return X,Y,sgn,R
def initial_height(W,H,sx,sy,X,Y,sgn,R):
    # all-horizontal tiling: vertical edges x = sx+1, sx+3, ... are crossed
    h=np.zeros_like(X); h[0,0]=R[0,0]
    crossedV=lambda i: ((i - 0) % 2 == 1)      # index i along x: vertical edge at X index i crossed iff i odd (interior)
    for j in range(1,H+1):   # go up the left boundary (x index 0, never crossed)
        h[0,j]=h[0,j-1]+sgn['py'][0,j-1]
    for i in range(1,W+1):
        # move +x along each row: horizontal edges are never crossed in the all-horizontal tiling
        h[i,:]=h[i-1,:]+sgn['px'][i-1,:]
    # vertical edges x-index odd interior are crossed: check consistency later via validity
    return h
def valid_edges(h,sgn):
    ok=True
    d=(h[1:,:]-h[:-1,:])*sgn['px'][:-1,:]; ok&=np.all((d==1)|(d==-3))
    d=(h[:,1:]-h[:,:-1])*sgn['py'][:,:-1]; ok&=np.all((d==1)|(d==-3))
    return ok
def sweep(h,sgn,rng,W,H):
    for a in (0,1):
        for b in (0,1):
            sl=(slice(1+a,W,2),slice(1+b,H,2))   # interior vertices of class (a,b)
            hp=h[sl]; prop=hp+4*rng.choice([-1,1],size=hp.shape)
            ia=np.arange(1+a,W,2); ib=np.arange(1+b,H,2)
            I,J=np.meshgrid(ia,ib,indexing='ij')
            ok=np.ones(hp.shape,bool)
            for (di,dj,key) in ((1,0,'px'),(-1,0,'mx'),(0,1,'py'),(0,-1,'my')):
                d=(h[I+di,J+dj]-prop)*sgn[key][I,J]
                ok&=(d==1)|(d==-3)
            hp[ok]=prop[ok]; h[sl]=hp
def vfield(h,X,Y,sgn,W,H):
    """v per cell (each cell gets v of its domino)"""
    # horizontal domino with lower-left vertex (i,j) (indices): vertical edge (i+1,j)->(i+1,j+1) crossed
    dV=(h[:,1:]-h[:,:-1])*sgn['py'][:,:-1]      # edge (i,j)->(i,j+1); crossed iff -3
    dH=(h[1:,:]-h[:-1,:])*sgn['px'][:-1,:]      # edge (i,j)->(i+1,j)
    v=np.zeros((W,H))
    eps=np.where(((h-(h%4))//4)%2==0,1,-1)
    xm=X%4; ym=Y%4
    # vectorised: horizontal dominoes at cell (i,j): crossed vertical edge at vertex index (i+1, j)
    Hd=np.zeros((W,H),bool); Hd[:-1,:]=(dV[1:W,:H]==-3)
    Vd=np.zeros((W,H),bool); Vd[:,:-1]=(dH[:W,1:H]==-3)
    sH=S_H[xm[:W,:H],ym[:W,:H]]; sV=S_V[xm[:W,:H],ym[:W,:H]]
    vH=sH*eps[:W,:H]; vV=sV*eps[:W,:H]
    v[Hd]=vH[Hd]; v[1:,:][Hd[:-1,:]]=vH[:-1,:][Hd[:-1,:]]
    v[Vd]=vV[Vd]; v[:,1:][Vd[:,:-1]]=vV[:,:-1][Vd[:,:-1]]
    return v,Hd,Vd
def run(W,H,sx,sy,burn,nmeas,gap,seed):
    rng=np.random.default_rng(seed)
    X,Y,sgn,R=setup(W,H,sx,sy); h=initial_height(W,H,sx,sy,X,Y,sgn,R)
    assert valid_edges(h,sgn) and np.all(h%4==R%4)
    for _ in range(burn): sweep(h,sgn,rng,W,H)
    vs=np.zeros((W,H)); es=np.zeros((W+1,H+1)); tots=[]; hc=[]
    for m in range(nmeas):
        for _ in range(gap): sweep(h,sgn,rng,W,H)
        v,Hd,Vd=vfield(h,X,Y,sgn,W,H); assert np.all(v!=0)
        vs+=v; es+=np.where(((h-R)//4)%2==0,1,-1); tots.append(v.sum()/2); hc.append(h[W//2,H//2])
    return dict(W=W,H=H,shift=(sx,sy),N=nmeas,vmean=vs/nmeas,emean=es/nmeas,tot=np.array(tots),hcenter=np.array(hc))
if __name__=='__main__':
    W=int(sys.argv[1]); H=int(sys.argv[2]); burn=int(sys.argv[3]); nm=int(sys.argv[4]); gap=int(sys.argv[5])
    t=time.time(); D=run(W,H,0,0,burn,nm,gap,W*100+H)
    pickle.dump(D,open(f'rect/r_{W}x{H}.pkl','wb'))
    print(W,H,'time %.0fs'%(time.time()-t),'sum v',sorted(set(D['tot'].tolist())),'Var h(c) %.2f'%D['hcenter'].var())
