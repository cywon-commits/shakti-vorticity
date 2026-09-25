"""Numba domino shuffling + height function + vorticity formula for large Aztec diamonds."""
import sys, time, pickle, numpy as np
from numba import njit
from vform2 import S_TAB, r4
SH=np.zeros((4,4),np.int64); SV=np.zeros((4,4),np.int64)
for x in range(4):
    for y in range(4): SH[x,y]=S_TAB[('H',x,y)]; SV[x,y]=S_TAB[('V',x,y)]
R4=np.zeros((2,2),np.int64)
for a in range(2):
    for b in range(2): R4[a,b]=r4((a,b))
@njit(cache=True)
def indiamond(i,j,m): return abs(2*i+1)+abs(2*j+1)<=2*m
@njit(cache=True)
def shuffle(n,seed):
    np.random.seed(seed)
    O=n+1; S=2*n+3
    T=np.zeros((S,S),np.int8)   # 1=N 2=S 3=E 4=W at anchor
    for k in range(n):
        # destruction
        for i in range(-k-1,k+1):
            for j in range(-k-1,k+1):
                t=T[i+O,j+O]
                if t==2 and T[i+O,j-1+O]==1: T[i+O,j+O]=0; T[i+O,j-1+O]=0
                elif t==3 and T[i+1+O,j+O]==4: T[i+O,j+O]=0; T[i+1+O,j+O]=0
        # sliding
        U=np.zeros((S,S),np.int8)
        for i in range(-k-1,k+1):
            for j in range(-k-1,k+1):
                t=T[i+O,j+O]
                if t==1: U[i+O,j+1+O]=1
                elif t==2: U[i+O,j-1+O]=2
                elif t==3: U[i+1+O,j+O]=3
                elif t==4: U[i-1+O,j+O]=4
        T=U
        m=k+1
        occ=np.zeros((S,S),np.bool_)
        for i in range(-m-1,m+1):
            for j in range(-m-1,m+1):
                t=T[i+O,j+O]
                if t==1 or t==2: occ[i+O,j+O]=True; occ[i+1+O,j+O]=True
                elif t==3 or t==4: occ[i+O,j+O]=True; occ[i+O,j+1+O]=True
        for j in range(m-1,-m-1,-1):
            for i in range(-m,m):
                if indiamond(i,j,m) and not occ[i+O,j+O]:
                    if np.random.random()<0.5:
                        T[i+O,j+O]=1; T[i+O,j-1+O]=2
                    else:
                        T[i+O,j-1+O]=4; T[i+1+O,j-1+O]=3
                    occ[i+O,j+O]=True; occ[i+1+O,j+O]=True; occ[i+O,j-1+O]=True; occ[i+1+O,j-1+O]=True
    return T
@njit(cache=True)
def partner_array(T,n):
    # part[i,j] = code of the partner cell direction: 0 right,1 left,2 up,3 down (cells indexed with offset O)
    O=n+1; S=2*n+3
    P=-np.ones((S,S),np.int8)
    for i in range(-n,n):
        for j in range(-n,n):
            t=T[i+O,j+O]
            if t==1 or t==2: P[i+O,j+O]=0; P[i+1+O,j+O]=1
            elif t==3 or t==4: P[i+O,j+O]=2; P[i+O,j+1+O]=3
    return P
@njit(cache=True)
def blk(i,j,shx,shy): return ((i+shx+j+shy)%2)==0
@njit(cache=True)
def height(T,n,shx,shy,R4):
    O=n+1; S=2*n+3
    P=partner_array(T,n)
    NONE=-10**9
    h=np.full((S+1,S+1),NONE,np.int64)   # vertex (x,y) -> h[x+O, y+O], x in [-n, n]
    # start at vertex (-1, -n) : lower corner of bottom cells (cell (-1,-n) is in the diamond)
    sx=-1; sy=-n
    ax=sx+shx; ay=sy+shy
    h[sx+O,sy+O]=R4[ax%2,ay%2]
    stack=np.zeros(((S+1)*(S+1),2),np.int64); top=0
    stack[0,0]=sx; stack[0,1]=sy; top=1
    dxs=(1,-1,0,0); dys=(0,0,1,-1)
    while top>0:
        top-=1; px=stack[top,0]; py=stack[top,1]; hp=h[px+O,py+O]
        for d in range(4):
            qx=px+dxs[d]; qy=py+dys[d]
            if qx<-n or qx>n or qy<-n or qy>n: continue
            if d==0: Lx,Ly,Rx,Ry=px,py,px,py-1
            elif d==1: Lx,Ly,Rx,Ry=qx,qy-1,qx,qy
            elif d==2: Lx,Ly,Rx,Ry=px-1,py,px,py
            else: Lx,Ly,Rx,Ry=qx,qy,qx-1,qy
            iL=indiamond(Lx,Ly,n); iR=indiamond(Rx,Ry,n)
            if not iL and not iR: continue
            if iL and iR:
                pl=P[Lx+O,Ly+O]
                # partner of L
                if pl==0: ox,oy=Lx+1,Ly
                elif pl==1: ox,oy=Lx-1,Ly
                elif pl==2: ox,oy=Lx,Ly+1
                else: ox,oy=Lx,Ly-1
                if ox==Rx and oy==Ry: continue
            val=hp+(1 if blk(Lx,Ly,shx,shy) else -1)
            if h[qx+O,qy+O]==NONE:
                h[qx+O,qy+O]=val; stack[top,0]=qx; stack[top,1]=qy; top+=1
    return h
@njit(cache=True)
def accumulate(n,N,seed0,shx,shy,SH,SV,R4):
    O=n+1; S=2*n+3
    vs=np.zeros((S,S)); es=np.zeros((S+1,S+1)); hc=np.zeros(N,np.int64); tot=np.zeros(N,np.int64)
    for t in range(N):
        T=shuffle(n,seed0+t)
        h=height(T,n,shx,shy,R4)
        s=0
        for i in range(-n,n):
            for j in range(-n,n):
                tp=T[i+O,j+O]
                if tp==0: continue
                x=i+shx; y=j+shy; hv=h[i+O,j+O]
                q=(hv-(hv%4))//4
                e=1 if q%2==0 else -1
                if tp==1 or tp==2:
                    v=SH[x%4,y%4]*e; vs[i+O,j+O]+=v; vs[i+1+O,j+O]+=v
                else:
                    v=SV[x%4,y%4]*e; vs[i+O,j+O]+=v; vs[i+O,j+1+O]+=v
                s+=v
        tot[t]=s
        for x in range(-n,n+1):
            for y in range(-n,n+1):
                hv=h[x+O,y+O]
                if hv>-10**8:
                    r=R4[(x+shx)%2,(y+shy)%2]
                    es[x+O,y+O]+= 1 if ((hv-r)//4)%2==0 else -1
        hc[t]=h[O,O]
    return vs/N, es/N, hc, tot
if __name__=='__main__':
    n=int(sys.argv[1]); N=int(sys.argv[2]); P=sys.argv[3] if len(sys.argv)>3 else 'I'; K=4
    shx=((n+1)%2) if P=='I' else (n%2); shy=0
    t=time.time(); blocks=[]
    for b in range(K):
        vm,em,hc,tot=accumulate(n,N//K,1000*n+7+b*N,shx,shy,SH,SV,R4); blocks.append((vm,em,hc,tot))
    vm=sum(b[0] for b in blocks)/K; em=sum(b[1] for b in blocks)/K; hc=np.concatenate([b[2] for b in blocks]); tot=np.concatenate([b[3] for b in blocks])
    D=dict(n=n,N=N,placement=P,shift=(shx,shy),vmean=vm,emean=em,hcenter=hc,tot=tot,vblocks=[b[0] for b in blocks],eblocks=[b[1] for b in blocks])
    pickle.dump(D,open(f'nb/a_{P}_{n}_{N}.pkl','wb'))
    print(n,P,N,'time %.1fs'%(time.time()-t),'Σv',sorted(set(tot.tolist())),'Var h(centre) %.3f'%hc.var(),'<eps(c)> %+.3f'%em[n+1,n+1],flush=True)
