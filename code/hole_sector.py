"""Fixed-sector Monte Carlo on an L x L board with a central H x H hole (local flips only; the sector = winding k is conserved).
Initial tilings in a chosen sector are built from a maximum matching plus alternating cycles that wind around the hole."""
import sys, time, pickle, numpy as np, networkx as nx
from collections import deque, Counter
from vform2 import S_TAB, r4, height_norm
SH=np.zeros((4,4)); SV=np.zeros((4,4))
for x in range(4):
    for y in range(4): SH[x,y]=S_TAB[('H',x,y)]; SV[x,y]=S_TAB[('V',x,y)]
def region(L,H):
    o=(L-H)//2; cell=np.ones((L,L),bool)
    if H>0: cell[o:o+H,o:o+H]=False
    return cell,o
def crossing(L,H,o,a,b):
    """+1 if the step a->b crosses the ray (vertical bond between rows jc-1, jc at x>=o+H) upwards, -1 downwards, else 0"""
    jc=L//2
    if a[0]==b[0] and a[0]>=o+H and {a[1],b[1]}=={jc-1,jc}: return 1 if b[1]==jc else -1
    return 0
def sector(L,H,o,M):
    k=0
    for a,b in M.items():
        if (a[0]+a[1])%2==0:   # count each domino once, oriented black->white
            k+=crossing(L,H,o,a,b)
    return k
def initial_matching(cell):
    L=cell.shape[0]; G=nx.Graph()
    for i in range(L):
        for j in range(L):
            if not cell[i,j]: continue
            for di,dj in ((1,0),(0,1)):
                q=(i+di,j+dj)
                if q[0]<L and q[1]<L and cell[q]: G.add_edge((i,j),q)
    B={c for c in G if (c[0]+c[1])%2==0}
    m=nx.bipartite.hopcroft_karp_matching(G,top_nodes=B)
    assert len(m)==cell.sum(); return m,G
def shift_sector(L,H,o,M,G,direction):
    """find an alternating cycle with winding = direction and flip it"""
    blacks=[c for c in G if (c[0]+c[1])%2==0]
    for b0 in sorted(blacks,key=lambda c:abs(c[0]-(o+H))+abs(c[1]-L//2)):
        start=(b0,0); prev={start:None}; Q=deque([start]); found=None
        while Q and found is None:
            b,w=Q.popleft()
            for wn in G[b]:
                if M[b]==wn: continue
                w1=w+crossing(L,H,o,b,wn)
                bn=M[wn]; w2=w1+crossing(L,H,o,wn,bn)
                if abs(w2)>2: continue
                st=(bn,w2)
                if st in prev: continue
                prev[st]=(b,w,wn)
                if bn==b0 and w2==direction: found=st; break
                Q.append(st)
        if found:
            # reconstruct cycle and flip
            path=[]; st=found
            while prev[st] is not None:
                b,w,wn=prev[st]; path.append((b,wn,st[0])); st=(b,w)
            verts=[b for b,_,_ in path]+[wn for _,wn,_ in path]
            if len(set(verts))!=len(verts): continue      # not a simple cycle
            for b,wn,bn in path:
                M[b]=wn; M[wn]=b
            assert all(M[M[c]]==c for c in M)
            return True
    return False
def run(L,H,ktarget,burn,nmeas,gap,seed):
    rng=np.random.default_rng(seed); cell,o=region(L,H)
    M,G=initial_matching(cell); k=sector(L,H,o,M)
    while k!=ktarget:
        ok=shift_sector(L,H,o,M,G,1 if ktarget>k else -1); assert ok,'no winding cycle'
        k2=sector(L,H,o,M); assert abs(k2-k)==1; k=k2
    cells={(i,j) for i in range(L) for j in range(L) if cell[i,j]}
    doms=list({tuple(sorted((a,b))) for a,b in M.items()})
    hd=height_norm(doms,cells)
    X,Y=np.meshgrid(np.arange(L+1),np.arange(L+1),indexing='ij'); h=np.zeros((L+1,L+1),np.int64)
    for p,val in hd.items(): h[p]=val
    blk=lambda i,j: ((i+j)%2==0)
    sgn={'px':np.where(blk(X,Y),1,-1),'mx':np.where(blk(X-1,Y-1),1,-1),'py':np.where(blk(X-1,Y),1,-1),'my':np.where(blk(X,Y-1),1,-1)}
    C=np.zeros((L+2,L+2),bool); C[1:-1,1:-1]=cell
    interior=C[0:L+1,0:L+1]&C[1:L+2,0:L+1]&C[0:L+1,1:L+2]&C[1:L+2,1:L+2]
    def sweep():
        for a in (0,1):
            for b in (0,1):
                sl=(slice(a,L+1,2),slice(b,L+1,2)); I,J=np.meshgrid(np.arange(a,L+1,2),np.arange(b,L+1,2),indexing='ij')
                hp=h[sl]; prop=hp+4*rng.choice([-1,1],size=hp.shape); ok=interior[sl].copy()
                for (di,dj,key) in ((1,0,'px'),(-1,0,'mx'),(0,1,'py'),(0,-1,'my')):
                    d=(h[np.clip(I+di,0,L),np.clip(J+dj,0,L)]-prop)*sgn[key][I,J]; ok&=(d==1)|(d==-3)
                hp[ok]=prop[ok]; h[sl]=hp
    for _ in range(burn): sweep()
    vs=np.zeros((L,L)); tots=[]
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
        vs+=v; tots.append(int(v.sum()//2))
    return dict(L=L,H=H,k=ktarget,N=nmeas,vmean=vs/nmeas,cell=cell,tot=np.array(tots))
if __name__=='__main__':
    L,H,k=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]); nm=int(sys.argv[4]) if len(sys.argv)>4 else 3000
    t=time.time(); D=run(L,H,k,20*L*L,nm,L*L//16,L*100+H*10+k+5)
    pickle.dump(D,open(f'hole/s_{L}_{H}_k{k}.pkl','wb'))
    print(L,H,'k=',k,'time %.0fs'%(time.time()-t),'Σv',dict(Counter(D['tot'].tolist())),flush=True)
