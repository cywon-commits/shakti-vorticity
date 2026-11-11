"""Finite-temperature selection of the boundary in a block-aligned straight-edged Shakti sample.
Model (validated by the open vertex-model ILP, edge_partial_ilp.py): the Aztec-II cells must be tiled; each of the
2(n+1) edge plaquettes is either free (complete vortex, no unhappy vertex) or paired into a tetromino; every tetromino
costs J (one unhappy vertex). Weight of a configuration: exp(-J N_dominoes / T) -> fugacity z = exp(J/2T) per free edge
plaquette.  Monte Carlo: monomer-dimer Metropolis with auxiliary interior monomer fugacity w; measurements only in
configurations without interior monomers (exact conditioning).  Observables: fraction of paired edge plaquettes and the
central amplitude |Psi|/Psi_sat (cells within distance 4 of the centre), with the sign referred to a fixed corner vertex.
Usage: python edge_mc.py n sweeps seed x1 x2 ...   (x = J/T)"""
import sys, math, cmath, numpy as np, numba as nb, io, contextlib, pickle, time
with contextlib.redirect_stdout(io.StringIO()):
    import open_ilp as O
from exact_v import aztec_cells, black, Region, start_vertex
from vform2 import r4
n=int(sys.argv[1]); SW=int(sys.argv[2]); seed=int(sys.argv[3]); XS=[float(a) for a in sys.argv[4:]]
A,_=aztec_cells(n,'II'); pts=[p for c in A for p in O.CELLS[c]]
x0=min(p[0] for p in pts); x1=max(p[0] for p in pts)+1; y0=min(p[1] for p in pts); y1=max(p[1] for p in pts)+1
box=O.shakti_rect(x0,y0,x1-x0,y1-y0); E=box-A
cells=sorted(box); idx={c:k for k,c in enumerate(cells)}; N=len(cells)
isedge=np.array([c in E for c in cells],dtype=np.int8)
nbr=-np.ones((N,4),dtype=np.int64)
for c,k in idx.items():
    for d,(a,b) in enumerate(((1,0),(-1,0),(0,1),(0,-1))):
        q=(c[0]+a,c[1]+b)
        if q in idx: nbr[k,d]=idx[q]
# 2x2 blocks (lower-left cell)
blocks=[]
for (x,y) in cells:
    q=[(x,y),(x+1,y),(x,y+1),(x+1,y+1)]
    if all(t in idx for t in q): blocks.append([idx[t] for t in q])
blocks=np.array(blocks,dtype=np.int64)
# initial: Aztec rows are even -> pair (x,y),(x,y+1) along y? find a perfect matching of A greedily by rows
part=-np.ones(N,dtype=np.int64)
rowsA={}
for (x,y) in A: rowsA.setdefault(x,[]).append(y)
for x,ys in rowsA.items():
    ys=sorted(ys); assert len(ys)%2==0 and ys==list(range(ys[0],ys[0]+len(ys)))
    for k in range(0,len(ys),2): a,b=idx[(x,ys[k])],idx[(x,ys[k+1])]; part[a]=b; part[b]=a
# height tree on the vertices of A, rooted at the bottom vertex of A
RA=Region(A); v0,h0=start_vertex(A)
from collections import deque
V=sorted(RA.V); vid={v:k for k,v in enumerate(V)}
par=-np.ones(len(V),dtype=np.int64); inc=np.zeros(len(V),dtype=np.int64); eL=-np.ones(len(V),dtype=np.int64); eR=-np.ones(len(V),dtype=np.int64)
order=[vid[v0]]; seen={v0}; Q=deque([v0])
while Q:
    u=Q.popleft()
    for d in ((1,0),(-1,0),(0,1),(0,-1)):
        v=(u[0]+d[0],u[1]+d[1])
        if v in RA.V and v not in seen:
            L,R=RA.edge_cells(u,v)
            if L in A or R in A:
                seen.add(v); Q.append(v); k=vid[v]; par[k]=vid[u]; order.append(k)
                inc[k]=1 if black(L) else -1
                eL[k]=idx.get(L,-1); eR[k]=idx.get(R,-1)
order=np.array(order,dtype=np.int64)
# central cells (distance <= 4 from the centre of A) and their corner vertices
cx=np.mean([c[0]+.5 for c in A]); cy=np.mean([c[1]+.5 for c in A])
cen=[c for c in A if math.hypot(c[0]+.5-cx,c[1]+.5-cy)<=4]
corner=np.array([[vid[(c[0]+a,c[1]+b)] for a,b in ((0,0),(1,0),(0,1),(1,1))] for c in cen],dtype=np.int64)
phase=np.array([[cmath.exp(1j*math.pi*r4((c[0]+a,c[1]+b))/4) for a,b in ((0,0),(1,0),(0,1),(1,1))] for c in cen])
@nb.njit(cache=True)
def sweep(part,nbr,isedge,blocks,z,w,nsteps):
    N=part.shape[0]
    for _ in range(nsteps):
        if np.random.random()<0.5:
            b=blocks[np.random.randint(blocks.shape[0])]
            a0,a1,a2,a3=b[0],b[1],b[2],b[3]   # (x,y),(x+1,y),(x,y+1),(x+1,y+1)
            if part[a0]==a1 and part[a2]==a3:
                part[a0]=a2; part[a2]=a0; part[a1]=a3; part[a3]=a1
            elif part[a0]==a2 and part[a1]==a3:
                part[a0]=a1; part[a1]=a0; part[a2]=a3; part[a3]=a2
            continue
        c=np.random.randint(N); cp=nbr[c,np.random.randint(4)]
        if cp<0: continue
        fc=z if isedge[c] else w; fcp=z if isedge[cp] else w
        if part[c]==-1 and part[cp]==-1:
            if np.random.random()<1.0/(fc*fcp): part[c]=cp; part[cp]=c
        elif part[c]==cp:
            if np.random.random()<fc*fcp: part[c]=-1; part[cp]=-1
        elif part[c]==-1 and part[cp]>=0:
            d=part[cp]; fd=z if isedge[d] else w
            if np.random.random()<fd/fc: part[d]=-1; part[c]=cp; part[cp]=c
        elif part[cp]==-1 and part[c]>=0:
            d=part[c]; fd=z if isedge[d] else w
            if np.random.random()<fd/fcp: part[d]=-1; part[cp]=c; part[c]=cp
@nb.njit(cache=True)
def interior_monomers(part,isedge):
    m=0
    for i in range(part.shape[0]):
        if part[i]==-1 and isedge[i]==0: m+=1
    return m
@nb.njit(cache=True)
def eps_all(part,order,par,inc,eL,eR,h0):
    h=np.zeros(order.shape[0],dtype=np.int64); h[order[0]]=h0
    for t in range(1,order.shape[0]):
        k=order[t]; s=inc[k]
        if eL[k]>=0 and eR[k]>=0 and part[eL[k]]==eR[k]: s=s-4 if s>0 else s+4
        h[k]=h[par[k]]+s
    e=np.empty(order.shape[0])
    for k in range(order.shape[0]): e[k]=1.0 if (h[k]//4)%2==0 else -1.0
    return e
np.random.seed(seed)
@nb.njit(cache=True)
def seed_nb(s): np.random.seed(s)
seed_nb(seed)
w=float(__import__("os").environ.get("WMON","0.01")); SAT=1.3066; nE=int(isedge.sum()); res={}
for x in XS:
    z=math.exp(x/2); t=time.time()
    sweep(part,nbr,isedge,blocks,z,w,SW//5*N)      # equilibration
    ne=0; fp=[]; acc=np.zeros(len(V)); nmeas=0; tries=0
    while nmeas<SW:
        sweep(part,nbr,isedge,blocks,z,w,N); tries+=1
        if interior_monomers(part,isedge)>0: continue
        nmeas+=1
        fp.append(sum(1 for i in np.where(isedge==1)[0] if part[i]>=0)/nE)
        acc+=eps_all(part,order,par,inc,eL,eR,h0)
        if tries>50*SW: break
    epsm=acc/max(nmeas,1)
    psi=np.array([0.5*np.sum(phase[i]*epsm[corner[i]]) for i in range(len(cen))])
    fp=np.array(fp); B=10; blk=fp[:len(fp)//B*B].reshape(B,-1).mean(1) if len(fp)>=B else fp
    res[x]=dict(paired=fp.mean(),paired_err=blk.std()/math.sqrt(len(blk)),amp=np.abs(psi).mean()/SAT,psi=psi,nmeas=nmeas,tries=tries)
    print(f'n={n} x=J/T={x:5.2f} T/J={1/x if x>0 else float("inf"):.3f}  paired fraction {fp.mean():.3f}+-{res[x]["paired_err"]:.3f}  amp {np.abs(psi).mean()/SAT:.3f}  valid {nmeas}/{tries}  {time.time()-t:.0f}s',flush=True)
pickle.dump(dict(n=n,res=res),open(f'exact/edge_mc_n{n}_s{seed}.pkl','wb'))
