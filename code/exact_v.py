"""Exact <eps(p)> and <v(cell)> for uniform domino tilings via Kasteleyn with a sign-flipped seam.
eps(p) = eps0(p) * (-1)^{N_gamma},  N_gamma = # dominoes crossing a lattice path gamma from a boundary vertex to p.
<(-1)^N> = det K_gamma / det K ;  <1_d (-1)^N> = <(-1)^N> * K_gamma(b,w) K_gamma^{-1}(w,b)."""
import sys, math, cmath, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from vform2 import S_TAB, r4
def aztec_cells(n,placement='I'):
    sh=((n+1)%2,0) if placement=='I' else (n%2,0)
    return {(i+sh[0],j+sh[1]) for i in range(-n,n) for j in range(-n,n) if abs(i+.5)+abs(j+.5)<=n}, sh
def black(c): return (c[0]+c[1])%2==0
class Region:
    def __init__(self,cells):
        self.cells=cells
        self.B=sorted(c for c in cells if black(c)); self.W=sorted(c for c in cells if not black(c))
        self.bi={c:k for k,c in enumerate(self.B)}; self.wi={c:k for k,c in enumerate(self.W)}
        rows=[];cols=[];vals=[]
        for b in self.B:
            for di,dj in ((1,0),(-1,0),(0,1),(0,-1)):
                w=(b[0]+di,b[1]+dj)
                if w in self.wi: rows.append(self.bi[b]); cols.append(self.wi[w]); vals.append(1.0 if dj==0 else 1j)
        self.K=sp.csc_matrix((vals,(rows,cols)),shape=(len(self.B),len(self.W)),dtype=complex)
        V=set()
        for (x,y) in cells: V|={(x,y),(x+1,y),(x,y+1),(x+1,y+1)}
        self.V=V
        self.ld0,self.ph0=self.logdet(self.K)
    @staticmethod
    def perm_sign(p):
        p=np.asarray(p); seen=np.zeros(len(p),bool); s=1
        for i in range(len(p)):
            if not seen[i]:
                j=i; L=0
                while not seen[j]: seen[j]=True; j=p[j]; L+=1
                if L%2==0: s=-s
        return s
    def logdet(self,A):
        lu=spl.splu(A.tocsc(),permc_spec='COLAMD',diag_pivot_thresh=1.0)
        d=lu.U.diagonal(); la=np.sum(np.log(np.abs(d))); ph=np.prod(d/np.abs(d))*self.perm_sign(lu.perm_r)*self.perm_sign(lu.perm_c)
        return la,ph
    def edge_cells(self,p,q):
        """the two cells on either side of lattice edge p-q, and the cell on the left of the step p->q"""
        dx,dy=q[0]-p[0],q[1]-p[1]
        if (dx,dy)==(1,0): L,R=(p[0],p[1]),(p[0],p[1]-1)
        elif (dx,dy)==(-1,0): L,R=(q[0],q[1]-1),(q[0],q[1])
        elif (dx,dy)==(0,1): L,R=(p[0]-1,p[1]),(p[0],p[1])
        else: L,R=(q[0],q[1]),(q[0]-1,q[1])
        return L,R
    def path(self,start,p):
        """BFS lattice path through region vertices along edges that touch at least one cell"""
        from collections import deque
        prev={start:None}; Q=deque([start])
        while Q:
            u=Q.popleft()
            if u==p: break
            for d in ((1,0),(-1,0),(0,1),(0,-1)):
                v=(u[0]+d[0],u[1]+d[1])
                if v in self.V and v not in prev:
                    L,R=self.edge_cells(u,v)
                    if L in self.cells or R in self.cells: prev[v]=u; Q.append(v)
        path=[p]
        while prev[path[-1]] is not None: path.append(prev[path[-1]])
        return path[::-1]
    def seam(self,start,h_start,p):
        """sign-flipped Kasteleyn matrix for the path start->p and the dimer-free height h0(p)"""
        pts=self.path(start,p); h0=h_start; flips=set()
        for a,b in zip(pts,pts[1:]):
            L,R=self.edge_cells(a,b)
            h0+= 1 if black(L) else -1
            if L in self.cells and R in self.cells:
                bc,wc=(L,R) if black(L) else (R,L); flips.add((self.bi[bc],self.wi[wc]))
        K=self.K.tolil(copy=True)
        for (i,j) in flips: K[i,j]=-K[i,j]
        return K.tocsc(),h0
    def eps0(self,h0,p): return 1 if ((h0-r4(p))//4)%2==0 else -1
    def expect_eps(self,start,h_start,p):
        K,h0=self.seam(start,h_start,p); la,ph=self.logdet(K)
        return (ph/self.ph0*math.exp(la-self.ld0)).real*self.eps0(h0,p)
    def expect_v_contrib(self,start,h_start,p,dominoes):
        """<1_d eps(p_d)> for dominoes d whose lower-left vertex is p"""
        K,h0=self.seam(start,h_start,p); out={}
        for d in dominoes:
            b,w=(d[0],d[1]) if black(d[0]) else (d[1],d[0])
            i,j=self.bi[b],self.wi[w]
            rows=np.r_[0:i,i+1:K.shape[0]]; cols=np.r_[0:j,j+1:K.shape[1]]
            M=K[rows,:][:,cols]
            try: la,ph=self.logdet(M)
            except RuntimeError: out[d]=0.0; continue     # cofactor: (-1)^{i+j} det M = det(K) K^{-1}(w,b)
            val=K[i,j]*((-1)**(i+j))*ph/self.ph0*math.exp(la-self.ld0)
            out[d]=val.real*self.eps0(h0,p)
        return out
def start_vertex(cells):
    V=set()
    for (x,y) in cells: V|={(x,y),(x+1,y),(x,y+1),(x+1,y+1)}
    s=min(V,key=lambda p:(p[1],p[0]))
    return s, r4(s)
def v_map(R,start,hs,window_cells):
    """exact <v> on the given cells"""
    need={}
    for c in window_cells:
        for d in (((c[0]-1,c[1]),c),(c,(c[0]+1,c[1])),((c[0],c[1]-1),c),(c,(c[0],c[1]+1))):
            if d[0] in R.cells and d[1] in R.cells: need.setdefault(d[0],set()).add(d)
    contrib={}
    for p,ds in need.items(): contrib.update(R.expect_v_contrib(start,hs,p,ds))
    V={}
    for c in window_cells:
        s=0
        for d in (((c[0]-1,c[1]),c),(c,(c[0]+1,c[1])),((c[0],c[1]-1),c),(c,(c[0],c[1]+1))):
            if d in contrib:
                typ='H' if d[1][0]!=d[0][0] else 'V'; s+=S_TAB[(typ,d[0][0]%4,d[0][1]%4)]*contrib[d]
        V[c]=s
    return V
