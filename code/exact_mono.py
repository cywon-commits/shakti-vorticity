"""Exact <v> with two fixed monomers (half vortices) and a straight branch cut (the string) between them."""
from exact_defect import *
from collections import deque
class MonoRegion(DefRegion):
    def __init__(self,cells,A,B):
        super().__init__(cells-{A,B},[])
        assert A[1]==B[1] and A[0]<B[0]
        # cut: vertical lattice edges x = A0+1..B0 at height y=A1 -> A1+1 (crossed by the segment A->B)
        self.cut={((x,A[1]),(x,A[1]+1)) for x in range(A[0]+1,B[0]+1)}
        self.cut|={(b,a) for a,b in self.cut}
        # Kasteleyn sign string between the two single-cell holes: negate the cell pairs straddling the cut
        K=self.K.tolil()
        for x in range(A[0]+1,B[0]):
            L,R=(x,A[1]-1),(x,A[1])
            if L in self.cells and R in self.cells:
                b,w=(L,R) if black(L) else (R,L); K[self.bi[b],self.wi[w]]*=-1
        self.K=K.tocsc(); self.ld0,self.ph0=self.logdet(self.K)
    def path(self,start,p):
        prev={start:None}; Q=deque([start])
        while Q:
            u=Q.popleft()
            if u==p: break
            for d in ((1,0),(-1,0),(0,1),(0,-1)):
                v=(u[0]+d[0],u[1]+d[1])
                if v in self.V and v not in prev and (u,v) not in self.cut:
                    L,R=self.edge_cells(u,v)
                    if L in self.cells or R in self.cells: prev[v]=u; Q.append(v)
        path=[p]
        while prev[path[-1]] is not None: path.append(prev[path[-1]])
        return path[::-1]
def brute_height(doms,region):
    part={}
    for a,b in doms: part[a]=b; part[b]=a
    s,hs=start_vertex(region.full); h={s:hs}; Q=[s]; cut=getattr(region,'cut',set())
    for p in Q:
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            q=(p[0]+d[0],p[1]+d[1])
            if q not in region.V or (p,q) in cut: continue
            L,R=region.edge_cells(p,q)
            if L not in region.cells and R not in region.cells: continue
            if L in part and part[L]==R: continue
            if q not in h: h[q]=h[p]+(1 if black(L) else -1); Q.append(q)
    return h
if __name__=='__main__':
    from vform2 import S_TAB
    _p=open('placement2.py').read(); exec(_p[_p.index('def enum_cells'):_p.index('for n in range(1,6)')])
    from collections import defaultdict
    for n,pl,dx in [(5,'I',3),(5,'II',3),(5,'I',5),(6,'II',5)]:
        cells,sh=aztec_cells(n,pl); A=(sh[0]-2,0); B=(A[0]+dx,0)
        R=MonoRegion(cells,A,B)
        tl=list(enum_cells(R.cells))
        acc=defaultdict(float)
        for doms in tl:
            h=brute_height(doms,R)
            for d in doms:
                d=tuple(sorted(d)); p=d[0]; typ='H' if d[1][0]!=d[0][0] else 'V'
                if p not in h:     # p only reachable across the cut: use a perimeter vertex
                    cyc,inc=R.perim(d); i=next(i for i,q in enumerate(cyc) if q in h); hp=h[cyc[i]]-sum(inc[:i])
                else: hp=h[p]
                v=S_TAB[(typ,p[0]%4,p[1]%4)]*(1 if ((hp-r4(p))//4)%2==0 else -1)
                for c in d: acc[c]+=v
        N=len(tl); s,hs=start_vertex(cells)
        if N==0: print(n,pl,'no tilings'); continue
        Vx=v_map(R,s,hs,list(R.cells))
        print(f'n={n} {pl}: monomers {A},{B}; tilings {N} (|det K| = {math.exp(R.ld0):.1f}); max |brute-exact| = {max(abs(acc[c]/N-Vx[c]) for c in R.cells):.2e}',flush=True)
