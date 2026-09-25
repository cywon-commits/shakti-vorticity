"""Exact <v> with fixed O tiles (full vortices): Kasteleyn on the region with the O cells removed."""
import io, contextlib
from exact_v import *
with contextlib.redirect_stdout(io.StringIO()):
    from core import partner
class DefRegion(Region):
    def __init__(self,cells,Os):
        self.Os=Os
        hole=set(c for O in Os for c in O)
        super().__init__(cells-hole)
        self.full=cells
        for c1,c2 in Os:     # corner shared by the two diagonal cells is inside the defect
            self.V.discard((max(c1[0],c2[0]),max(c1[1],c2[1])))
    def perim(self,d):
        p=d[0]
        if d[1][0]!=d[0][0]: cyc=[p,(p[0]+1,p[1]),(p[0]+2,p[1]),(p[0]+2,p[1]+1),(p[0]+1,p[1]+1),(p[0],p[1]+1)]
        else: cyc=[p,(p[0]+1,p[1]),(p[0]+1,p[1]+1),(p[0]+1,p[1]+2),(p[0],p[1]+2),(p[0],p[1]+1)]
        inc=[]
        for a,b in zip(cyc,cyc[1:]+cyc[:1]):
            L,R=self.edge_cells(a,b); inc.append(1 if black(L) else -1)
        return cyc,inc
    def expect_v_contrib(self,start,h_start,p,dominoes):
        if p in self.V: return super().expect_v_contrib(start,h_start,p,dominoes)
        out={}
        for d in dominoes:
            cyc,inc=self.perim(d); i=next(i for i,q in enumerate(cyc) if q in self.V); q=cyc[i]
            K,h0=self.seam(start,h_start,q); h0p=h0-sum(inc[:i])
            b,w=(d[0],d[1]) if black(d[0]) else (d[1],d[0]); ii,jj=self.bi[b],self.wi[w]
            rows=np.r_[0:ii,ii+1:K.shape[0]]; cols=np.r_[0:jj,jj+1:K.shape[1]]
            try: la,ph=self.logdet(K[rows,:][:,cols])
            except RuntimeError: out[d]=0.0; continue
            val=K[ii,jj]*((-1)**(ii+jj))*ph/self.ph0*math.exp(la-self.ld0)
            out[d]=val.real*self.eps0(h0p,p)
        return out
def O_at(c): return tuple(sorted((c,partner(c))))
if __name__=='__main__':
    import sys
    _src=open('otet_formula.py').read().split('res=Counter()')[0]; exec(_src)
    from otet import enum_with_defects
    from vform2 import S_TAB
    from collections import defaultdict
    for n,pl in ([(5,'II'),(6,'II')] if len(sys.argv)>1 else [(4,'I'),(5,'I'),(4,'II'),(5,'II')]):
        cells,sh=aztec_cells(n,pl)
        c1=(sh[0],0) if black((sh[0],0)) else (sh[0]-1,0)
        O1=O_at(c1)
        c2=next(c for c in ([(c1[0]+3,c1[1])] if len(sys.argv)>1 else [])+[(c1[0]+1,c1[1]),(c1[0]-1,c1[1]),(c1[0],c1[1]+1),(c1[0],c1[1]-1),(c1[0]+2,c1[1]+1),(c1[0]-2,c1[1]-1)] if c in cells and O_at(c)[0] in cells and O_at(c)[1] in cells and not set(O_at(c))&set(O1))
        Os=[O1,O_at(c2)]
        # brute force over tilings containing exactly these O tiles
        rest=cells-set(Os[0])-set(Os[1]); exec(open('placement2.py').read().split('for n in range(1,6)')[0].split('def enum_cells')[1].join(['def enum_cells','']) if False else '')
        from mono_tile import random_tiling
        _p=open('placement2.py').read(); _i=_p.index('def enum_cells'); _j=_p.index('for n in range(1,6)'); exec(_p[_i:_j])
        with contextlib.redirect_stdout(io.StringIO()):
            tl=list(enum_cells(rest))
        acc=defaultdict(float); N=0; badh=0
        for doms in tl:
            tiles=[('D',)+tuple(sorted(d)) for d in doms]+[('O',)+O for O in Os]
            h,bad=height_mod8(tiles,cells); badh+=bad
            s0,hs0=start_vertex(cells); off=h[s0]-hs0
            for d in doms:
                d=tuple(sorted(d)); p=d[0]; typ='H' if d[1][0]!=d[0][0] else 'V'
                v=S_TAB[(typ,p[0]%4,p[1]%4)]*eps(corner_height(d,h)-off,p)
                for c in d: acc[c]+=v
            N+=1
        if N==0: print(n,pl,'no tilings'); continue
        R=DefRegion(cells,Os); s,hs=start_vertex(cells)
        win=[c for c in rest]
        Vx=v_map(R,s,hs,win)
        err=max(abs(acc[c]/N-Vx[c]) for c in win)
        print(f'n={n} {pl}: O tiles {Os}, tilings {N}, height inconsistencies {badh}, max |brute - exact| over {len(win)} cells = {err:.2e}',flush=True)
