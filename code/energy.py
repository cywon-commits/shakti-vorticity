"""Island spins of the whole sample from the tetromino templates; vertex charges; boundary excess charge."""
import sys, random
from vort_tiling import *
from placement2 import enum_cells
from vform2 import height_norm, v_formula
def island_config(doms):
    I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
    spin={}
    for k,A in enumerate(I):
        for e,(b,kind) in A['tp'].items():
            if kind=='ext': continue
            s=v[k]*b
            if e in spin: assert spin[e]==s
            spin[e]=s
    return spin, v, I
def charges(spin):
    q=Counter(); z=Counter()
    for (a,b),s in spin.items():
        # s=+1: arrow a->b : out of a, into b
        q[b]+=s; q[a]-=s; z[a]+=1; z[b]+=1
    return q,z
def full_z(p):
    x,y=p; return sum(is_island(((x,y),(x+dx,y+dy))) for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)])
def boundary_stats(doms):
    spin,v,I=island_config(doms); q,z=charges(spin)
    exc=0; hist=Counter(); nb=0
    for p in z:
        zp=z[p]; zf=full_z(p); aq=abs(q[p])
        extra=(aq-(zp%2))//2          # units of excess charge above the minimum allowed by the local coordination
        if zp<zf: nb+=1; hist[(zf,zp,aq)]+=1
        else: assert extra==0 and (zp!=4 or aq==0), (p,zp,aq)
        exc+=extra
    return exc, nb, hist, len(spin)
if __name__=='__main__':
    for n in range(1,int(sys.argv[1])+1):
        for label,shift in (('I',((n+1)%2,0)),('II',(n%2,0))):
            cells,_=aztec(n,shift); E=Counter(); H=Counter(); NB=set(); NI=set()
            it=enum_cells(cells) if n<=5 else None
            if it is None:
                rng=random.Random(n)
                it=([((a[0]+shift[0],a[1]+shift[1]),(b[0]+shift[0],b[1]+shift[1])) for a,b in shuffle(n,rng)] for _ in range(300))
            for doms in it:
                exc,nb,hist,ni=boundary_stats(doms); E[exc]+=1; NB.add(nb); NI.add(ni); H.update(hist)
            tot=sum(E.values())
            print(f'n={n:2d} {label:2s} boundary vertices {sorted(NB)} islands {sorted(NI)}  excess-charge units per tiling {dict(sorted(E.items()))}  mean {sum(k*c for k,c in E.items())/tot:.2f}')
            print('      boundary vertex types (full z, sample z, |q|): counts/tiling', {k:round(c/tot,2) for k,c in sorted(H.items())})
