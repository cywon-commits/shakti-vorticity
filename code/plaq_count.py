from shakti_torus import *
DEF=('z4q2','z3q3','z2q','z4II','z4q4')
def corners(T,r): return {T.red((x+dx,y+dy)) for (x,y) in r for dx in (0,1) for dy in (0,1)}
def plaq_account(T,spins):
    labs={p:classify(*T.vstate(p,spins)) for p in T.verts}
    res=[]
    for i,r in enumerate(T.plaqs):
        cs=corners(T,r)
        nu=sum(labs.get(p)=='z3unhappy' for p in cs)
        nc=tuple(sorted(labs[p] for p in cs if labs.get(p) in DEF))
        res.append((nu,nc))
    return labs,res
if __name__=='__main__':
    T=Torus(8); cnt=Counter()
    for p in T.verts:
        z=len(T.inc[p]); k=sum(1 for r in T.plaqs if p in corners(T,r)); cnt[(z,k)]+=1
    print('(z, #plaquettes touching)',cnt)
