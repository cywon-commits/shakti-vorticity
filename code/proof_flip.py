import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from proof_local import *
R=Counter()
for x in range(4):
    for y in range(4):
        H=[((x,y),(x+1,y)),((x,y+1),(x+1,y+1))]; V=[((x,y),(x,y+1)),((x+1,y),(x+1,y+1))]
        hH=heights_on(H,(x,y)); hV=heights_on(V,(x,y))
        bnd=[(x,y),(x+1,y),(x+2,y),(x+2,y+1),(x+2,y+2),(x+1,y+2),(x,y+2),(x,y+1)]
        assert all(hH[p]==hV[p] for p in bnd)
        sH=sum(vform(d,hH) for d in H); sV=sum(vform(d,hV) for d in V)
        R[(sH==sV, 'ferro-type' if abs(sH)==2 else 'antiferro-type', (x+y)%2)]+=1
print(dict(R))
