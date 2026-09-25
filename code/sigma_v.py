import io, contextlib
from vform2 import *
out=[]
for n in range(1,61):
    cells,sh=aztec(n); doms=[]; used=set()
    for c in sorted(cells):
        if c in used: continue
        doms.append((c,(c[0]+1,c[1]))); used|={c,(c[0]+1,c[1])}
    assert used==cells
    h=height_norm(doms,cells)
    vs={tuple(sorted(d)):v_formula(d,h) for d in doms}
    spin,free,ring,_=analyse(doms,cells)
    ref=[vs[tuple(sorted(d))] for d in doms if short_cell(*d)==ring[0]][0]
    out.append((n,sum(vs.values())*ref))
print(out)
