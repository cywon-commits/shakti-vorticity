import io, contextlib, random
with contextlib.redirect_stdout(io.StringIO()):
    from placement2 import *
    import energy2 as E2
from vform2 import height_norm, v_formula
# Σv law from brick tiling
law={}
for L in range(2,62,2):
    for sh in ((0,0),(1,0)):
        cells={(sh[0]+i,sh[1]+j) for i in range(L) for j in range(L)}
        doms=[((sh[0]+i,sh[1]+j),(sh[0]+i+1,sh[1]+j)) for j in range(L) for i in range(0,L,2)]
        h=height_norm(doms,cells); law[(L,sh)]=int(sum(v_formula(d,h) for d in doms))
print('Sigma v (square, brick tiling):',{L:(law[(L,(0,0))],law[(L,(1,0))]) for L in range(2,62,2)})
# boundary energy, both classes
for L in (2,4,6,8,10):
    for sh in ((0,0),(1,0)):
        cells={(sh[0]+i,sh[1]+j) for i in range(L) for j in range(L)}
        R=Counter(); nisl=set()
        it=enum_cells(cells)
        for k,doms in enumerate(it):
            if k>=60: break
            b1,r1=E2.relaxed(doms,'exc'); b2,r2=E2.relaxed(doms,'q2'); R[(b1,r1,b2,r2)]+=1
            sp,_,_=E2.island_config(doms); nisl.add(len(sp))
        print(f'L={L:2d} corner cell {"black" if sum(sh)%2==0 else "white"}: islands {sorted(nisl)}  (excess template, relaxed | Σq² template, relaxed): {dict(R)}')
