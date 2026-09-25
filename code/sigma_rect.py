import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from vform2 import height_norm, v_formula
from exact_v import start_vertex
for corner,(sx,sy) in (('black',(0,0)),('white',(1,0))):
    print(corner,'rows H, cols W (W even, all-horizontal tiling)')
    for H in range(1,17):
        row=[]
        for W in range(2,26,2):
            cells={(sx+i,sy+j) for i in range(W) for j in range(H)}
            doms=[((sx+i,sy+j),(sx+i+1,sy+j)) for j in range(H) for i in range(0,W,2)]
            h=height_norm(doms,cells); s,hs=start_vertex(cells); off=h[s]-hs; hh={p:x-off for p,x in h.items()}
            row.append(int(round(sum(v_formula(d,hh) for d in doms))))
        print(f'H={H:2d}',' '.join(f'{v:+d}' for v in row))
