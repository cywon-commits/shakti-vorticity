import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from vform2 import height_norm, v_formula
from exact_v import start_vertex
for corner,(sx,sy) in (('black',(0,0)),('white',(1,0))):
    out=[]
    for L in range(2,65,2):
        cells={(sx+i,sy+j) for i in range(L) for j in range(L)}
        doms=[((sx+i,sy+j),(sx+i+1,sy+j)) for j in range(L) for i in range(0,L,2)]
        h=height_norm(doms,cells); s,hs=start_vertex(cells); off=h[s]-hs; hh={p:x-off for p,x in h.items()}
        out.append((L,int(round(sum(v_formula(d,hh) for d in doms)))))
    print(corner,' '.join(f'{L}:{V:+d}' for L,V in out))
