"""Sigma v evaluated on the all-horizontal (brick-row) tiling of the Aztec diamond, both placements, n = 1..64;
and the symmetry formula for the centre height <h(c)> = (h(p)+h(Rp))/2."""
import io, contextlib, math
with contextlib.redirect_stdout(io.StringIO()):
    from vform2 import height_norm, v_formula, r4
from exact_v import aztec_cells, start_vertex
res={}
for pl in ('I','II'):
    out=[]
    for n in range(1,65):
        cells,sh=aztec_cells(n,pl); doms=[]
        for j in sorted({c[1] for c in cells}):
            row=sorted(c for c in cells if c[1]==j)
            for k in range(0,len(row),2): doms.append((row[k],row[k+1]))
        h=height_norm(doms,cells)
        # normalise to the same convention as the exact code: h(start)=r4(start)
        s,hs=start_vertex(cells); off=h[s]-hs; hh={p:x-off for p,x in h.items()}
        V=int(round(sum(v_formula(d,hh) for d in doms)))
        # symmetry prediction for the centre
        cx,cy=sh[0],0; Rp=(cx-(s[1]-cy),cy+(s[0]-cx))
        out.append((n,V,(hh[s]+hh[Rp])/2))
    res[pl]=out
    print(pl,'Sigma v (brick tiling):',' '.join(f'{n}:{V:+d}' for n,V,_ in out))
    print(pl,'centre height from symmetry (h(p)+h(Rp))/2:',' '.join(f'{n}:{c:g}' for n,_,c in out[:16]),'...')
