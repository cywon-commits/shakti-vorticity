import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from placement2 import *
from energy import boundary_stats
import itertools
res=defaultdict(Counter)
for L in (2,4,6):
    for sx in range(4):
        for sy in range(2):
            cells={(sx+i,sy+j) for i in range(L) for j in range(L)}
            tots=set(); cons=Counter(); fm=Counter(); nt=0
            for doms in enum_cells(cells):
                nt+=1
                if nt>3000: break
                I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
                cons[('ok' if not (bad or fr) else 'FRUSTRATED'),comps]+=1
                h=height_norm(doms,cells); f=[v_formula(d,h) for d in doms]
                fm[len({v[k]*f[k] for k in range(len(doms))})==1]+=1; tots.add(int(sum(f)))
            print(f'L={L} shift=({sx},{sy}) tilings checked={min(nt,3000)} {dict(cons)} formula-match {dict(fm)} Σv {sorted(tots)}')
