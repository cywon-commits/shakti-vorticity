"""Ground state of the open placement-II box: do the extra edge plaquettes carry a closed circulation?"""
import pickle, sys
import open_ilp as O
from exact_v import aztec_cells
out={}
for n in (5,6):
    for pl in ('II','I'):
        A,_=aztec_cells(n,pl); pts=[p for c in A for p in O.CELLS[c]]
        x0=min(p[0] for p in pts); x1=max(p[0] for p in pts)+1; y0=min(p[1] for p in pts); y1=max(p[1] for p in pts)+1
        box=O.shakti_rect(x0,y0,x1-x0,y1-y0); E,labs=O.solve(x0,y0,x1-x0,y1-y0); L=O.solve.last
        isl=L['isl']; sp=L['spins']; idx={e:k for k,e in enumerate(isl)}
        def circ(c):
            """+1/-1 if every perimeter island of the plaquette points along one circulation, 0 otherwise"""
            r=O.CELLS[c]; xs=[q[0] for q in r]; ys=[q[1] for q in r]
            X0,X1,Y0,Y1=min(xs),max(xs)+1,min(ys),max(ys)+1
            loop=[(x,Y0) for x in range(X0,X1)]+[(X1,y) for y in range(Y0,Y1)]+[(x,Y1) for x in range(X1,X0,-1)]+[(X0,y) for y in range(Y1,Y0,-1)]
            sg=set()
            for a,b in zip(loop,loop[1:]+loop[:1]):
                e=(a,b) if (a,b) in idx else ((b,a) if (b,a) in idx else None)
                if e is None: continue
                k=idx[e]; along= sp[k]==1 if e==(a,b) else sp[k]==-1   # spin +1 points start->end
                sg.add(1 if along else -1)
            return sg.pop() if len(sg)==1 else 0
        ex=sorted(box-A); cv=[circ(c) for c in ex]
        print(f'n={n} {pl}: E={E:.0f}, unhappy={sum(1 for l in labs.values() if l=="z3unhappy")}, extra edge plaquettes {len(ex)}, complete circulation {sum(1 for v in cv if v)}/{len(ex)}',flush=True)
        out[(n,pl)]=dict(box=box,A=A,bbox=(x0,y0,x1,y1),isl=isl,spins=sp,labs=labs)
pickle.dump(out,open('exact/open_gs.pkl','wb'))
