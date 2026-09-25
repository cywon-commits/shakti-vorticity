import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from placement2 import *
for n in (3,4,5):
    for label,shift in (('I',((n+1)%2,0)),('II',(n%2,0))):
        cells,_=aztec(n,shift)
        bcells=[c for c in cells if any((c[0]+dx,c[1]+dy) not in cells for dx,dy in [(1,0),(-1,0),(0,1),(0,-1)]) or partner(c) not in cells]
        vals=defaultdict(set); P=Counter()
        for doms in enum_cells(cells):
            h=height_norm(doms,cells)
            for d in doms:
                vv=v_formula(d,h)
                for c in d:
                    if c in bcells: vals[c].add(vv)
            P[tuple(sorted(short_cell(*d) for d in doms if any(c in bcells for c in d) and short_cell(*d) in bcells))]+=1
        frozen=sum(len(s)==1 for s in vals.values())
        print(f'n={n} {label}: boundary plaquettes {len(bcells)}, with the same v in every tiling: {frozen};  distinct sets of boundary strong arms: {len(P)}')
