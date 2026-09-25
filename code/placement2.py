import sys
from vort_tiling import *
from vform2 import height_norm, v_formula
def enum_cells(cells):
    cl=sorted(cells,key=lambda c:(c[1],c[0]))
    def rec(used,ds):
        free=next((c for c in cl if c not in used),None)
        if free is None: yield list(ds); return
        for dd in [(1,0),(0,1)]:
            o=(free[0]+dd[0],free[1]+dd[1])
            if o in cells and o not in used: yield from rec(used|{free,o},ds+[(free,o)])
    return rec(frozenset(),[])
for n in range(1,6):
    for label,shift in (('I',((n+1)%2,0)),('II',(n%2,0))):
        cells,_=aztec(n,shift); C=Counter(); tots=set(); fmatch=Counter()
        for doms in enum_cells(cells):
            I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
            C[('pair-conflict ' if bad else '')+('frustrated' if fr else 'consistent'), f'{comps} comp']+=1
            h=height_norm(doms,cells); f=[v_formula(d,h) for d in doms]
            fmatch[len({v[k]*f[k] for k in range(len(doms))})==1]+=1
            tots.add(abs(sum(f)))
        print(f'n={n} placement {label}: {dict(C)}  formula matches rule: {dict(fmatch)}  |Σv| values: {sorted(tots)}')
