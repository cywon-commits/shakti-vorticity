import sys
from vort_tiling import *
from stats_run import enum
n=int(sys.argv[1]); it,cells=enum(n)
pat=Counter(); agree=Counter(); boundary_isl=Counter()
for doms in it:
    I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
    spin,free,ring,_=analyse(doms,cells)
    sc={I[k]['S'] and short_cell(*doms[k]):k for k in range(len(doms))}
    ref=v[sc[ring[0]]]
    pat[tuple(v[sc[c]]*ref for c in ring)]+=1
    # compare with old sigma (ring-fixed)
    s0=spin[ring[0]]
    for c,s in spin.items(): agree[(s*s0)==(v[sc[c]]*ref)]+=1
    # islands on the outer boundary: owned by exactly one tetromino and lying on region perimeter
print('n',n,'distinct ring vorticity patterns',len(pat),list(pat.items())[:3])
print('agreement with old sigma on ring-linked short arms',dict(agree))
